# app.py
import os
import json
import tempfile
import shutil
import subprocess
import base64
from pathlib import Path
from typing import Dict, List

import streamlit as st


# -------------------------------------------------------------------
# 1. Marker definitions
# -------------------------------------------------------------------
MARKERS: Dict[str, Dict[str, str]] = {
    "BAT25": {
        "seq1": "TCGCCTCCAAGAATGTAA",
        "seq2": "ACTATGGCTCTAAAATGCTCTGT",
        "motif": "T",
    },
    "BAT26": {
        "seq1": "TGACTACTTTTGACTTCAGCC",
        "seq2": "AACCATTCAACATTTTTAACC",
        "motif": "A",
    },
    "D2S123": {
        "seq1": "AAACAGGATGCCTGCCTTTA",
        "seq2": "GGACTTTCCACCTATGGGAC",
        "motif": "AC",
    },
    "D5S346": {
        "seq1": "ACTCACTCTAGTGATAAATCGGG",
        "seq2": "AGCAGATAAGACAGTATTACTAGTT",
        "motif": "CA",
    },
    "D17S250": {
        "seq1": "GGAAGAATCAAATAGACAAT",
        "seq2": "GCTGGCCATATATATATTTAAACC",
        "motif": "AC",
    },
}


# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
def infer_stub_from_filename(filename: str) -> str:
    """
    Infer a 'stub' ID from a FASTQ filename, e.g.
    BVSBWG_3_500x.fastq.gz -> BVSBWG_3
    """
    base = os.path.basename(filename)

    # Strip common FASTQ extensions
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break

    parts = base.split("_")
    if len(parts) > 1:
        return "_".join(parts[:-1])
    return base


def parse_group_map_text(text: str) -> Dict[str, str]:
    """
    Parse mapping like:

        BVSBWG_3 = Sample1
        BVSBWG_5 = Sample2

    into {"BVSBWG_3": "Sample1", "BVSBWG_5": "Sample2"}.
    """
    mapping: Dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            stub, name = line.split("=", 1)
            stub = stub.strip()
            name = name.strip()
            if stub and name:
                mapping[stub] = name
    return mapping


def build_runtime_config(
    marker: str,
    fastq_paths: List[str],
    min_similarity: float,
    anchor_units: int,
    user_group_map_text: str,
) -> Dict:
    """
    Build the marker JSON config in memory based on:
    - selected marker
    - uploaded FASTQ files
    - numeric parameters
    - optional user-specified group_map
    """
    if marker not in MARKERS:
        raise ValueError(f"Marker '{marker}' is not defined in MARKERS dict.")

    marker_def = MARKERS[marker]

    # Parse optional user-specified mapping text
    user_map = parse_group_map_text(user_group_map_text or "")

    # Build group_map based on stubs from filenames
    group_map: Dict[str, str] = {}
    for fq in fastq_paths:
        stub = infer_stub_from_filename(fq)
        sample_name = user_map.get(stub, stub)
        group_map[stub] = sample_name

    cfg = {
        "markers": {
            marker: {
                "seq1": marker_def["seq1"],
                "seq2": marker_def["seq2"],
                "motif": marker_def["motif"],
                "group_map": group_map,
                "fastq_files": fastq_paths,
            }
        },
        "min_similarity": float(min_similarity),
        "anchor_units": int(anchor_units),
        "fastq_files": fastq_paths,
    }

    return cfg


def init_session_state():
    """Initialize keys in session_state, if needed."""
    if "results" not in st.session_state:
        st.session_state["results"] = None  # will hold logs + file list + fastqs
    if "tmpdir" not in st.session_state:
        st.session_state["tmpdir"] = None


def clear_results():
    """Clear results and (optionally) delete temp directory."""
    tmpdir = st.session_state.get("tmpdir")
    # Optional: shutil.rmtree(tmpdir, ignore_errors=True)
    st.session_state["results"] = None
    st.session_state["tmpdir"] = None


def embed_svg(path: str, caption: str):
    """Display an SVG file inline using HTML."""
    try:
        with open(path, "rb") as f:
            svg_bytes = f.read()
    except OSError:
        return
    b64 = base64.b64encode(svg_bytes).decode("utf-8")
    html = f'<img src="data:image/svg+xml;base64,{b64}" alt="{caption}" />'
    st.markdown(html, unsafe_allow_html=True)


# Local hg38 reference directory, vendored into this repo
EXAMPLES_HG38_DIR = Path(__file__).resolve().parent / "hg38"


def get_reference_for_marker(marker: str) -> str:
    """
    Return the path to the reference FASTA for the given marker,
    by searching in the local hg38/ directory.
    """
    if not EXAMPLES_HG38_DIR.is_dir():
        raise FileNotFoundError(
            f"Reference directory not found: {EXAMPLES_HG38_DIR}. "
            "Make sure the hg38 folder is included in your repo."
        )

    marker_lower = marker.lower()

    for filename in os.listdir(EXAMPLES_HG38_DIR):
        if filename.lower().endswith((".fa", ".fasta")):
            if marker_lower in filename.lower():
                return str(EXAMPLES_HG38_DIR / filename)

    raise FileNotFoundError(
        f"No reference FASTA found in {EXAMPLES_HG38_DIR} for marker '{marker}'."
    )


# -------------------------------------------------------------------
# Core runner (does the heavy work once, stores in session_state)
# -------------------------------------------------------------------
def run_msianalyzer_once(
    marker: str,
    fastq_files,
    min_similarity: float,
    anchor_units: int,
    run_tests: bool,
    skip_variant_summary: bool,
    threads: int,
    group_map_text: str,
):
    """Run MSIanalyzer in a new temporary working directory and store results in session_state."""

    # New temp dir for this run
    tmpdir = tempfile.mkdtemp(prefix="msianalyzer_")
    st.session_state["tmpdir"] = tmpdir

    st.write(f"Working directory: `{tmpdir}`")

    # 1) Save FASTQs into a subfolder
    fastq_dir = os.path.join(tmpdir, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)

    fastq_paths: List[str] = []
    for fobj in fastq_files:
        fq_path = os.path.join(fastq_dir, fobj.name)
        with open(fq_path, "wb") as fh:
            fh.write(fobj.read())
        fastq_paths.append(fq_path)

    if not fastq_paths:
        st.error("No FASTQ files were saved. Please re-upload your FASTQ files.")
        st.stop()

    st.markdown("**FASTQ files saved:**")
    for p in fastq_paths:
        st.write(os.path.basename(p))

    # 2) Build runtime config in memory
    try:
        cfg = build_runtime_config(
            marker=marker,
            fastq_paths=fastq_paths,
            min_similarity=min_similarity,
            anchor_units=anchor_units,
            user_group_map_text=group_map_text,
        )
    except Exception as e:
        st.error("Failed to build marker configuration.")
        st.exception(e)
        st.stop()

    # 3) Write config to JSON in temp dir
    manifest_path = os.path.join(tmpdir, "manifest_runtime.json")
    with open(manifest_path, "w") as f:
        json.dump(cfg, f, indent=2)

    st.markdown("**Runtime manifest written to:**")
    st.code(manifest_path)

    # 4) Build command
    st.subheader("Command to be executed")
    cmd = ["msianalyzer", "run-marker", marker, manifest_path]
    if run_tests:
        cmd.append("--run-tests")
    if skip_variant_summary:
        cmd.append("--skip-variant-summary")
    if threads and threads > 1:
        cmd.extend(["--threads", str(threads)])

    st.code(" ".join(cmd), language="bash")

    # 5) Run MSIanalyzer
    try:
        result = subprocess.run(
            cmd,
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour; adjust if needed
        )
    except FileNotFoundError as e:
        st.error(
            "Failed to run `msianalyzer`. Make sure it is installed and in PATH "
            "(via `requirements.txt`)."
        )
        st.exception(e)
        st.stop()
    except subprocess.TimeoutExpired as e:
        st.error("MSIanalyzer run timed out.")
        st.exception(e)
        st.stop()

    # 6) Collect logs and files into session_state
    files_info = []
    for root, dirs, files in os.walk(tmpdir):
        rel_root = os.path.relpath(root, tmpdir)
        for name in files:
            path = os.path.join(root, name)
            relpath = os.path.join(rel_root, name) if rel_root != "." else name
            files_info.append(
                {
                    "path": path,
                    "relpath": relpath,
                    "name": name,
                }
            )

    st.session_state["results"] = {
        "stdout": result.stdout or "",
        "stderr": result.stderr or "",
        "returncode": result.returncode,
        "files": files_info,
        "fastq_paths": fastq_paths,  # store full paths for pileup
        "marker": marker,
    }


def refresh_results_files():
    """Re-scan tmpdir and refresh the file list in st.session_state['results']."""
    tmpdir = st.session_state.get("tmpdir")
    res = st.session_state.get("results")
    if not tmpdir or not res:
        return
    if not os.path.isdir(tmpdir):
        return

    files_info = []
    for root, dirs, files in os.walk(tmpdir):
        rel_root = os.path.relpath(root, tmpdir)
        for name in files:
            path = os.path.join(root, name)
            relpath = os.path.join(rel_root, name) if rel_root != "." else name
            files_info.append(
                {
                    "path": path,
                    "relpath": relpath,
                    "name": name,
                }
            )
    res["files"] = files_info
    st.session_state["results"] = res


def run_pileup_and_show(marker: str, fastq_path: str):
    """
    Run the `msianalyzer pileup` command for the given FASTQ and reference.
    Uses the same tmpdir as the main analysis and updates results.
    """

    tmpdir = st.session_state.get("tmpdir")
    if not tmpdir or not os.path.isdir(tmpdir):
        st.error("No working directory found. Please run the main analysis first.")
        return

    try:
        ref_path = get_reference_for_marker(marker)
    except FileNotFoundError as e:
        st.error(str(e))
        return

    st.write(f"Pileup working directory: `{tmpdir}`")

    if not os.path.isfile(fastq_path):
        st.error(f"FASTQ file not found on disk for pileup: {fastq_path}")
        return

    # Build the pileup command
    cmd = ["msianalyzer", "pileup", fastq_path, ref_path]
    st.code(" ".join(cmd), language="bash")

    # Execute the pileup
    try:
        result = subprocess.run(
            cmd,
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=1800,
        )
    except FileNotFoundError:
        st.error("Failed to run `msianalyzer pileup`. Make sure it’s installed in PATH.")
        return
    except subprocess.TimeoutExpired:
        st.error("Pileup command timed out.")
        return

    st.subheader("Pileup command output")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**stdout (pileup)**")
        st.text_area("pileup_stdout", result.stdout or "(no stdout)", height=200)
    with col2:
        st.markdown("**stderr (pileup)**")
        st.text_area("pileup_stderr", result.stderr or "(no stderr)", height=200)

    if result.returncode != 0:
        st.error(f"`msianalyzer pileup` exited with code {result.returncode}")
        return

    # Update files list in results so new plots appear in the general "Generated files" section
    refresh_results_files()
    st.success("Pileup completed. New files have been added to the generated files list below.")

    # Optionally, also show any image files created by pileup directly
    plot_files = [
        finfo["path"]
        for finfo in st.session_state["results"]["files"]
        if finfo["name"].lower().endswith((".png", ".jpg", ".jpeg", ".svg"))
    ]

    if plot_files:
        st.subheader("Pileup Plot Results")
        for filepath in plot_files:
            name = os.path.basename(filepath)
            if name.lower().endswith(".svg"):
                st.markdown(f"**SVG: {name}**")
                embed_svg(filepath, caption=name)
            else:
                st.image(filepath, caption=name)
            with open(filepath, "rb") as fh:
                data = fh.read()
            st.download_button(
                label=f"Download {name}",
                data=data,
                file_name=name,
                key=f"pileup_{name}",
            )


# -------------------------------------------------------------------
# UI rendering of stored results (no heavy work)
# -------------------------------------------------------------------
def render_results():
    """Render logs and generated files from st.session_state['results']."""
    res = st.session_state.get("results")
    if not res:
        return

    st.subheader("MSIanalyzer output")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**stdout**")
        st.text_area("stdout", res["stdout"] or "(no stdout)", height=200)
    with col2:
        st.markdown("**stderr**")
        st.text_area("stderr", res["stderr"] or "(no stderr)", height=200)

    if res["returncode"] != 0:
        st.error(f"MSIanalyzer exited with non-zero code: {res['returncode']}")
    else:
        st.success("MSIanalyzer finished successfully.")

    st.subheader("Generated files")

    for finfo in res["files"]:
        path = finfo["path"]
        relpath = finfo["relpath"]
        name = finfo["name"]

        # Preview images: PNG/JPG
        if name.lower().endswith((".png", ".jpg", ".jpeg")):
            st.image(path, caption=f"Preview: {relpath}")

        # Preview SVG
        if name.lower().endswith(".svg"):
            st.markdown(f"**SVG preview: {relpath}**")
            embed_svg(path, caption=relpath)

        # Download button for every file
        try:
            with open(path, "rb") as fh:
                data = fh.read()
        except OSError:
            continue

        st.download_button(
            label=f"Download {relpath}",
            data=data,
            file_name=relpath,
            key=path,  # must be unique per file
        )


# -------------------------------------------------------------------
# Main UI
# -------------------------------------------------------------------
def main():
    st.set_page_config(page_title="MSIanalyzer Web App", layout="wide")
    init_session_state()

    st.title("MSIanalyzer Web Interface")

    st.markdown(
        """
        This web app provides a browser-based interface to the
        [**MSIanalyzer**](https://github.com/NagelLabHub/MSIanalyzer) command-line tool.

        **Workflow:**

        1. Select the Bethesda panel marker (BAT25, BAT26, D2S123, D5S346, D17S250).
        2. Optionally adjust `min_similarity` and `anchor_units`.
        3. Upload one or more FASTQ files (Nanopore reads for that marker).
        4. Optionally provide a mapping from FASTQ-derived IDs to sample/group names.
        5. Click **Run analysis**.
        6. Optionally generate a pileup plot for one of the FASTQs.
        7. You can download any output files without losing the others.
        8. Click **Clear results** if you want to reset and start again.

        Questions or issues: open an issue on the GitHub repo.
        """
    )

    # Top-level controls
    st.header("1. Marker and parameters")

    marker = st.selectbox(
        "Marker",
        options=list(MARKERS.keys()),
        index=0,
        help="Marker name passed to `msianalyzer run-marker`.",
    )

    min_similarity = st.number_input(
        "Minimum similarity (`min_similarity`)",
        min_value=0.0,
        max_value=1.0,
        value=0.85,
        step=0.01,
        help="Similarity threshold for anchor-extension matching.",
    )

    anchor_units = st.number_input(
        "Anchor units (`anchor_units`)",
        min_value=1,
        max_value=20,
        value=3,
        step=1,
        help="Number of repeat units used as the anchor.",
    )

    run_tests = st.checkbox(
        "Run MSIanalyzer tests (`--run-tests`)",
        value=False,
        help="If checked, `--run-tests` is added to the command.",
    )

    skip_variant_summary = st.checkbox(
        "Skip variant summary (`--skip-variant-summary`)",
        value=False,
        help="If checked, `--skip-variant-summary` is added to the command.",
    )

    threads = st.number_input(
        "Threads (`--threads`)",
        min_value=1,
        max_value=32,
        value=4,
        step=1,
        help="Number of threads to use for MSIanalyzer (if >1, `--threads` is used).",
    )

    st.header("2. Upload FASTQ files")

    fastq_files = st.file_uploader(
        "Upload FASTQ file(s) for this marker",
        type=["fastq", "fq", "fastq.gz", "fq.gz"],
        accept_multiple_files=True,
        help="You can select multiple FASTQ files at once (e.g. multiple samples or replicates).",
    )

    st.header("3. Optional: specify sample names")

    st.markdown(
        """
        By default, sample IDs are inferred from FASTQ filenames by removing the extension
        and dropping the last underscore-delimited token (e.g., `BVSBWG_3_500x.fastq` → stub `BVSBWG_3`),
        and the stub itself is used as the sample name.

        You can override this by providing mappings like:

        BVSBWG_3 = Sample1

        One mapping per line, `stub = SampleName`.
        """
    )

    group_map_text = st.text_area(
        "Optional stub → sample name mapping",
        value="",
        height=120,
    )

    st.header("4. Run / Clear")

    cols = st.columns(2)
    with cols[0]:
        run_clicked = st.button("Run analysis")
    with cols[1]:
        clear_clicked = st.button("Clear results")

    if clear_clicked:
        clear_results()
        st.info("Results cleared. You can upload new files and run again.")

    if run_clicked:
        if not fastq_files:
            st.error("Please upload at least one FASTQ file.")
        else:
            with st.spinner(
                "Running MSIanalyzer... this may take a few minutes depending on input size."
            ):
                # Running will overwrite any previous results and tmpdir
                run_msianalyzer_once(
                    marker=marker,
                    fastq_files=fastq_files,
                    min_similarity=float(min_similarity),
                    anchor_units=int(anchor_units),
                    run_tests=run_tests,
                    skip_variant_summary=skip_variant_summary,
                    threads=int(threads),
                    group_map_text=group_map_text,
                )

    st.header("5. Pileup Plot")

    res = st.session_state.get("results")
    if not res or not res.get("fastq_paths"):
        st.info("Run the main analysis first to enable pileup plotting.")
    else:
        fastq_paths = res["fastq_paths"]
        fastq_basenames = [os.path.basename(p) for p in fastq_paths]

        pileup_fastq = st.selectbox(
            "Select a FASTQ file to generate pileup",
            options=fastq_basenames,
        )

        if st.button("Generate pileup plot"):
            try:
                idx = fastq_basenames.index(pileup_fastq)
                fq_full_path = fastq_paths[idx]
            except (ValueError, IndexError):
                st.error("Selected FASTQ not found. Please rerun the analysis.")
            else:
                with st.spinner("Running `msianalyzer pileup`..."):
                    # Use the marker stored with the results (in case user changed marker in UI)
                    marker_for_pileup = res.get("marker", marker)
                    run_pileup_and_show(marker=marker_for_pileup, fastq_path=fq_full_path)

    # Always render results if present (survive download button clicks)
    render_results()


if __name__ == "__main__":
    main()
