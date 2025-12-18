# app.py
import os
import json
import tempfile
import subprocess

import streamlit as st


def run_msianalyzer(marker: str,
                    manifest_file,
                    fastq_files,
                    run_tests: bool,
                    threads: int):
    """Run MSIanalyzer in a temporary working directory and return results."""

    # Use a temporary working directory for this run
    with tempfile.TemporaryDirectory() as tmpdir:
        st.write(f"Working directory: `{tmpdir}`")

        # 1) Save manifest JSON
        manifest_path = os.path.join(tmpdir, manifest_file.name)
        with open(manifest_path, "wb") as f:
            f.write(manifest_file.read())

        # 2) Save FASTQs into a subfolder
        fastq_dir = os.path.join(tmpdir, "fastq")
        os.makedirs(fastq_dir, exist_ok=True)

        fastq_paths = []
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

        # 3) Load and update the manifest JSON
        try:
            with open(manifest_path, "r") as f:
                cfg = json.load(f)
        except json.JSONDecodeError as e:
            st.error("Uploaded manifest file is not valid JSON.")
            st.exception(e)
            st.stop()

        # Overwrite fastq_files to use the uploaded FASTQs
        cfg["fastq_files"] = fastq_paths

        updated_manifest_path = os.path.join(tmpdir, "manifest_runtime.json")
        with open(updated_manifest_path, "w") as f:
            json.dump(cfg, f, indent=2)

        st.markdown("**Updated manifest (runtime copy) written to:**")
        st.code(updated_manifest_path)

        # 4) Build command
        st.subheader("Command to be executed")
        cmd = ["msianalyzer", "run-marker", marker, updated_manifest_path]
        if run_tests:
            cmd.append("--run-tests")
        if threads and threads > 1:
            cmd.extend(["--threads", str(threads)])

        st.code(" ".join(cmd), language="bash")

        # 5) Run MSIanalyzer
        try:
            # Add a timeout to avoid hanging forever if something goes wrong
            result = subprocess.run(
                cmd,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour; adjust as needed
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

        # 6) Show logs
        st.subheader("MSIanalyzer output")

        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**stdout**")
            st.text_area("stdout", result.stdout or "(no stdout)", height=200)
        with col2:
            st.markdown("**stderr**")
            st.text_area("stderr", result.stderr or "(no stderr)", height=200)

        if result.returncode != 0:
            st.error(f"MSIanalyzer exited with non-zero code: {result.returncode}")
            # We still show any files that might have been produced below.
        else:
            st.success("MSIanalyzer finished successfully.")

        # 7) Expose all output files for download, preview images
        st.subheader("Generated files")

        for root, dirs, files in os.walk(tmpdir):
            rel_root = os.path.relpath(root, tmpdir)
            for name in files:
                path = os.path.join(root, name)
                relpath = os.path.join(rel_root, name) if rel_root != "." else name

                # Show image preview for PNG/JPEG
                if name.lower().endswith((".png", ".jpg", ".jpeg")):
                    st.image(path, caption=f"Preview: {relpath}")

                # Download button for every file
                with open(path, "rb") as fh:
                    data = fh.read()
                st.download_button(
                    label=f"Download {relpath}",
                    data=data,
                    file_name=relpath,
                    key=path,  # unique key per file
                )


def main():
    # This must be the first Streamlit call
    st.set_page_config(page_title="MSIanalyzer Web App", layout="wide")

    st.title("MSIanalyzer Web Interface")

    st.markdown(
        """
This web app runs **MSIanalyzer** on user-provided input.

**Workflow:**

1. Enter the marker name (e.g., `BAT25`).
2. Upload a JSON marker/manifest file for that marker.
3. Upload one or more FASTQ files (Nanopore reads for the marker).
4. The app will:
   - Save the uploaded files to a temporary working directory,
   - Rewrite the `fastq_files` field in the JSON to point to the uploaded FASTQs,
   - Run `msianalyzer run-marker`,
   - Show logs and make output files available for download.

For simplicity, the contents of `fastq_files` in your uploaded JSON are **ignored** and replaced by the uploaded FASTQs.
"""
    )

    st.header("1. Input parameters")

    marker = st.text_input(
        "Marker name",
        value="BAT25",
        help="This is the marker passed to `msianalyzer run-marker`.",
    )

    run_tests = st.checkbox(
        "Run MSIanalyzer tests (`--run-tests`)",
        value=False,
        help="If checked, `--run-tests` is added to the command.",
    )

    threads = st.number_input(
        "Threads (`--threads`)",
        min_value=1,
        max_value=32,
        value=4,
        step=1,
        help="Number of threads to use for MSIanalyzer (if >1, `--threads` is used).",
    )

    st.header("2. Upload inputs")

    manifest_file = st.file_uploader(
        "Upload MSIanalyzer JSON marker/manifest file",
        type=["json"],
        help="This JSON describes the marker (primers, reference, etc.). "
             "`fastq_files` will be overwritten.",
    )

    fastq_files = st.file_uploader(
        "Upload FASTQ file(s) for this marker",
        type=["fastq", "fq", "fastq.gz", "fq.gz"],
        accept_multiple_files=True,
        help="You can select multiple FASTQ files at once (e.g. multiple samples or replicates).",
    )

    st.header("3. Run MSIanalyzer")

    if st.button("Run analysis"):
        # Basic validation
        if manifest_file is None:
            st.error("Please upload a JSON marker/manifest file.")
        elif not fastq_files:
            st.error("Please upload at least one FASTQ file.")
        elif not marker.strip():
            st.error("Please provide a marker name.")
        else:
            with st.spinner(
                "Running MSIanalyzer... this may take a few minutes depending on input size."
            ):
                run_msianalyzer(
                    marker=marker.strip(),
                    manifest_file=manifest_file,
                    fastq_files=fastq_files,
                    run_tests=run_tests,
                    threads=int(threads),
                )


if __name__ == "__main__":
    # This is executed when Streamlit runs `app.py`
    main()
