# MSIanalyzer Web App

This Streamlit web interface runs MSIanalyzer on user-uploaded JSON + FASTQ files.

Deploys to Streamlit Community Cloud.

## Usage

1. Enter marker name.
2. Upload the marker JSON.
3. Upload your FASTQ(s).
4. Click “Run analysis”.

Outputs are presented for download.

## Dev

This repo contains:

- `app.py` — the Streamlit app
- `requirements.txt` — Python deps
- `packages.txt` — system deps for minimap2 & samtools
