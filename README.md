# MSIanalyzer Web App

This Streamlit web interface runs MSIanalyzer on user-uploaded FASTQ files.

Link: 

https://msianalyzer-web-app-tzhai.streamlit.app/ 

## Usage

1. Enter marker name.
2. Select options for your analysis.
3. Upload your FASTQ(s).
4. Click “Run analysis”.

Outputs are presented for download.

## Dev

This repo contains:

- `app.py` — the Streamlit app
- `requirements.txt` — Python deps
- `packages.txt` — system deps for minimap2 & samtools
- hg38 reference files for pileup plot. 
