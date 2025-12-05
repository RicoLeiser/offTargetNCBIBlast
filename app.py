import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import time
from io import StringIO

KMER_SIZE = 21
DELAY = 1.5  # seconds

st.title("21-mer NCBI BLAST Pipeline (Cloud Compatible)")

# ------------------- File Upload -------------------
uploaded_file = st.file_uploader(
    "Upload a FASTA file (drag-and-drop supported)",
    type=["fasta", "fa", "fna"]
)

# ------------------- Functions -------------------
def read_fasta(file) -> str:
    """Read a FASTA file and return concatenated sequence(s)"""
    seq = ""
    for line in file:
        line = line.decode().strip() if isinstance(line, bytes) else line.strip()
        if line.startswith(">"):
            continue
        seq += line
    return seq

def generate_kmers(sequence: str, k=21):
    """Return list of all kmers of length k"""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def blast_short(seq: str) -> str:
    """Submit a short sequence to NCBI BLAST and return XML"""
    return NCBIWWW.qblast("blastn", "nt", seq, format_type="XML").read()

def extract_species_from_xml(xml_data) -> list:
    """
    Parse BLAST XML and extract species with 100% identity and coverage.
    Handles cases where the title contains 'PREDICTED:'.
    """
    hits = []
    blast_record = NCBIXML.read(xml_data)
    
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            query_cov = (hsp.align_length / blast_record.query_length) * 100
            identity = (hsp.identities / hsp.align_length) * 100
            if identity == 100.0 and query_cov == 100.0:
                title = alignment.title
                parts = title.split()
                
                # Look for 'PREDICTED:' in the title
                if "PREDICTED:" in parts:
                    idx = parts.index("PREDICTED:")
                    if len(parts) > idx + 2:
                        species = f"{parts[idx+1]} {parts[idx+2]}"
                        hits.append(species)
                else:
                    # Fallback: take the first two words after the accession (usually parts[1] and parts[2])
                    if len(parts) >= 3:
                        species = f"{parts[1]} {parts[2]}"
                        hits.append(species)
    return hits
# ------------------- Main Pipeline -------------------
if uploaded_file is not None:
    st.info("Reading uploaded FASTA...")
    sequence = read_fasta(uploaded_file)
    st.success("FASTA loaded!")

    st.info(f"Generating {KMER_SIZE}-mers...")
    kmers = generate_kmers(sequence, KMER_SIZE)
    st.success(f"{len(kmers)} kmers generated!")

    st.info("Starting BLAST for all kmers (may take some time)...")
    progress_bar = st.progress(0)
    status_text = st.empty()
    species_counts = {}

    for idx, kmer in enumerate(kmers):
        status_text.text(f"BLASTing kmer {idx+1}/{len(kmers)}")
        try:
            xml_string = blast_short(kmer)
            xml_data = StringIO(xml_string)
            hits = extract_species_from_xml(xml_data)
            for sp in hits:
                species_counts[sp] = species_counts.get(sp, 0) + 1
        except Exception as e:
            st.warning(f"Error BLASTing kmer {idx+1}: {e}")
        time.sleep(DELAY)  # respect NCBI

        progress_bar.progress((idx + 1) / len(kmers))

    st.success("BLAST complete!")

    if species_counts:
        df = pd.DataFrame(
            list(species_counts.items()),
            columns=["Latin species name", "Count"]
        )
        st.dataframe(df)

        # Excel download
        excel_file = "blast_output.xlsx"
        df.to_excel(excel_file, index=False)
        with open(excel_file, "rb") as f:
            st.download_button(
                label="Download Excel",
                data=f,
                file_name="blast_output.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    else:
        st.warning("No hits with 100% identity and coverage found.")
