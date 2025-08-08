import streamlit as st
from Bio import SeqIO
import numpy as np
import pandas as pd
import requests
from io import StringIO, BytesIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, Image as RLImage
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from PIL import Image

st.set_page_config(page_title="🧬 CRISPR + AI Simulator", layout="centered")

st.markdown("<h1 style='text-align: center;'>CRISPR + AI Simulator</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: grey;'>by Patience Bambu | v6 with Gene Fetch, Simulation, Visualization, PDF Export</p>", unsafe_allow_html=True)

# Sidebar with local logo
try:
    sidebar_logo = Image.open("assets/logo.png")
    st.sidebar.image(sidebar_logo, width=120)
except:
    st.sidebar.warning("Logo not found in assets/")

st.sidebar.title("Options")
mode = st.sidebar.radio("Select Input Mode", ["🧬 Paste DNA", "📁 Upload FASTA", "🔍 Fetch Real Gene (Ensembl)"])
vector = st.sidebar.selectbox("Delivery Vector", ["Lipid Nanoparticles", "AAV", "Electroporation"])
simulate_immune = st.sidebar.checkbox("Simulate Immune Risk", True)

# Simulation model
def simulate_outcomes(vector_choice, gRNA_count):
    success, immune_risk = 90, 10
    if vector_choice == "AAV":
        success, immune_risk = 75, 25
    elif vector_choice == "Electroporation":
        success, immune_risk = 65, 20
    success -= gRNA_count
    immune_risk += gRNA_count
    return max(success, 0), min(immune_risk, 100)

# gRNA prediction
def predict_gRNA_binding_sites(dna_seq, pam="NGG"):
    candidates = []
    for i in range(len(dna_seq) - 23):
        gRNA = dna_seq[i:i+20]
        pam_seq = dna_seq[i+20:i+23]
        if pam_seq.endswith("GG"):
            candidates.append((i, gRNA, pam_seq))
    return candidates

# AI-based gRNA scoring
def score_gRNA(gRNA):
    gc = (gRNA.count("G") + gRNA.count("C")) / len(gRNA)
    penalty = abs(0.5 - gc)
    return round(100 - penalty * 100, 2)

# Fetch gene sequence from Ensembl
def fetch_gene_sequence_ensembl(gene_name):
    try:
        url_lookup = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
        res = requests.get(url_lookup, headers={"Content-Type": "application/json"})
        if res.status_code != 200 or not res.json():
            return None, "Gene not found."
        gene_id = res.json()[0]['id']
        url_seq = f"https://rest.ensembl.org/sequence/id/{gene_id}?type=genomic"
        res_seq = requests.get(url_seq, headers={"Content-Type": "application/json"})
        if res_seq.status_code == 200:
            return res_seq.json().get("seq"), f"Gene ID: {gene_id}"
        else:
            return None, "Sequence not available."
    except Exception as e:
        return None, str(e)

# Generate PDF report
def generate_pdf_report(df, vector, immune, success, gene_name):
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    content = []
    try:
        logo_path = "assets/logo.png"
        content.append(RLImage(logo_path, width=1.5*inch, height=1.5*inch))
    except:
        pass
    content += [
        Paragraph("CRISPR + AI Gene Editing Report", styles["Title"]),
        Paragraph(f"Gene: {gene_name}", styles["Heading2"]),
        Spacer(1, 12)
    ]
    data = [["Position", "gRNA", "PAM", "Score"]]
    for _, row in df.iterrows():
        data.append([str(row["Position"]), row["Guide RNA"], row["PAM"], str(row["GC Score"])])
    content.append(Table(data))
    content.append(Spacer(1, 12))
    content.append(Paragraph(f"Vector: {vector}", styles["Normal"]))
    content.append(Paragraph(f"Simulated Success: {success}%", styles["Normal"]))
    if simulate_immune:
        content.append(Paragraph(f"Immune Risk: {immune}%", styles["Normal"]))
    doc.build(content)
    buffer.seek(0)
    return buffer

# Display results
def display_results(dna_seq, gene_name="N/A"):
    st.subheader("2. Candidate gRNAs + AI Score")
    gRNAs = predict_gRNA_binding_sites(dna_seq)
    if not gRNAs:
        st.warning("No valid gRNA candidates found.")
        return
    data = []
    for i, gRNA, pam in gRNAs[:10]:
        data.append({"Position": i, "Guide RNA": gRNA, "PAM": pam, "GC Score": score_gRNA(gRNA)})
    df = pd.DataFrame(data)
    st.dataframe(df)

    st.subheader("3. Delivery + Immune Simulation")
    success, immune_risk = simulate_outcomes(vector, len(df))
    st.success(f"Delivery Success: {success}%")
    if simulate_immune:
        st.error(f"Immune Risk: {immune_risk}%")

    st.subheader("4. Simulated Gene Edit (Preview)")
    cut_index = gRNAs[0][0] + 10
    edited_seq = dna_seq[:cut_index] + "---CUT---" + dna_seq[cut_index+3:]
    st.code(edited_seq)

    st.subheader("5. Export Full Report")
    pdf = generate_pdf_report(df, vector, immune_risk, success, gene_name)
    st.download_button("📄 Download PDF Report", pdf, "CRISPR_AI_Report.pdf", mime="application/pdf")

# Main logic
if mode == "🧬 Paste DNA":
    st.subheader("1. Paste DNA Sequence")
    user_input = st.text_area("Enter a DNA sequence (A, T, G, C only):", height=200)
    if user_input:
        display_results(user_input.upper())

elif mode == "📁 Upload FASTA":
    st.subheader("1. Upload FASTA File")
    fasta_file = st.file_uploader("Choose a .fasta file", type=["fasta", "fa"])
    if fasta_file is not None:
        try:
            fasta_string = fasta_file.read().decode("utf-8")
            fasta_io = StringIO(fasta_string)
            record = SeqIO.read(fasta_io, "fasta")
            st.success(f"Loaded: " + record.id)
            display_results(str(record.seq).upper(), gene_name=record.id)
        except Exception as e:
            st.error(f"Error reading FASTA file: {e}")

elif mode == "🔍 Fetch Real Gene (Ensembl)":
    st.subheader("1. Search by Gene Name")
    gene_name = st.text_input("Enter gene name (e.g. BRCA1):")
    if gene_name:
        with st.spinner("Fetching gene data..."):
            gene_seq, info = fetch_gene_sequence_ensembl(gene_name.strip())
        if gene_seq:
            st.success(info)
            st.text_area("Sequence Preview", gene_seq[:2000] + "..." if len(gene_seq) > 2000 else gene_seq, height=150)
            display_results(gene_seq.upper(), gene_name=gene_name)
        else:
            st.error(f"Failed: {info}")
