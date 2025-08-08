
import streamlit as st
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import requests
from io import StringIO, BytesIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter

st.set_page_config(page_title="🧬 CRISPR + AI Simulator", layout="centered")

st.markdown("<h1 style='text-align: center;'>CRISPR + AI Simulator</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: grey;'>by Patience Bambu | v5 with Experiment Simulation + PDF</p>", unsafe_allow_html=True)

# Sidebar
st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/e/e3/CRISPR-Cas9.svg", width=120)
st.sidebar.title("Options")
mode = st.sidebar.radio("Select Input Mode", ["🧬 Paste DNA", "📁 Upload FASTA", "🔍 Fetch Real Gene (Ensembl)"])
vector = st.sidebar.selectbox("Delivery Vector", ["Lipid Nanoparticles", "AAV", "Electroporation"])
simulate_immune = st.sidebar.checkbox("Simulate Immune Risk", True)

# Simulation model for immune & delivery
def simulate_outcomes(vector_choice, gRNA_count):
    if vector_choice == "Lipid Nanoparticles":
        success = 90
        immune_risk = 10
    elif vector_choice == "AAV":
        success = 75
        immune_risk = 25
    else:
        success = 65
        immune_risk = 20
    success -= gRNA_count
    immune_risk += gRNA_count
    return max(success, 0), min(immune_risk, 100)

# gRNA predictors
def predict_gRNA_binding_sites(dna_seq, pam="NGG"):
    candidates = []
    for i in range(len(dna_seq) - 23):
        gRNA = dna_seq[i:i+20]
        pam_seq = dna_seq[i+20:i+23]
        if pam_seq.endswith("GG"):
            candidates.append((i, gRNA, pam_seq))
    return candidates

def score_gRNA(gRNA):
    gc = (gRNA.count("G") + gRNA.count("C")) / len(gRNA)
    penalty = abs(0.5 - gc)
    score = round(100 - penalty * 100, 2)
    return score

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

# Generate PDF Report
def generate_pdf_report(df, vector, immune, success, gene_name):
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    content = [Paragraph("CRISPR + AI Gene Editing Report", styles["Title"]),
               Paragraph(f"Gene: {gene_name}", styles["Heading2"]),
               Spacer(1, 12)]

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

# Process Sequence
def display_results(dna_seq, gene_name="N/A"):
    st.subheader("2. Candidate gRNAs + AI Score")
    gRNAs = predict_gRNA_binding_sites(dna_seq)
    if not gRNAs:
        st.warning("No valid gRNA candidates found.")
        return

    data = []
    for i, gRNA, pam in gRNAs[:10]:
        score = score_gRNA(gRNA)
        data.append({"Position": i, "Guide RNA": gRNA, "PAM": pam, "GC Score": score})

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

# Input Modes
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
