import streamlit as st
from Bio import SeqIO
import numpy as np
import pandas as pd
import requests
from io import StringIO, BytesIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
import matplotlib.pyplot as plt

st.set_page_config(page_title="🧬 CRISPR + AI Simulator", layout="centered")
st.markdown("<h1 style='text-align: center;'>CRISPR + AI Simulator</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: grey;'>by Patience Bambu | v7: Charts + CpG + Codon Bias + Off-target</p>", unsafe_allow_html=True)

st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/e/e3/CRISPR-Cas9.svg", width=120)
st.sidebar.title("Options")
mode = st.sidebar.radio("Select Input Mode", ["🧬 Paste DNA", "📁 Upload FASTA", "🔍 Fetch Real Gene (Ensembl)"])
vector = st.sidebar.selectbox("Delivery Vector", ["Lipid Nanoparticles", "AAV", "Electroporation"])
simulate_immune = st.sidebar.checkbox("Simulate Immune Risk", True)

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
    cpg_density = gRNA.count("CG") / len(gRNA)
    codon_bias = np.random.uniform(0.8, 1.2)  # Simulated codon bias factor
    base_score = 100 - abs(0.5 - gc) * 100
    adjusted_score = base_score * (1 - 0.2 * cpg_density) * codon_bias
    return round(adjusted_score, 2), round(cpg_density, 2), round(codon_bias, 2)

def simulate_outcomes(vector_choice, gRNA_count, avg_score):
    base_success = 90 if vector_choice == "Lipid Nanoparticles" else 75 if vector_choice == "AAV" else 65
    immune_risk = 10 if vector_choice == "Lipid Nanoparticles" else 25 if vector_choice == "AAV" else 20
    penalty = gRNA_count * (1 - avg_score/100)
    success = max(base_success - penalty, 10)
    immune_risk = min(immune_risk + penalty/2, 90)
    return round(success), round(immune_risk)

def plot_scores(df):
    fig, ax = plt.subplots()
    ax.plot(df["Position"], df["GC Score"], label="GC Score")
    ax.plot(df["Position"], df["CpG"], label="CpG Density")
    ax.plot(df["Position"], df["Codon Bias"], label="Codon Bias")
    ax.set_xlabel("Position")
    ax.set_ylabel("Score")
    ax.set_title("gRNA Quality Metrics")
    ax.legend()
    st.pyplot(fig)

def display_results(dna_seq, gene_name="N/A"):
    st.subheader("2. Candidate gRNAs + AI Score")
    gRNAs = predict_gRNA_binding_sites(dna_seq)
    if not gRNAs:
        st.warning("No valid gRNA candidates found.")
        return
    data = []
    for i, gRNA, pam in gRNAs[:10]:
        score, cpg, bias = score_gRNA(gRNA)
        data.append({"Position": i, "Guide RNA": gRNA, "PAM": pam, "GC Score": score, "CpG": cpg, "Codon Bias": bias})
    df = pd.DataFrame(data)
    st.dataframe(df)
    plot_scores(df)

    st.subheader("3. Delivery + Immune Simulation")
    avg_score = df["GC Score"].mean()
    success, immune_risk = simulate_outcomes(vector, len(df), avg_score)
    st.success(f"Delivery Success: {success}%")
    if simulate_immune:
        st.error(f"Immune Risk: {immune_risk}%")

    st.subheader("4. Simulated Gene Edit (Preview)")
    cut_index = gRNAs[0][0] + 10
    edited_seq = dna_seq[:cut_index] + "---CUT---" + dna_seq[cut_index+3:]
    st.code(edited_seq)

    st.subheader("5. Export Full Report")
    st.markdown("PDF export includes gene name, vector, score table. Off-target % will be in v8.")

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
            try:
                url_lookup = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
                res = requests.get(url_lookup)
                gene_id = res.json()[0]['id']
                url_seq = f"https://rest.ensembl.org/sequence/id/{gene_id}?type=genomic"
                seq_res = requests.get(url_seq)
                gene_seq = seq_res.json().get("seq")
                st.success(f"Gene ID: {gene_id}")
                display_results(gene_seq.upper(), gene_name=gene_name)
            except:
                st.error("Gene not found or API error.")

