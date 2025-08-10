import os
from io import StringIO
import pandas as pd
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt

from utils import (
    fetch_ensembl_sequence, fetch_ncbi_sequence,
    find_spcas9_sites, score_guides, simulated_cut_preview,
    compute_delivery_success, compute_immune_risk,
    build_pdf_report
)

# ------------- Page & Style -------------
st.set_page_config(page_title="🧬 CRISPR + AI Simulator", layout="wide")
if os.path.exists("static/style.css"):
    with open("static/style.css", "r") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

st.markdown("# CRISPR + AI Simulator")
st.caption("by Patience Bambu | v7: Charts + CpG + Codon Bias + Off-target + PDF")

# ------------- Sidebar -------------
st.sidebar.header("Options")
mode = st.sidebar.radio("Select Input Mode", ["🧬 Paste DNA", "📁 Upload FASTA", "🧫 Ensembl (human)", "🧪 NCBI (any)"])
vector = st.sidebar.selectbox("Delivery Vector", ["Lipid Nanoparticles", "AAV", "Electroporation"])
simulate_immune = st.sidebar.checkbox("Simulate Immune Risk", True)

st.sidebar.markdown("---")
st.sidebar.subheader("Report details (optional)")
user_name = st.sidebar.text_input("User name")
user_age = st.sidebar.text_input("User age")
logo_path = "static/logo.png" if os.path.exists("static/logo.png") else None

# ------------- Input + Fetch -------------
seq = None
gene_info = "Sequence"

if mode == "🧬 Paste DNA":
    st.subheader("1. Paste DNA")
    _txt = st.text_area("DNA sequence (A/T/G/C):", height=160)
    if _txt:
        seq = "".join([c for c in _txt.upper() if c in "ATGC"])
        if not seq:
            st.error("No valid bases found.")
elif mode == "📁 Upload FASTA":
    st.subheader("1. Upload FASTA")
    f = st.file_uploader("Choose a .fasta or .fa", type=["fasta", "fa"])
    if f:
        try:
            lines = f.read().decode("utf-8").splitlines()
            header = lines[0] if lines and lines[0].startswith(">") else ">FASTA"
            gene_info = header
            seq = "".join(line.strip().upper() for line in lines if not line.startswith(">"))
        except Exception as e:
            st.error(f"FASTA read error: {e}")
elif mode == "🧫 Ensembl (human)":
    st.subheader("1. Search (Ensembl gene symbol)")
    symbol = st.text_input("Enter gene symbol (e.g., BRCA1)")
    if symbol:
        with st.spinner("Fetching from Ensembl..."):
            seq, info = fetch_ensembl_sequence(symbol.strip())
        if seq:
            st.success(info)
            gene_info = info
        else:
            st.error(info or "Ensembl error.")
elif mode == "🧪 NCBI (any)":
    st.subheader("1. Search (NCBI)")
    query = st.text_input("Free text (e.g., 'HSV-1 complete genome', 'UL30', 'HIV pol')")
    if query:
        with st.spinner("Searching NCBI..."):
            seq, info = fetch_ncbi_sequence(query.strip())
        if seq:
            st.success(info)
            gene_info = info
        else:
            st.error(info or "NCBI error.")

if not seq:
    st.stop()

# ------------- Preview -------------
st.subheader("Sequence Preview")
st.text_area("First 1200 bp (truncated for display)", value=seq[:1200], height=180)

# ------------- Guide Prediction & Scoring -------------
st.subheader("2. Candidate gRNAs + Scores")
rows = score_guides(seq, max_rows=10)
if not rows:
    st.warning("No SpCas9 NGG sites found in the first window.")
    st.stop()

df = pd.DataFrame(rows)
st.dataframe(df, use_container_width=True)

# ------------- Charts (Matplotlib) -------------
st.subheader("3. Visualization")
fig, ax = plt.subplots()
x = np.arange(len(df))
ax.plot(x, df["GC Score"], marker="o", label="GC%")
ax.plot(x, df["CpG"]*100, marker="o", label="CpG (×100)")
ax.plot(x, df["Off-target (↑=better)"], marker="o", label="Off-target")
ax.set_xticks(x)
ax.set_xticklabels(df["Position"].astype(str), rotation=0)
ax.set_xlabel("Guide position")
ax.set_ylabel("Score")
ax.legend()
st.pyplot(fig)

# ------------- Simulation -------------
avg_gc = float(df["GC Score"].mean())
avg_cpg = float(df["CpG"].mean())
success = compute_delivery_success(vector, avg_gc, avg_cpg)
immune = compute_immune_risk(simulate_immune, avg_cpg, vector)

st.subheader("4. Delivery + Immune Simulation")
st.success(f"Delivery Success (sim): **{success}%**")
if simulate_immune:
    st.error(f"Immune Risk (sim): **{immune}%**")

# ------------- Simulated Cut Preview -------------
st.subheader("5. Simulated Gene Edit (Preview)")
first_pos = int(df.iloc[0]["Position"]) if len(df) else None
cut_txt = simulated_cut_preview(seq, first_pos)
st.code(cut_txt[:2000])

# ------------- PDF Export -------------
st.subheader("6. Export Full Report")
title = f"{gene_info}"
pdf_buffer = build_pdf_report(
    df=df,
    vector=vector,
    immune_risk=immune,
    success=success,
    gene_title=title,
    cut_preview=cut_txt,
    user_name=user_name,
    user_age=user_age,
    logo_path=logo_path
)
st.download_button(
    "📄 Download PDF Report",
    data=pdf_buffer,
    file_name="CRISPR_AI_Report.pdf",
    mime="application/pdf"
)

st.caption("Educational use only. Not for clinical decisions.")
