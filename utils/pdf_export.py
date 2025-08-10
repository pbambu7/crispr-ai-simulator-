from __future__ import annotations
from io import BytesIO
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from datetime import datetime

def build_pdf_report(
    df,
    vector: str,
    immune_risk: int,
    success: int,
    gene_title: str,
    cut_preview: str,
    user_name: str | None = None,
    user_age: str | None = None,
    logo_path: str | None = None,
) -> BytesIO:
    buf = BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter, topMargin=36, bottomMargin=36, leftMargin=36, rightMargin=36)
    styles = getSampleStyleSheet()
    items = []

    if logo_path:
        try:
            items.append(Image(logo_path, width=1.2*inch, height=1.2*inch))
        except Exception:
            pass

    heading = "CRISPR + AI Simulator Report"
    items.append(Paragraph(heading, styles["Title"]))
    items.append(Spacer(1, 6))
    meta = f"Generated: {datetime.now().strftime('%Y-%m-%d %I:%M %p')}"
    items.append(Paragraph(meta, styles["Normal"]))
    if user_name or user_age:
        items.append(Paragraph(f"User: {user_name or '—'}  |  Age: {user_age or '—'}", styles["Normal"]))
    items.append(Spacer(1, 12))

    items.append(Paragraph(gene_title, styles["Heading2"]))
    items.append(Paragraph(f"Delivery Vector: {vector}", styles["Normal"]))
    items.append(Paragraph(f"Delivery Success (sim): {success}%", styles["Normal"]))
    items.append(Paragraph(f"Immune Risk (sim): {immune_risk}%", styles["Normal"]))
    items.append(Spacer(1, 12))

    # table of guides
    data = [list(df.columns)] + [list(map(str, r)) for r in df.values[:10]]
    items.append(Table(data, repeatRows=1))
    items.append(Spacer(1, 12))
    items.append(Paragraph("Simulated Cut Preview:", styles["Heading3"]))
    items.append(Paragraph(cut_preview[:2500], styles["Code"]))
    doc.build(items)
    buf.seek(0)
    return buf
