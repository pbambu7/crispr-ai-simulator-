import requests

ENSEMBL_HEADERS = {"Content-Type": "application/json"}

def _lookup_gene_id(symbol: str) -> str | None:
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}?content-type=application/json"
    r = requests.get(url, headers=ENSEMBL_HEADERS, timeout=30)
    if r.status_code != 200:
        return None
    data = r.json()
    if not data:
        return None
    return data[0].get("id")

def _fetch_sequence_by_id(ensembl_id: str) -> str | None:
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?type=genomic"
    r = requests.get(url, headers=ENSEMBL_HEADERS, timeout=30)
    if r.status_code != 200:
        return None
    j = r.json()
    return j.get("seq")

def fetch_ensembl_sequence(gene_symbol: str) -> tuple[str|None, str]:
    """
    Returns (sequence, info_text). info_text includes Ensembl ID if resolved.
    """
    try:
        gene_symbol = gene_symbol.strip()
        ensembl_id = _lookup_gene_id(gene_symbol)
        if not ensembl_id:
            return None, "Ensembl: gene not found."
        seq = _fetch_sequence_by_id(ensembl_id)
        if not seq:
            return None, "Ensembl: sequence unavailable."
        return seq.upper(), f"Ensembl ID: {ensembl_id}"
    except Exception as e:
        return None, f"Ensembl error: {e}"

