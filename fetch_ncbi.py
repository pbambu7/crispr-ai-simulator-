import requests

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def fetch_ncbi_sequence(search_term: str) -> tuple[str|None, str]:
    """
    Search NCBI nucleotide by free-text and fetch the first FASTA.
    Returns (sequence, info_text).
    """
    try:
        params = {
            "db": "nucleotide",
            "retmode": "json",
            "term": search_term
        }
        s = requests.get(f"{NCBI_EUTILS}/esearch.fcgi", params=params, timeout=30)
        s.raise_for_status()
        ids = s.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None, "NCBI: no hits. Try a different term."

        nuccore_id = ids[0]
        f = requests.get(
            f"{NCBI_EUTILS}/efetch.fcgi",
            params={"db":"nucleotide","id":nuccore_id,"rettype":"fasta","retmode":"text"},
            timeout=60
        )
        f.raise_for_status()
        raw = f.text.strip().splitlines()
        header = raw[0] if raw else ""
        seq = "".join([line.strip() for line in raw[1:]]).upper()
        if not seq:
            return None, "NCBI: empty sequence."
        return seq, f"NCBI nuccore ID: {nuccore_id} | {header}"
    except Exception as e:
        return None, f"NCBI error: {e}"

