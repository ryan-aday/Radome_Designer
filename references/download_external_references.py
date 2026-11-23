"""Optional: retrieve publicly hosted PDFs locally for your own research.
Run from package root: python references/download_external_references.py
This is not executed by the Streamlit app. Publisher usage terms apply.
"""
from pathlib import Path
from urllib.request import Request, urlopen
import json

base = Path(__file__).resolve().parent
catalog = json.loads((base / "catalog.json").read_text(encoding="utf-8"))
target = base / "external_downloads"
target.mkdir(exist_ok=True)
for ref in catalog["references"]:
    if not ref.get("pdf_url") or ref.get("local_pdf"):
        continue
    url = ref["pdf_url"]
    output = target / f"{ref['id']}_source.pdf"
    try:
        req = Request(url, headers={"User-Agent": "RadomeStudio-ReferenceDownloader/1.0"})
        with urlopen(req, timeout=25) as response:
            data = response.read(40_000_001)
        if len(data) > 40_000_000 or not data.startswith(b"%PDF-"):
            raise ValueError("Response did not look like a PDF of acceptable size")
        output.write_bytes(data)
        print(f"Saved {output} ({len(data)} bytes)")
    except Exception as exc:
        print(f"Could not retrieve {ref['id']} from {url}: {exc}")
print("Check publishers' terms before sharing downloaded copies.")
