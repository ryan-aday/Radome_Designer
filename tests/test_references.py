"""Source-reference integrity tests; run with pytest -q from project root."""
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]


def test_all_original_sources_listed_and_have_valid_urls():
    refs=json.loads((BASE/'references/catalog.json').read_text(encoding='utf-8'))['references']
    assert [r['id'] for r in refs] == [f'R{i:02d}' for i in range(1,8)]
    for ref in refs:
        assert ref['url'].startswith('https://')
        assert ref['relevance'] and ref['implementation'] and ref['category']
    assert refs[0]['url'].endswith('/Nose_cone_design')
    assert 'swra705.pdf' in refs[1]['pdf_url']
    assert 'mancini++18-radome-part-i.pdf' in refs[2]['pdf_url']
    assert '13-03-Kouroupis.pdf' in refs[3]['pdf_url']


def test_user_uploaded_pdfs_bundled_and_verified():
    refs=json.loads((BASE/'references/catalog.json').read_text(encoding='utf-8'))['references']
    bundled=[r for r in refs if r['local_pdf']]
    assert len(bundled) == 2
    for ref in bundled:
        p=BASE/ref['local_pdf']
        content=p.read_bytes()
        assert content.startswith(b'%PDF-')
        assert len(content) > 100_000
        assert hashlib.sha256(content).hexdigest()==ref['sha256']


def test_readme_app_and_citations_present():
    readme=(BASE/'README.md').read_text(encoding='utf-8')
    app=(BASE/'app.py').read_text(encoding='utf-8')
    catalog=(BASE/'references/catalog.json').read_text(encoding='utf-8')
    for item in json.loads(catalog)['references']:
        assert item['id'] in readme
    assert 'References & papers' in app
    assert 'references/catalog.json' in readme
    assert '_refs = json.loads' in app
