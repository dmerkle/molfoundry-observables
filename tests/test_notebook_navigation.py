"""Notebook links must reach actual cells or the separately served website."""
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parents[1]


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.links.append(dict(attrs))


def test_notebook_section_links_reach_cells_and_open_web_pages_externally():
    for path in (ROOT/'notebooks/examples').glob('[0-9][0-9]_*.ipynb'):
        book = json.loads(path.read_text(encoding="utf-8"))
        navigation = next(c for c in book['cells'] if c['id'] == 'section-nav-links')
        parser = Links()
        parser.feed(''.join(navigation['source']))
        for link in parser.links:
            url = urlsplit(link['href'])
            if url.scheme:
                assert url.scheme == 'https' and url.netloc == 'stochmech.algochem.techfak.de'
                assert link['target'] == '_blank'
                assert 'noopener' in link['rel']
            else:
                assert url.fragment.startswith('cell-id='), (path, link)
                target = json.loads((path.parent/url.path).read_text(encoding="utf-8")) if url.path else book
                identifier = url.fragment.removeprefix('cell-id=')
                assert any(c['id'] == identifier and ''.join(c['source']).strip()
                           for c in target['cells']), (path, link)
