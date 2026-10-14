#!/usr/bin/env python3
"""Check that saved notebook links and inspection-page assets are in Git.

Run after staging notebook outputs and their referenced evidence. New local
runs remain ignored; only evidence linked by the published notebooks is needed.
"""
import json
from html.parser import HTMLParser
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit

PACKAGE = Path(__file__).resolve().parents[1]


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = set()

    def handle_starttag(self, tag, attrs):
        self.urls.update(value for key, value in attrs
                         if key in ("href", "src") and value)


def references(path):
    parser = Links()
    if path.suffix == ".ipynb":
        for cell in json.loads(path.read_text(encoding="utf-8"))["cells"]:
            if cell["cell_type"] == "markdown":
                source = "".join(cell["source"])
                parser.urls.update(re.findall(r"\]\(([^)]+)\)", source))
                parser.feed(source)
            for output in cell.get("outputs", []):
                for mime in ("text/html", "image/svg+xml"):
                    parser.feed("".join(output.get("data", {}).get(mime, "")))
    else:
        parser.feed(path.read_text(encoding="utf-8"))
    for url in parser.urls:
        parts = urlsplit(url)
        if parts.scheme or parts.netloc or not parts.path:
            continue
        yield (path.parent / unquote(parts.path)).resolve()


def linked_files(notebook_paths=None):
    notebooks = PACKAGE / "notebooks"
    pending = [notebooks / "paper_examples.ipynb", *sorted((notebooks / "examples").glob("[0-9][0-9]_*.ipynb")),
               *sorted((notebooks / "visual/comparisons").glob("*.ipynb"))]
    if notebook_paths is not None:
        pending = list(notebook_paths)
    seen, targets = set(), set()
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        for target in references(path):
            targets.add(target)
            # Follow linked notebooks and their exported inspection pages.
            if (target.is_file() and (
                    (target.suffix == ".html" and target.is_relative_to(notebooks / "visual/inspection"))
                    or (target.suffix == ".ipynb" and target.is_relative_to(notebooks)))):
                pending.append(target)
    return targets


def main():
    REPO = Path(subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"], cwd=PACKAGE, text=True).strip())
    tracked = {REPO / name for name in subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=REPO, text=True).split("\0") if name}
    targets = linked_files()
    missing = sorted(path for path in targets if not path.is_file())
    untracked = sorted(targets - tracked - set(missing))
    for label, paths in (("Missing", missing), ("Not in Git", untracked)):
        for path in paths:
            print(f"{label}: {path.relative_to(REPO) if path.is_relative_to(REPO) else path}")
    print(f"{len(targets)} linked files; {len(missing)} missing; {len(untracked)} not in Git.")
    raise SystemExit(bool(missing or untracked))


if __name__ == "__main__":
    main()
