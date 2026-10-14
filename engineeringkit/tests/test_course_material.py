"""Documentation contract: the tutorials, the exercise solutions and the
quickstart notebook run, and every printed result matches the output written
under it.

The docs build runs the same ``.. plot::`` blocks (so a broken block fails the
``-W`` build), but it does not compare what they print. Here a block's
expected output is the run of ``# `` comment lines directly under a code
line; a comment after a blank line is ordinary prose. Each page runs top to
bottom in one namespace, reset where a block says ``:context: reset``.

Skipped outside a source checkout (e.g. an installed wheel without docs/).
"""

import contextlib
import io
import json
import re
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "source"
NOTEBOOK = ROOT / "notebooks" / "quickstart.ipynb"
PAGES = sorted((DOCS / "tutorials").glob("*.rst")) + sorted((DOCS / "exercises").glob("*.rst"))

pytestmark = pytest.mark.skipif(not (DOCS / "tutorials").is_dir(), reason="needs a source checkout")

_PLOT = re.compile(r"^(?P<indent> *)\.\. plot::\s*$")


def _plot_blocks(text):
    """Yield ``(options, code)`` for every ``.. plot::`` directive in an rst page."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        match = _PLOT.match(lines[i])
        i += 1
        if not match:
            continue
        indent = len(match.group("indent"))
        body = []
        while i < len(lines) and (not lines[i].strip() or len(lines[i]) - len(lines[i].lstrip()) > indent):
            body.append(lines[i])
            i += 1
        options = {}
        while body and body[0].strip().startswith(":"):
            key, _, value = body.pop(0).strip()[1:].partition(":")
            options[key] = value.strip()
        width = min(len(line) - len(line.lstrip()) for line in body if line.strip())
        yield options, "\n".join(line[width:] for line in body).strip("\n") + "\n"


def _expected_output(code):
    """The ``# `` comment runs that directly follow a code line."""
    expected, previous = [], "blank"
    for line in code.splitlines():
        stripped = line.strip()
        if not stripped:
            previous = "blank"
        elif stripped.startswith("#"):
            if previous in ("code", "output"):
                expected.append(stripped[2:] if stripped.startswith("# ") else stripped[1:])
                previous = "output"
            else:
                previous = "comment"
        else:
            previous = "code"
    return expected


def _run(code, namespace):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        exec(compile(code, "<doc block>", "exec"), namespace)
    plt.close("all")
    return out.getvalue().splitlines()


def test_pages_exist():
    names = {p.stem for p in PAGES}
    assert len([p for p in PAGES if p.parent.name == "tutorials" and p.stem != "index"]) >= 4
    from engineeringkit.tests.test_package import DOMAINS

    assert set(DOMAINS) <= names, f"missing exercise sets: {sorted(set(DOMAINS) - names)}"


@pytest.mark.parametrize("page", [p for p in PAGES if p.stem != "index"], ids=lambda p: f"{p.parent.name}/{p.stem}")
def test_page_runs_and_prints_what_it_says(page):
    blocks = list(_plot_blocks(page.read_text()))
    assert blocks, f"{page.name} has no executed code"
    assert blocks[0][0].get("context") == "reset", f"{page.name}: the first block must reset the context"
    namespace = {}
    for number, (options, code) in enumerate(blocks, 1):
        assert "context" in options, f"{page.name} block {number}: every block shares the page context"
        if options["context"] == "reset":
            namespace = {}
        printed = _run(code, namespace)
        assert printed == _expected_output(code), f"{page.name} block {number}"


def test_quickstart_notebook_runs():
    notebook = json.loads(NOTEBOOK.read_text())
    assert "colab.research.google.com" in "".join(notebook["cells"][0]["source"])
    namespace = {}
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            code = "".join(line for line in cell["source"] if not line.lstrip().startswith(("%", "!")))
            _run(code, namespace)
        assert not cell.get("outputs"), "commit the notebook without outputs"
