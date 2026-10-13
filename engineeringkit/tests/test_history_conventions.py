"""Documentation contract: every domain's history page lists at least 15
breakthroughs, each linked to its own gallery example(s), and every gallery
example belongs to exactly one breakthrough.

Skipped outside a source checkout (e.g. an installed wheel without docs/).
"""

import re
from pathlib import Path

import pytest

from engineeringkit.tests.test_package import DOMAINS

ROOT = Path(__file__).resolve().parents[2]
HISTORY = ROOT / "docs" / "source" / "history"
EXAMPLES = ROOT / "examples"

pytestmark = pytest.mark.skipif(not HISTORY.is_dir() or not EXAMPLES.is_dir(), reason="needs a source checkout")

MINIMUM_BREAKTHROUGHS = 15
_HEADING = re.compile(r"^(?P<title>[^\n]+)\n-{3,}\n", re.M)


def _breakthroughs(domain):
    text = (HISTORY / f"{domain}_breakthroughs.rst").read_text()
    matches = list(_HEADING.finditer(text))
    sections = []
    for i, match in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        body = text[match.end() : end]
        links = [p for line in re.findall(r"^\.\. minigallery::(.*)$", body, re.M) for p in line.split()]
        if match.group("title").strip() != "See Also":
            sections.append((match.group("title").strip(), links))
    return sections


@pytest.mark.parametrize("domain", DOMAINS)
def test_history_has_enough_breakthroughs_each_with_a_gallery(domain):
    sections = _breakthroughs(domain)
    assert len(sections) >= MINIMUM_BREAKTHROUGHS, f"{domain}: only {len(sections)} breakthroughs"
    for title, links in sections:
        assert links, f"{domain}: '{title}' has no minigallery example"


@pytest.mark.parametrize("domain", DOMAINS)
def test_every_example_is_linked_from_exactly_one_breakthrough(domain):
    linked = []
    for _, links in _breakthroughs(domain):
        linked += [(ROOT / "docs" / "source" / link).resolve() for link in links]
    assert len(linked) == len(set(linked)), f"{domain}: an example is linked from two breakthroughs"
    for path in linked:
        assert path.is_file(), f"{domain}: minigallery target {path} does not exist"
    on_disk = {p.resolve() for p in (EXAMPLES / domain).rglob("plot_*.py")}
    assert on_disk == set(linked), f"{domain}: unlinked examples {sorted(map(str, on_disk - set(linked)))}"
