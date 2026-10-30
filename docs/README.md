# Documentation architecture

The structure follows MathematicsKit at commit
`7c075588a3208f5fdaed9cb7dfd41bf0cd01e9e4`.

| Entry point | Contents |
|---|---|
| `source/index.rst` | Engineering overview, hero plots and grouped domain cards |
| `source/subpackages/index.rst` | Twelve domain hubs with History / Examples / API doors |
| `source/history/index.rst` | Grouped cards and `<domain>_breakthroughs.rst` pages |
| `source/examples/index.rst` | Executed galleries plus tutorials, exercises and model guides |
| `source/tutorials/` | Cross-domain tutorials, one problem each |
| `source/exercises/` | One exercise set per domain, solutions in drop-downs |
| `source/api/index.rst` | Domain API cards, shared utilities and validation pages |

`SUBPACKAGES` in `source/conf.py` is the single source of truth for domain names,
categories and descriptions. A builder hook creates `_generated/vars.rst`,
category grids, individual hub pages and per-domain cross-link strips. All generated
source is ignored by git. API and history pages include the same cross-link strip.

Scripts live in `examples/<domain>/<topic>/plot_NN_<breakthrough>.py`, with a
`README.rst` header for each domain and each topic subsection. Sphinx-Gallery runs
the scripts and produces `source/api/gallery/<domain>/`, including figures,
thumbnails, downloadable Python files and notebooks. Each history page lists at
least 15 dated breakthroughs, and each breakthrough uses `minigallery` links to its
own scripts; every script is titled after its breakthrough, ends in a check of the
result it illustrates and is linked from exactly one breakthrough
(`engineeringkit/tests/test_history_conventions.py` enforces the linking).
Failures abort the build. The twelve domain galleries contain 213 scripts.

The PyData theme, grouped Sphinx-Design cards, NumPy docstrings, viewcode,
autodoc type descriptions and domain-first navigation follow MathematicsKit.
EngineeringKit also retains its explicit SI/validity guides.

Tutorials and exercise solutions are matplotlib `.. plot::` blocks sharing one
namespace per page (`:context: reset` on the first block, `:context: close-figs`
after it). The build runs them, so an exception fails the `-W` build. Under
each `print`, write the output as `# ` comment lines with no blank line between;
`engineeringkit/tests/test_course_material.py` runs every block and checks
those lines. A comment after a blank line is ordinary prose.

`notebooks/quickstart.ipynb` (committed without outputs) and every gallery
notebook are served in JupyterLite under `lite/`. Their first cell installs
engineeringkit from PyPI, so they need a published release.

## Build from source

```bash
python -m pip install -e '.[docs]'
cd docs && make html
```

Open `docs/build/html/index.html`. The offline documentation archive places
`index.html` at its root; extract the entire archive before opening it so CSS,
JavaScript, plots and downloadable examples remain available.
