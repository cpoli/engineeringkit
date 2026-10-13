# Documentation architecture

The structure follows MathematicsKit at commit
`7c075588a3208f5fdaed9cb7dfd41bf0cd01e9e4`.

| Entry point | Contents |
|---|---|
| `source/index.rst` | Engineering overview, hero plots and grouped domain cards |
| `source/subpackages/index.rst` | Twelve domain hubs with History / Examples / API doors |
| `source/history/index.rst` | Grouped cards and `<domain>_breakthroughs.rst` pages |
| `source/examples/index.rst` | Executed galleries plus tutorials and model guides |
| `source/api/index.rst` | Domain API cards, shared utilities and validation pages |

`SUBPACKAGES` in `source/conf.py` is the single source of truth for domain names,
categories and descriptions. A builder hook creates `_generated/vars.rst`,
category grids, individual hub pages and per-domain cross-link strips. All generated
source is ignored by git. API and history pages include the same cross-link strip.

Scripts live in `examples/<domain>/plot_<domain>_models.py`, with a `README.rst`
header for each domain. Sphinx-Gallery runs the scripts and produces
`source/api/gallery/<domain>/`, including figures, thumbnails, downloadable Python
files and notebooks. History pages use `minigallery` links to those scripts.
Failures abort the build. The twelve domain galleries contain 24 scripts demonstrating representative
models; they do not yet walk through every public API feature.

The PyData theme, grouped Sphinx-Design cards, NumPy docstrings, viewcode,
autodoc type descriptions and domain-first navigation follow MathematicsKit.
EngineeringKit also retains its explicit SI/validity guides. JupyterLite is not
configured: the package has not been published, so browser notebooks cannot yet
install a public release as MathematicsKit's notebooks do.

## Build from source

```bash
python -m pip install -e '.[docs]'
MPLBACKEND=Agg python -m sphinx -E -W --keep-going -b html docs/source docs/_build/html
```

Open `docs/_build/html/index.html`. The offline documentation archive places
`index.html` at its root; extract the entire archive before opening it so CSS,
JavaScript, plots and downloadable examples remain available.
