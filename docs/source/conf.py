"""MathematicsKit-style domain hubs, cross-links and executed example galleries."""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from engineeringkit import __version__

project = "engineeringkit"
author = "Charles Poli"
release = __version__
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.doctest",
    "sphinx.ext.viewcode",
    "sphinx.ext.extlinks",
    "sphinx_autodoc_typehints",
    "sphinx_gallery.gen_gallery",
    "sphinx_design",
    "matplotlib.sphinxext.plot_directive",
    "jupyterlite_sphinx",
    "sphinx_sitemap",
    "sphinxext.opengraph",
]
napoleon_google_docstring = False
napoleon_numpy_docstring = True
autodoc_default_options = {"members": True, "undoc-members": False, "show-inheritance": True}
autodoc_typehints = "description"
templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["engineeringkit.css"]
html_theme = "pydata_sphinx_theme"

# Search-engine metadata: sphinx-sitemap writes sitemap.xml with absolute
# URLs from html_baseurl, and sphinxext-opengraph adds a <meta
# name="description"> plus og:* link-preview tags from each page's text.
# html_baseurl also gives every page a <link rel="canonical">.
html_baseurl = "https://cpoli.github.io/engineeringkit/"
html_title = project
sitemap_url_scheme = "{link}"
sitemap_excludes = [
    "search.html",
    "genindex.html",
    "py-modindex.html",
    "*sg_execution_times.html",
    "_modules/*",
    # Include-only fragments (card grids, nav strips), not standalone pages.
    "_generated/grid_*",
    "_generated/vars.html",
    "_generated/nav/*",
]
ogp_site_url = html_baseurl
ogp_site_name = project
ogp_image = html_baseurl + "_static/images/readme_hero.png"
ogp_social_cards = {"enable": False}

# Page-view counts with GoatCounter (no cookies, no personal data). One
# site, cpoli.goatcounter.com, covers every cpoli.github.io project; the
# URL path tells the kits apart.
html_js_files = [
    (
        "https://gc.zgo.at/count.js",
        {"async": "async", "data-goatcounter": "https://cpoli.goatcounter.com/count"},
    ),
]
html_logo = "_static/images/engineeringkit_logo_transparent.png"
html_theme_options = {
    "icon_links": [
        {
            "name": "Charles Poli",
            "url": "https://cpoli.github.io/",
            "icon": "fa-solid fa-user",
            "type": "fontawesome",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "show_toc_level": 2,
    "navigation_with_keys": True,
    "navigation_depth": 2,
    "logo": {
        "alt_text": "engineeringkit logo",
    },
}
exclude_patterns = ["_build", "build", "_generated/nav/*", "_generated/grid_*.rst", "_generated/vars.rst"]

SUBPACKAGES = [
    {
        "name": "structural_engineering",
        "category": "Structures and Machines",
        "blurb": (
            "Members, sections, beams on any supports, columns, stress and strain at a point, "
            "pressure vessels, energy methods, trusses, frames, cables and arches."
        ),
    },
    {
        "name": "mechanical_engineering",
        "category": "Structures and Machines",
        "blurb": "Rigid-body motion, Lagrange's equations, machine elements, mechanisms, cams and multi-degree-of-freedom vibration.",
    },
    {
        "name": "thermodynamics",
        "category": "Energy and Flow",
        "blurb": (
            "Ideal and real gases, mixtures, the two laws, steady-flow devices, engine and refrigeration cycles, "
            "exergy, steam properties, psychrometrics and combustion."
        ),
    },
    {
        "name": "fluid_engineering",
        "category": "Energy and Flow",
        "blurb": "Hydrostatics, flow measurement, momentum, pipe friction, drag, similarity, open channels, pumps and networks.",
    },
    {"name": "heat_transfer", "category": "Energy and Flow", "blurb": "Conduction, convection, radiation, fins and heat exchangers."},
    {"name": "aerodynamics", "category": "Energy and Flow", "blurb": "Airfoil and wing baselines, drag polars and compressible flow."},
    {
        "name": "electrical_engineering",
        "category": "Circuits, Signals and Control",
        "blurb": "DC and AC circuits, transients, electronics, machines and power systems.",
    },
    {
        "name": "control_systems",
        "category": "Circuits, Signals and Control",
        "blurb": "Classical stability, compensator design, PID with anti-windup, nonlinear loops, state feedback, LQR and estimation.",
    },
    {
        "name": "signal_processing",
        "category": "Circuits, Signals and Control",
        "blurb": "Sampling and rate conversion, the FFT and DCT, spectra and AR models, IIR and FIR design, adaptive filters.",
    },
    {
        "name": "communications",
        "category": "Circuits, Signals and Control",
        "blurb": "Modulation and pulse shaping, noise, capacity, detection, PCM, source and channel coding, fading and link budgets.",
    },
    {
        "name": "materials",
        "category": "Materials and Reliability",
        "blurb": (
            "Crystals and diffusion, elastic constants, hardness, yield criteria, composites and "
            "laminates, fatigue and notches, fracture, creep and dislocation plasticity."
        ),
    },
    {
        "name": "reliability",
        "category": "Materials and Reliability",
        "blurb": "Lifetimes and life data, accelerated testing, networks, structural reliability, FMEA, fault trees and maintenance.",
    },
]


def _jupyterlite_install_cell(notebook_content, notebook_filename):
    """Prepend a ``%pip install`` cell so the Pyodide kernel fetches engineeringkit from PyPI."""
    notebook_content["cells"].insert(
        0,
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": ["%pip install engineeringkit"],
        },
    )


sphinx_gallery_conf = {
    "examples_dirs": [f"../../examples/{s['name']}" for s in SUBPACKAGES],
    "gallery_dirs": [f"api/gallery/{s['name']}" for s in SUBPACKAGES],
    "filename_pattern": r"/plot_",
    "download_all_examples": True,
    "within_subsection_order": "FileNameSortKey",
    "backreferences_dir": "gen_modules/backreferences",
    "doc_module": ("engineeringkit",),
    "abort_on_example_error": True,
    "matplotlib_animations": (True, "jshtml"),
    # "Launch JupyterLite" button: runs the notebook in the browser on
    # Pyodide, with no server.
    "jupyterlite": {
        "use_jupyter_lab": True,
        "notebook_modification_function": _jupyterlite_install_cell,
    },
    # Run the examples in one worker process per core; ENGINEERINGKIT_GALLERY_JOBS=1
    # runs them one at a time in the Sphinx process instead.
    "parallel": int(os.environ.get("ENGINEERINGKIT_GALLERY_JOBS", os.cpu_count() or 1)),
}

# The quickstart notebook also ships in the JupyterLite site, so the "Try it
# online" links open on it; sphinx-gallery appends its own notebooks to this
# list. Its first cell installs engineeringkit with %pip, which works on Pyodide.
jupyterlite_contents = ["../../notebooks/quickstart.ipynb"]

# The tutorials and the exercise solutions are `.. plot::` blocks: the build
# runs them (an exception is a warning, so -W fails), shows their source and
# embeds their figures. Each page shares one namespace through `:context:`.
# engineeringkit/tests/test_course_material.py runs the same blocks and
# checks the output written under each print.
plot_include_source = True
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = [("png", 110)]
plot_pre_code = ""

# :doi:`10.xxxx/yyyy` renders as the plain "doi:10.xxxx/yyyy" text the
# history pages cite, but as a link, so `make linkcheck` checks it.
extlinks = {"doi": ("https://doi.org/%s", "doi:%s")}

# make linkcheck (.github/workflows/linkcheck.yml, weekly). linkcheck_ignore
# matches the URL as written, so a publisher that answers bots with 403 is
# skipped by its DOI registrant prefix; DOIs from other publishers are
# followed through doi.org to the article. Anchors on the generated pages
# are not real headings, so only each page itself is checked.
linkcheck_ignore = [
    r"(\.\./)*lite/",  # the JupyterLite site, written by the html build
    r"https://doi\.org/10\.1002/",  # Wiley
    r"https://doi\.org/10\.1021/",  # ACS
    r"https://doi\.org/10\.1061/",  # ASCE
    r"https://doi\.org/10\.1098/",  # Royal Society
    r"https://doi\.org/10\.1145/",  # ACM
]
linkcheck_anchors = False
linkcheck_timeout = 30
linkcheck_retries = 3
linkcheck_rate_limit_timeout = 120

# The JupyterLite hook above makes sphinx_gallery_conf unpicklable, which only
# disables Sphinx's config cache.
suppress_warnings = ["config.cache"]


def _cards(entries):
    text = ".. grid:: 1 2 3 3\n   :gutter: 2\n\n"
    for title, target, description in entries:
        text += f"   .. grid-item-card:: {title}\n      :link: {target}\n      :link-type: doc\n\n      {description}\n\n"
    return text


def _generate_domain_docs(app):
    generated = Path(app.srcdir) / "_generated"

    def write(relative, text):
        target = generated / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)

    write("vars.rst", f".. |num_subpackages| replace:: {len(SUBPACKAGES)}\n")
    targets = {
        "hub": lambda n: f"/_generated/subpackages/{n}",
        "history": lambda n: f"/history/{n}_breakthroughs",
        "examples": lambda n: f"/api/gallery/{n}/index",
        "api": lambda n: f"/api/{n}",
    }
    categories = list(dict.fromkeys(s["category"] for s in SUBPACKAGES))
    for kind, target in targets.items():
        chunks = []
        for category in categories:
            chunks.append(f".. rubric:: {category}\n\n")
            chunks.append(
                _cards([(f"engineeringkit.{s['name']}", target(s["name"]), s["blurb"]) for s in SUBPACKAGES if s["category"] == category])
            )
        write(f"grid_{kind}.rst", "\n".join(chunks))
    for domain in SUBPACKAGES:
        name = domain["name"]
        title = f"engineeringkit.{name}"
        entries = [
            ("History", targets["history"](name), "Foundational ideas linked to the implementation."),
            ("Examples", targets["examples"](name), "Executed scripts, plots, source and notebook downloads."),
            ("API reference", targets["api"](name), "Public classes, functions and result objects."),
        ]
        write(
            f"subpackages/{name}.rst",
            f"{title}\n{'=' * len(title)}\n\n{domain['blurb']}\n\n"
            + _cards(entries)
            + f"\n:doc:`Units, assumptions and numerical validation </guides/{name}>`.\n",
        )
        nav = [
            f":doc:`Domain hub <{targets['hub'](name)}>`",
            f":doc:`History <{targets['history'](name)}>`",
            f":doc:`Examples <{targets['examples'](name)}>`",
            f":doc:`API reference <{targets['api'](name)}>`",
        ]
        write(f"nav/{name}.rst", ".. container:: subpkg-nav\n\n   " + " · ".join(nav) + "\n")


def setup(app):
    app.connect("builder-inited", _generate_domain_docs)
