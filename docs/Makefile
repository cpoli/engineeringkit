# Minimal makefile for Sphinx documentation

SPHINXOPTS    ?= -E -W --keep-going
# Prefer the project's .venv, which has the docs extra (joblib included),
# over whatever sphinx-build is first on PATH.
SPHINXBUILD   ?= $(if $(wildcard ../.venv/bin/sphinx-build),../.venv/bin/sphinx-build,sphinx-build)
SOURCEDIR     = source
BUILDDIR      = build

export MPLBACKEND ?= Agg

.PHONY: help html-noplot Makefile

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Builds the pages without running the gallery examples, to check text,
# cross-references and toctrees in seconds.
html-noplot:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) -D plot_gallery=0 $(O)

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
