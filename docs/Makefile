# Build the LineStacker documentation.
#   make html   -> the website, built in ../build/docs/html and copied into this
#                  folder (docs/, served by GitHub Pages)
#   make pdf    -> docs/LineStacker.pdf (needs LaTeX: latexmk, pdflatex)
#   make clean  -> remove the build files (the website stays)
# Requirements: pip install -r docs/requirements.txt

SPHINXBUILD ?= sphinx-build
SOURCE      = source
BUILD       = ../build/docs
PYTHON      ?= python

.PHONY: html pdf clean

html:
	$(SPHINXBUILD) -b html $(SOURCE) $(BUILD)/html -d $(BUILD)/doctrees
	$(PYTHON) publish.py $(BUILD)/html

pdf:
	$(SPHINXBUILD) -b latex $(SOURCE) $(BUILD)/latex -d $(BUILD)/doctrees
	$(MAKE) -C $(BUILD)/latex all-pdf LATEXMKOPTS="-interaction=nonstopmode -halt-on-error"
	cp $(BUILD)/latex/LineStacker.pdf LineStacker.pdf

clean:
	rm -rf $(BUILD)
