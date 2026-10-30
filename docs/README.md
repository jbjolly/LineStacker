# LineStacker documentation

This folder is the website https://jbjolly.github.io/LineStacker/ (GitHub
Pages serves `docs/` of the main branch). The pages are written in
`source/` (Markdown, plus `api.rst` for the reference generated from the
docstrings); the HTML files next to this README are built from them.

To rebuild after changing the pages or the docstrings:

    pip install -r docs/requirements.txt
    make -C docs html          # the website, into docs/
    make -C docs pdf           # docs/LineStacker.pdf (needs LaTeX)

then commit `docs/` and push: the site updates within a few minutes.

- The examples notebook page is a copy of
  `LineStacker/Example/LineStacker_examples.ipynb`, shown with its saved
  outputs (it is not run during the build): run the notebook first if the
  examples changed.
- `changes.md` includes `PYTHON3_MIGRATION.md` from the repository root.
- The uv-domain pages appear automatically when `LineStacker.uv` is in the
  package (with `UV_STACKING.md` at the repository root).
- The old page names of the 1.x site (`flowCharts.html`, `testos.html`,
  `example.html`, `getting_started.html`) redirect to their new pages.
- CASA is not needed to build the documentation. `make html` replaces the
  built pages in this folder (`publish.py`) and keeps `source/`, this README,
  the Makefile and the PDF.
- Logos: `source/_static/LineStacker-logo.svg` (large: sidebar, README, and
  as `LineStacker-logo.pdf` on the PDF title page) and
  `source/_static/LineStacker-icon.svg` (small: favicon, sidebar on narrow
  screens; `favicon-32.png` and `apple-touch-icon.png` are PNG copies for
  browsers without SVG icons).
