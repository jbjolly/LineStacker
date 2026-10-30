# Sphinx configuration of the LineStacker documentation.
#
# Build the website (copied into docs/, the folder GitHub Pages serves):
#     make -C docs html
# and the PDF (docs/LineStacker.pdf):
#     make -C docs pdf
# Requirements: pip install -r docs/requirements.txt

import importlib.util
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)

import LineStacker  # noqa: E402

project = 'LineStacker'
author = 'Jean-Baptiste Jolly'
copyright = '2019-2026, Jean-Baptiste Jolly'
release = LineStacker.__version__
version = release

extensions = [
    'myst_nb',                     # Markdown pages and the examples notebook
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]

# -- generated pages ---------------------------------------------------------
# The examples notebook lives in LineStacker/Example; LineStacker.uv (2.1) gets
# its pages only when it is part of the package being documented.
GENERATED = os.path.join(HERE, '_generated')
os.makedirs(GENERATED, exist_ok=True)
shutil.copyfile(os.path.join(ROOT, 'LineStacker', 'Example', 'LineStacker_examples.ipynb'),
                os.path.join(GENERATED, 'examples_notebook.ipynb'))


def _uv_guide():
    for path in (os.path.join(ROOT, 'UV_STACKING.md'), os.path.join(ROOT, 'LineStacker', 'uv', 'UV_STACKING.md'),
                 os.path.join(ROOT, 'docs', 'UV_STACKING.md')):
        if os.path.isfile(path):
            return path
    return None


HAS_UV = importlib.util.find_spec('LineStacker.uv') is not None
if not HAS_UV and os.path.exists(os.path.join(GENERATED, 'uv.md')):
    os.remove(os.path.join(GENERATED, 'uv.md'))
with open(os.path.join(GENERATED, 'uv_toc.md'), 'w') as fh:
    if HAS_UV:
        guide = _uv_guide()
        with open(os.path.join(GENERATED, 'uv.md'), 'w') as page:
            if guide:
                page.write('```{{include}} {0}\n```\n\n'.format(os.path.relpath(guide, GENERATED)))
            else:
                page.write('# uv-domain stacking: LineStacker.uv\n\n')
            page.write('## LineStacker.uv reference\n\n```{eval-rst}\n.. automodule:: LineStacker.uv\n'
                       '   :members:\n```\n')
        fh.write('```{toctree}\n:maxdepth: 1\n:caption: uv-domain stacking\n\n_generated/uv\n```\n')

# -- general -----------------------------------------------------------------
source_suffix = {'.rst': 'restructuredtext', '.md': 'myst-nb', '.ipynb': 'myst-nb'}
exclude_patterns = ['_redirects', '_generated/uv_toc.md']
templates_path = ['_templates']   # layout.html: PNG icons next to the SVG favicon
pygments_style = 'sphinx'
nitpicky = False

myst_enable_extensions = ['colon_fence', 'deflist', 'dollarmath']
myst_heading_anchors = 3
nb_execution_mode = 'off'           # the notebook is saved with its outputs
nb_merge_streams = True

# CASA and the uv-domain dependencies are not needed to build the documentation.
autodoc_mock_imports = ['casatools', 'casatasks', 'casaconfig', 'casacore', 'numba']
autodoc_member_order = 'bysource'
autodoc_typehints = 'none'
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# -- HTML --------------------------------------------------------------------
html_theme = 'alabaster'
html_static_path = ['_static']
html_extra_path = ['_redirects']    # old page names of the 1.x site
html_title = 'LineStacker {0}'.format(release)
html_favicon = '_static/LineStacker-icon.svg'   # the small logo
html_theme_options = {
    'logo': 'LineStacker-logo.svg',                # the large logo, at the top of the sidebar
    'description': 'Stacking of spectral lines in data cubes and spectra',
    'github_user': 'jbjolly',
    'github_repo': 'LineStacker',
    'github_button': False,
    'fixed_sidebar': True,
    'page_width': '1080px',
    'sidebar_width': '230px',
    'font_family': "-apple-system, 'Segoe UI', 'Helvetica Neue', Arial, sans-serif",
    'head_font_family': "-apple-system, 'Segoe UI', 'Helvetica Neue', Arial, sans-serif",
    'code_font_size': '0.85em',
}
html_sidebars = {'**': ['smallicon.html', 'about.html', 'navigation.html', 'relations.html', 'searchbox.html']}
html_show_sourcelink = False
html_copy_source = False

# -- LaTeX (PDF) -------------------------------------------------------------
latex_documents = [('index', 'LineStacker.tex', 'LineStacker documentation', author, 'manual')]
latex_logo = '_static/LineStacker-logo.pdf'   # vector copy of LineStacker-logo.svg
latex_elements = {'papersize': 'a4paper', 'pointsize': '10pt',
                  'fontpkg': r'\usepackage{mathptmx}\usepackage[scaled=.9]{helvet}\usepackage{courier}',
                  # characters of the pages that pdflatex does not know
                  'preamble': r'\usepackage{newunicodechar}\newunicodechar{σ}{\ensuremath{\sigma}}'
                              r'\newunicodechar{λ}{\ensuremath{\lambda}}\newunicodechar{≈}{\ensuremath{\approx}}'
                              r'\newunicodechar{±}{\ensuremath{\pm}}\newunicodechar{°}{\ensuremath{^\circ}}'
                              r'\newunicodechar{√}{\ensuremath{\surd}}\newunicodechar{²}{\ensuremath{^2}}'
                              r'\newunicodechar{×}{\ensuremath{\times}}'}
