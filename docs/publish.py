"""Copy a Sphinx HTML build into this folder (docs/, served by GitHub Pages),
replacing the previously built pages and keeping the files listed in KEEP.
Used by `make html`:  python publish.py <build/html>"""
import os
import shutil
import sys

KEEP = {'source', 'Makefile', 'README.md', 'requirements.txt', 'publish.py', '.nojekyll', 'LineStacker.pdf'}

build = os.path.abspath(sys.argv[1])
here = os.path.dirname(os.path.abspath(__file__))
if not os.path.isfile(os.path.join(build, 'index.html')):
    sys.exit('no index.html in {0}: build the HTML first'.format(build))
for name in os.listdir(here):
    if name in KEEP:
        continue
    path = os.path.join(here, name)
    shutil.rmtree(path) if os.path.isdir(path) and not os.path.islink(path) else os.remove(path)
shutil.copytree(build, here, dirs_exist_ok=True)
open(os.path.join(here, '.nojekyll'), 'a').close()
print('website updated in {0}'.format(here))
