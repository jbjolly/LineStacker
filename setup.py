import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="LineStacker",
    version='1.5.2',
    author="Jean-Baptiste Jolly, Lukas Lindroos",
    author_email="jean.jolly@chalmers.se",
    description="Python module to stack spectra in the image domain.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jbjolly/LineStacker/releases",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
)
