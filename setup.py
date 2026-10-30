import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="LineStacker",
    version='2.0.0',  # Python 3 / CASA 6 port
    author="Jean-Baptiste Jolly, Lukas Lindroos",
    author_email="jbjolly@mpe.mpg.de",
    description="Python module to stack spectra in the image domain.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jbjolly/LineStacker",
    project_urls={
        "Documentation": "https://jbjolly.github.io/LineStacker/",
        "Releases": "https://github.com/jbjolly/LineStacker/releases",
    },
    packages=setuptools.find_packages(exclude=["tests", "tests.*"]),
    install_requires=[
        "numpy>=1.22",
        "scipy>=1.7",
        # astropy < 7.2 breaks with NumPy >= 2.4 (np.in1d removed); both
        # need Python >= 3.11, so older Pythons keep older versions.
        "astropy>=7.2; python_version>='3.11'",
        "astropy>=5.0; python_version<'3.11'",
    ],
    extras_require={
        # Only needed to read/write CASA images (FITS cubes need no CASA).
        # `pip install LineStacker[casa]`
        "casa": [
            "casatools>=6.6",
            "casatasks>=6.6",
        ],
        "plot": ["matplotlib"],
        "test": ["pytest"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)
