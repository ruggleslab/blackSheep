import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

files = []

setuptools.setup(
    name="transposcope",
    version="0.1.2",
    author="Mark Grivainis",
    author_email="mark.grivainis@fenyolab.org",
    description="A package for visualizing read coverage in areas surrounding novel mobile element insertions.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FenyoLab/transposcope",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Licencse :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Windows",
        "Operating System :: Unix",
    ],
    # packages=["transposcope", "transposcope.viewer", "transposcope.parsers"],
    packages=setuptools.find_packages(),
    package_data={"transposcope": files},
    entry_points={"console_scripts": ["transposcope = transposcope.cli:main"]},
)
