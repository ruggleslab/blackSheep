import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blksheep",
    version="0.0.7",
    author="Ruggles Lab",
    author_email="ruggleslab@gmail.com",
    description="A package for differential extreme values analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ruggleslab/blackSheep/",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    install_requires=[
        "matplotlib >= 3.3.4",
        "numpy >= 1.20.1",
        "pandas >= 1.2.2",
        "scipy >= 1.6.0",
        "seaborn >= 0.11.1",
        "statsmodels >= 0.12.2"
    ],
    packages=setuptools.find_packages(
        exclude=["*.tests", "*.tests.*", "tests.*", "tests"]
    ),
    entry_points={"console_scripts": ["blacksheep = blacksheep.cli:_main"]},
)
