import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blksheep",
    version="0.0.3",
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
        "pandas >= 0.24.2",
        "numpy >= 1.16.4",
        "scipy >= 1.2.1",
        "matplotlib >= 3.1.0",
        "seaborn >= 0.9.0",
        "scikit-learn >= 0.21.2",
        "statsmodels >= 0.10.0"
    ],
    packages=setuptools.find_packages(
        exclude=["*.tests", "*.tests.*", "tests.*", "tests"]
    ),
    entry_points={"console_scripts": ["blacksheep = blacksheep.cli:_main"]},
)
