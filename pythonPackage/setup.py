import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

files = ["*.csv", "*.tsv"]

setuptools.setup(
    name="blacksheep-outliers",
    version="0.0.1",
    author="Ruggles Lab",
    author_email="lmb529@nyumc.org",
    description="A package for differential expression analysis using outlier values.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ruggleslab/blackSheep/pythonPackage",
    # project_urls={
    #         "Bug Tracker": "https://bugs.example.com/HelloWorld/",
    #         "Documentation": "https://docs.example.com/HelloWorld/",
    #         "Source Code": "https://code.example.com/HelloWorld/",
    #     },
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    install_requires=[
        'pandas >= 0.24.2',
        'numpy >= 1.16.4',
        'scipy >= 1.2.1',
        'matplotlib >= 3.1.0',
        'seaborn >= 0.9.0',
        'pickle',
        'catheat'
    ],
    dependency_links =[
        'git+git://github.com/schlegelp/catheat@master'
    ],
    packages=setuptools.find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_data={"tests": files},
    entry_points={"console_scripts": ["BlackSheep = blacksheep.cli:main"]},
)
