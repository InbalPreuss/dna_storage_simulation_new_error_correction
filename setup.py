import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dna_storage",
    version="0.0.0",
    author="Inbal Preuss",
    author_email="inbalpreuss@gmail.com",
    description="DNA Storage using shortmer algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/InbalPreuss/DnaStorage",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
       'numpy',
       'pandas',
        'python-Levenshtein',
        'biopython'
    ]
)
