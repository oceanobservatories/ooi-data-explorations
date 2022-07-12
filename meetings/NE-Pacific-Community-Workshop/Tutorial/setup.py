import setuptools

with open(README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyOOI",
    version="0.0.1a",
    author="Andrew Reed",
    authoer_email="areed@whoi.edu",
    descriptio"="A package for interacting with the Ocean Observatories Initiative (OOI) Machine-2-Machine (M2M) API.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.9.7"
)
