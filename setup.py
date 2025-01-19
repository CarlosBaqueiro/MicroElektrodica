from setuptools import setup, find_packages

setup(
    name="melektrodica",
    version="1.0.0",
    author="C. Baqueiro Basto, M. Secanell, L.C. Ordoñez",
    author_email="carlosbaqueirob@gmail.com",
    description="A Python Electrochemistry Toolbox for Modeling Microkinetic Electrocatalytic Reactions",
    license="CC BY-NC-SA 4.0",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    ],
    packages=find_packages(where="melektrodica"),
    package_dir={"": "melektrodica"},
    include_package_data=True,
    install_requires=[
        "pytest",
        "scipy",
        "numpy",
        "pandas",
        "matplotlib",
        "networkx",
        "copy",
        "warnings",
        "os",
        "re",
    ],
    python_requires=">=3.6",
)
