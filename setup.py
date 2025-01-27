from setuptools import setup, find_packages

setup(
    name="melektrodica",
    version="1.0.0",
    author="C. Baqueiro Basto, M. Secanell, L.C. Ordoñez",
    author_email="carlosbaqueirob@gmail.com",
    url="https://github.com/CarlosBaqueiro/MicroElektrodica",
    description="A Python Electrochemistry Toolbox for Modeling Microkinetic Electrocatalytic Reactions",
    license="CC BY-NC-SA 4.0",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(exclude=["tests", "tests.*"]),
    package_dir={"melektrodica": "melektrodica"},
    include_package_data=True,
    install_requires=[
        "pytest",
        "scipy",
        "numpy",
        "pandas",
        "matplotlib",
        "networkx",
        "colorlog",
        "tabulate",
        "IPython",
    ],
    python_requires=">=3.6",
)
