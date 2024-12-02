# setup.py

from setuptools import setup, find_packages

setup(
    name='microelektrodica',
    version='1.0.0',
    author='C. Baqueiro Basto, M. Secanell, L.C. Ordoñez',
    author_email='carlosbaqueirob@gmail.com',
    description='A Python Tool for Modeling Microkinetic Electrocatalytic Reactions',
    license='CC BY-NC-SA 4.0',
    packages=find_packages(where='Source'),
    include_package_data=True,
    install_requires=[
        'pytest',
        'scipy',
        'numpy',
        'pandas',
        'matplotlib',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License',
    ],
    python_requires='>=3.6',
)
