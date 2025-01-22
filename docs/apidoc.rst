.. _api-doc:

API Documentation
==================

μElektrodica is a Python electrochemistry toolbox using the SciPy, NumPy, Networkx, re, os, and Matplotlib libraries. μElektrodica algorithm has been designed to adapt to the wide variety of kinetic data sets available in the literature and to be a tool to simplify the complex task of modeling electrocatalytic reactions. Although it was designed to facilitate the handling of microkinetic models, μElektrodica can also support macrokinetic models; however, methods to estimate the determining step of the reaction have not yet been implemented.

μElektrodica is a Python package built by combining functional and object-oriented programming. The codebase utilizes a modular architecture to promote ease of understanding, adaptability, and future development efforts. It comprises several modules with specific functions for pre-processing, processing, and post-processing of data.

.. toctree::
   :maxdepth: 1

   Collector <apidoc-pages/collector>
   Kpynetic <apidoc-pages/kpynetic>
   Calculator <apidoc-pages/calculator>
   Fitter <apidoc-pages/fitter>
   Coordinator <apidoc-pages/coordinator>
   Grapher <apidoc-pages/grapher>
   Writer <apidoc-pages/writer>