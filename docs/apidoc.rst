.. _api-doc:

API Documentation
==================

*μ*\Elektrodica is a Python electrochemistry toolbox using the SciPy, NumPy, Networkx, re, os, and Matplotlib libraries.
*μ*\Elektrodica algorithm has been designed to adapt to the wide variety of kinetic data sets available in the
literature and to be a tool to simplify the complex task of modeling electrocatalytic reactions.
Although it was designed to facilitate the handling of microkinetic models, *μ*\Elektrodica can also support macrokinetic
models; however, methods to estimate the determining step of the reaction have not yet been implemented.

*μ*\Elektrodica is a Python package built by combining functional and object-oriented programming.
The codebase utilizes a modular architecture to promote ease of understanding, adaptability, and future development efforts.
It comprises several modules with specific functions for pre-processing, processing, and post-processing of data.

.. raw:: latex

    \begin{figure}[htbp]
    \centering
    \begin{tikzpicture}[node distance=2cm]
    \node (start) [startstop] {\textbf{\textit{\textmu{}}Elektrodica}};

    \node (collector) [process, below of=start] {\textbf{Collector}};

    \node (collector2) [process, left of=collector, xshift=-2cm] {
        \begin{minipage}{3cm}
            \centering
            \small{DataParameter} \\
            \small{DataSpecies} \\
            \small{DataReactions}
        \end{minipage}
    };

    \node (kpynetic) [process, below of=collector] {\textbf{Kpynetic}};
    \node (kpynetic2) [process, left of=kpynetic, xshift=-2cm] {
        \begin{minipage}{3cm}
            \centering
            \small{FreeEnergy} \\
            \small{RateConstants} \\
            \small{ReactionRate}
        \end{minipage}
    };
    \node (split) [circle, minimum size=2mm, inner sep=0pt, draw=black, fill=black, below of=kpynetic, yshift=0.7cm] {};
    \node (fitter) [process, below of=kpynetic, yshift=-0.5cm] {\textbf{Fitter}};
    \node (coordinator) [process, below right of=kpynetic, xshift=3cm, yshift=-1cm] {\textbf{Coordinator}};
    \node (calculator) [process, below left of=kpynetic, xshift=-3cm, yshift=-1cm] {
        \textbf{Calculator}};
    \node (merge) [circle, minimum size=2mm, inner sep=0pt, draw=black, fill=black, below of=fitter, yshift=0.7cm] {};
    \node (grapher) [startstop, below left of=fitter, xshift=-1cm, yshift=-1cm] {\textbf{Grapher}};
    \node (writer) [startstop, below right of=fitter, xshift=1cm, yshift=-1cm] {\textbf{Writer}};
    \draw [arrow] (start) -- (collector) node[midway, right] {directory};
    \draw [arrow] (collector) -- (kpynetic) node[midway, right] {data};
    \draw [arrow] (collector2) -- (collector);
    \draw [arrow] (kpynetic2) -- (kpynetic);
    \draw [arrow] (kpynetic) -- (split) node[midway, right] {Kpy};
    \draw [arrow] (split) -- (coordinator);
    \draw [arrow] (split) -- (calculator);
    \draw [arrow] (split) -- (fitter);
    \draw [arrow] (calculator) -- (fitter);
    \draw [arrow] (coordinator) -- (fitter);
    \draw [arrow] (coordinator) -- (merge);
    \draw [arrow] (calculator) -- (merge);
    \draw [arrow] (fitter) -- (merge);
    \draw [arrow] (merge) -- (grapher) node[midway, yshift=-0.1cm, right] {results};
    \draw [arrow] (merge) -- (writer);
    \end{tikzpicture}
    \caption{\textit{\textmu{}}Elektrodica class diagram. Some classes omitted for readability.}
    \label{fig:diagrama_flux}
    \end{figure}

.. toctree::
   :maxdepth: 1
   :caption: Modules:

   Collector <apidoc-pages/collector>
   Kpynetic <apidoc-pages/kpynetic>
   Calculator <apidoc-pages/calculator>
   Fitter <apidoc-pages/fitter>
   Coordinator <apidoc-pages/coordinator>
   Grapher <apidoc-pages/grapher>
   Writer <apidoc-pages/writer>
   Tools <apidoc-pages/tools>