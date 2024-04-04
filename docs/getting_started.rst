Getting Started
===============

This page details how to get started with montecarlo_jg. montecarlo_jg is a package which was developed for the CHEM 3684 Quantum Software class
at Virginia Tech.

Installation
------------
To install molecool, you will need an environment with the following packages:

* Python 3.11
* NumPy
* Matplotlib
* networkx

Once you have these packages installed, you can install molecool in the same environment using
::

    pip install -e .


Usage
-------
Once installed, you can use the package. This example creates a new bit string, which represents the spins of each node in the graph.
.. code-block:: python


    import montecarlo_jg

    bitstring = montecarlo_jg.BitString(10)

