.. interlab documentation master file, created by
   sphinx-quickstart on Fri Aug 24 10:01:44 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Interlab: A python module for consensus analysis in interlaboratory studies
===========================================================================

`David A. Sheen <mailto:david.sheen@nist.gov>`_

*National Institute of Standards and Technology*

*Last Revision* |today|

`Download this software from GitHub <https://github.com/usnistgov/interlab_py>`_

Welcome to the home page for analysis software of Interlab: A python module for consensus analysis in interlaboratory studies. This software is a Python package that will perform consensus analysis on spectral data such as NMR, GC-MS and LC-MS. Use of this code allows researchers to identify laboratories producing data closest to the consensus values, thereby ensuring that untargeted studies are using the most precise data available to them. The software was originally developed for analyzing NMR data [1, 2] but can be applied to any array data, including Raman or FTIR spectroscopy and GC-MS or LC-MS. Details on the implementation of the code can be found in Ref. [1].
The input for the code consists of a set of sample labels identifying the physical objects measured in the interlaboratory study, facility labels that identify the facility of origin of the measurements, and the data themselves. It is the responsibility of the user to format the data and metadata so that the code can read it. In addition, the user must specify the distance function that will be used to compare the spectra and the statistical distribution that these distances will be fit to. One example is given in the `example notebook <https://pages.nist.gov/interlab_py/analysis_demo.html>`_.
Given the input data, the code will perform the following tasks:

  * Calculates the interspectral distances
  * Fits the project’s distribution function to the distance data and calculate the corresponding scores.
  * Identifies outliers within each spectral population
  * Conducts a principal components analysis on the scores and compute the projected statistical distance
  * Uses the projected statistical distance to determine the data set outliers.

The software cannot be used out of the box. Users must create an interface to their own software, and that interface will be specific to the user’s application. The `example notebook <https://pages.nist.gov/interlab_py/analysis_demo.html>`_, demonstrates one such possible interface.


Contents
++++++++

.. toctree::
   :maxdepth: 2
   :includehidden:

   Project
   Subobjects
   Metrics
   analysis_demo.ipynb
..   analysis_demo
   :caption: Contents:

Legal
+++++

This software is subject to the `NIST Software License <https://www.nist.gov/director/licensing>`_ (revised as of August 2018).


Contact
+++++++

David Sheen

* `Email <mailto:david.sheen@nist.gov>`_
* `NIST Staff Page <https://www.nist.gov/people/david-sheen>`_
* `GitHub profile <https://github.com/davidasheen>`_

Links
+++++

`NIST GitHub Organization <https://github.com/usnistgov>`_

`NIST MetQual program <https://www.nist.gov/programs-projects/metabolomics-quality-assurance-and-quality-control-materials-metqual-program>`_

`NIST Chemical Informatics Research Group <https://www.nist.gov/mml/csd/chemical-informatics-research-group>`_

`NIST home page <http://nist.gov>`_


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
