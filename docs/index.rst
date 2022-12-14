.. pyzdcf documentation master file, created by
   sphinx-quickstart on Fri Oct 14 12:07:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyZDCF documentation!
================================

**pyZDCF** is a Python module that emulates a widely used Fortran
program called ZDCF (Z-transformed Discrete Correlation Function,
`Alexander 1997 <https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract>`__).
It is used for robust estimation of cross-correlation function of sparse
and unevenly sampled astronomical time-series. This Python
implementation also introduces sparse matrices in order to significantly
reduce RAM usage when running the code on large time-series (> 3000
points).

pyZDCF is based on the original Fortran code fully developed by 
Prof. Tal Alexander from Weizmann Institute of Science, Israel 
(see :ref:`Acknowledgements <ackn>` and :ref:`References <refer>` for details and further reading).

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   pyZDCF/pyZDCF.rst
   Examples/examples.ipynb
   API-docs/API.rst

   

   
   
   
   
   

