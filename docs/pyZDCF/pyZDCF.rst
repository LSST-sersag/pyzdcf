pyZDCF
======

Description
-----------

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


Motivation
----------

Development of pyZDCF module was motivated by the long and successful usage of 
the original ZDCF Fortran code in the analysis of light curves of active galactic
nuclei by our research group (see `Kovacevic et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014AdSpR..54.1414K/abstract>`__,
`Shapovalova et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4790S/abstract>`__,
and reference therein). One of the science cases we investigate is photometric reverberation
mapping in the context of Legacy Survey of Space and Time (LSST) survey strategies (see `Jankov et
al. 2022 <https://ui.adsabs.harvard.edu/abs/2022AN....34310090J/abstract>`__). However, this module is **general** and is meant to be used for cross-correlation of spectroscopic or photometric light curves, same as the original Fortran version.

Installation
------------

pyZDCF can be installed using pip:

.. code:: sh

   pip install pyzdcf

Dependencies
~~~~~~~~~~~~

   ::

      python = ">=3.8,<3.11"
      numpy = ">=1.16.5,<1.23.0"
      pandas = "^1.3"
      scipy = "^1.7.3"


How to use
----------

Input files
~~~~~~~~~~~

This code requires user-provided plain text files as input. CSV files
are accepted by default, but you can use any other delimited file, as
long as you provide the ``sep`` keyword argument when calling ``pyzdcf``
function. The input light curve file should be in 3 columns: time
(ordered), flux/magnitude and absolute error on flux/magnitude. Make
sure to exclude the header (column names) from the input files.

First few lines of the example input file accepted by default (CSV):

::

   0.0,0.9594479339474323,0.0019188958678948648
   1.0,0.9588196871078336,0.0019176393742156672
   2.0,0.9637198686651904,0.0019274397373303808
   3.0,0.9622807967282166,0.0019245615934564328
   
.. note::
   pyZDCF is tested only with input files having whole numbers (integers) for 
   time column. If you have decimal numbers (e.g., you have a light curve with 
   several measurments in the same night expressed as fractions of a day 
   instead of minutes), just convert them into a time format with integer 
   values (e.g., minutes instead of days). On the other hand, you could round 
   the values from the same day (e.g. 5.6 --> 5, 5.8 --> 5, etc.) and the 
   algorithm will take in the information and average the flux for that day. 
   
Input parameters
~~~~~~~~~~~~~~~~

If you use interactive mode (``intr = True``), then pyZDCF will ask you to
enter all input parametars interactively, similarly to original ZDCF interface.
There is also a manual mode (``intr = False``) where you can provide input
parameters using a dictionary and passing it to ``parameters`` keyword argument.

Available input parameters (keys in the ``parameters`` dictionary) are:

- ``autocf`` - if ``True``, perform auto-correlation, otherwise do the cross-correlation.
- ``prefix`` - provide a name for the output .dcf file.
- ``uniform_sampling`` - if ``True``, set flag to perform uniform sampling of the light curve.
- ``omit_zero_lags`` - if True, omit zero lag points.
- ``minpts`` - minimal number of points per bin.
- ``num_MC`` - number of Monte Carlo simulations for error estimation.
- ``lc1_name`` - Name of the first light curve file
- ``lc2_name`` - Name of the second light curve file (required only if we do cross-correlation)

For more information on the correct syntax, see "Running the code" subsection and the "API Documentation" section.

Output
~~~~~~

The return value of the ``pyzdcf`` function is a ``pandas.DataFrame`` object
displaying the results in 7 columns:

.. csv-table::
   :header: index,tau,-sig(tau),+sig(tau),dcf,-err(dcf),+err(dcf),#bin
   :widths: 5,10,12,12,12,12,12,5
   
   0,-991,4,0,0.13598,0.361559,0.342224,10
   1,-988,2,0,-0.217733,0.279988,0.301034,13
   2,-985,2,0,-0.0614938,0.266546,0.27135,16
   3,-982,2,0,0.239601,0.237615,0.223317,19
   4,-979,2,0,0.331415,0.208171,0.192523,22


The columns are: time-lag, negative time-lag std, positive time-lag std,
zdcf, negative zdcf sampling error, positive zdcf sampling error, number
of points per bin. For more information on how these values are
calculated see `Alexander 1997 <https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract>`__.

The code will also generate an output .dcf file in a specified folder on your
computer with same 7 columns containing the results. It is allowed to
name these files however you want using the ``prefix`` parameter (see
example in the next subsection).

Optionally, by adding keyword argument ``savelc = True``, ``pyzdcf`` can
create and save light curve files used as input after averaging points
with identical times.

Running the code
~~~~~~~~~~~~~~~~

An example for calculating cross-correlation between two light curves:

.. code:: python

   from pyzdcf import pyzdcf

   input = './input/'           # Path to the input data
   output = './output/'         # Path to the directory for saving the results

   # Light curve names
   lc1 = 'lc_name1'
   lc2 = 'lc_name2'

   # Parameters are passed to the pyZDCF as a dictionary

   params = dict(autocf            =  False, # Autocorrelation (T) or cross-correlation (F)
                 prefix            = 'ccf',  # Output files prefix
                 uniform_sampling  =  False, # Uniform sampling?
                 omit_zero_lags    =  True,  # Omit zero lag points?
                 minpts            =  0,     # Min. num. of points per bin (0 is a flag for default value of 11)
                 num_MC            =  100,   # Num. of Monte Carlo simulations for error estimation
                 lc1_name          =  lc1,   # Name of the first light curve file
                 lc2_name          =  lc2    # Name of the second light curve file (required only if we do CCF)
                )

   # Here we use non-interactive mode (intr=False)
   dcf_df = pyzdcf(input_dir  = input, 
                   output_dir = output, 
                   intr       = False, 
                   parameters = params, 
                   sep        = ',', 
                   sparse     = 'auto', 
                   verbose    = True)

   # To run the program in interactive mode (like the original Fortran code):
   dcf_df = pyzdcf(input_dir  = input, 
                   output_dir = output, 
                   intr       = True, 
                   sep        = ',', 
                   sparse     = 'auto', 
                   verbose    = True
                   )

-  For more examples see `example
   notebook <https://github.com/LSST-sersag/pyzdcf/blob/main/notebooks/examples.ipynb>`__.
-  Additionally, you can also check out code description of the original
   Fortran version because the majority of input parameters and all
   output files are the same as in pyZDCF. You can download the fortran
   source code
   `here <https://www.weizmann.ac.il/particle/tal/research-activities/software>`__.
..
   -  To see pyZDCF application in the context of a science case involving
      photometric reverberation mapping of active galactic nuclei see `this
      tutorial <https://github.com/LSST-sersag/dle-private/blob/main/AstroDataLab_NB/%5BDRAFT%5D%20PhotoRM%20NB%20PyZDCF.ipynb>`__.


Features
--------

-  **Sparse matrix implementation** for reduced
   RAM usage when working with long light curves (>3000 points). 
      - The main benefit is that we can now run these demanding calculations
       on our own personal computers (8 GB of RAM is enough for light curves 
       containing up to 15000 points), making the usage of this algorithm more
       convinient than ever. You can turn this on/off by specifying ``sparse`` 
       keyword argument to ``True`` or ``False``. Default value is ``'auto'``, 
       where sparse marices are utilized when there are more than 3000 points
       per light curve. Note that by reducing RAM usage, we pay in increased 
       program running time.
-  **Interactive mode**: program specifically asks the user to provide
   necessary parameters (similar to original Fortran version);
-  **Manual mode**: user can provide all parameters in one dictionary.
-  **Fixed bugs** from original ZDCF (v2.3) written in Fortran 95.

.. note::
   The module was tested (i.e., compared with original ZDCF v2.3 output, with fixed bugs) for various parameter combinations on a set of 100 AGN light curve candidates (*g* and *r* bands). 
   The list of object ids and coordinates was taken from a combined catalogue of known AGNs (`Sánchez-Sáez et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021AJ....162..206S/abstract>`__).



License
-------

Distributed under the MIT License.

Contact
-------

   | **Isidora Jankov (main)** - ijankov@proton.me
   | Andjelka Kovačević - andjelka@matf.bg.ac.rs
   | Dragana Ilić - dilic@matf.bg.ac.rs

You are welcome to write to us:
	- if there are any problems running the code on your system;
	- to give suggestions for code improvements.

If you want to report a bug, please open an Issue on GitHub.

Project Link: https://github.com/LSST-sersag/pyzdcf

Citation
--------

If you use pyZDCF for scientific work leading to a publication, 
please consider acknowledging it using the following DOI: `<https://doi.org/10.5281/zenodo.7253034>`__


Acknowledgments
---------------
.. _ackn:

-  The pyZDCF module is based on the original Fortran code developed by
   Prof. Tal Alexander (Weizmann Institute of Science, Israel). You can download original Fortran version from professor’s 
   `page <https://www.weizmann.ac.il/particle/tal/research-activities/software>`__.
-  For theoretical details regarding the ZDCF algorithm see this publication: `Alexander 1997 <https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract>`__.
-  Huge thanks to my closest collegues and mentors Dr. Andjelka
   Kovačević and Dr. Dragana Ilić, as well as to Dr. Paula Sánchez-Sáez
   and Dr. Robert Nikutta for invaluable input during the development
   and testing of this python module.  
-  Many thanks to Prof. Eli Waxman, Amir Bar On and former students of
   Prof. Tal Alexander for their kind assistance regarding the
   development of pyZDCF module and for its acknowledgment as part of the
   legacy behind late Prof. Alexander.


References
----------
.. _refer:

-  `Alexander, T. 1997, in: Astronomical Time Series, eds. D. Maoz, A.
   Sternberg, & E. M. Leibowitz, Vol. 218, Springer, Is AGN Variability
   Correlated with Other AGN Properties? ZDCF Analysis of Small Samples
   of Sparse Light
   Curves <https://ui.adsabs.harvard.edu/abs/1997ASSL..218..163A/abstract>`__
-  `Jankov, I.; Kovačević A. B.; Ilić, D.; et al. 2022, Astronomische
   Nachrichten, 343,
   e210090 <https://ui.adsabs.harvard.edu/abs/2022AN....34310090J/abstract>`__
-  `Kovačević, A.; Popović, L. Č.; Shapovalova, A. I.; et al. 2014,
   Advances in Space Research, 54,
   1414-1428 <https://ui.adsabs.harvard.edu/abs/2014AdSpR..54.1414K/abstract>`__
-  `Sánchez-Sáez, P.; Lira, H.; Martí, L.; et al. 2021,
   The Astronomical Journal, 162,
   id.206 <https://ui.adsabs.harvard.edu/abs/2021AJ....162..206S/abstract>`__
-  `Shapovalova, A. I.; Popović, L. Č.; Afanasiev, V. L.; et al. 2019,
   MNRAS, 485,
   4790-4803 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4790S/abstract>`__
   

   
   
   
   
   

