#!/usr/bin/env python

__author__ = "Isidora Jankov"
__credits__ = ["Andjelka Kovacevic", "Dragana Ilic", "Paula Sanchez Saez", "Robert Nikutta"]
__maintainer__ = "Isidora Jankov"
__contact__ = "isidora_jankov@matf.bg.ac.rs"
__license__ = "MIT License"
__version__ = "1.0.0"

# Imports

import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp


# Constants

ENOUGH = 11
ISEED = 123456
RSLUTN = 0.001
EPS = 1e-7


# Module functions

def read_obs(
    fname,
    input_dir="./",
    output_dir="./",
    out_name="light_curve.lc",
    sep=",",
    savelc=False):
    """
    Read the light curve data. The file should contain three columns: time, 
    flux/magnitude and flux/magnitude error. The file should not contain header
    (column names).

    Parameters
    ----------
    fname : str
        Name of the file.
    input_dir : str, optional
        Path to the folder containing the file. Defaults to './'.
    output_dir : str, optional
        Path to the folder where the results will be written. Defaults to './'.
    out_name : str, optional
        Name of the condensed light curve file. Defaults to 'light_curve.lc'.
    sep: str, optional
        Delimiter to use when reading the file. Defaults to ','.
    savelc : bool, optional
        If True, save the condensed light curve, otherwise skip this step.
        Defaults to False.

    Returns
    -------
    lc : pandas.DataFrame
        Light curve data decribed by three parameters: time, flux and flux error.
    """

    # Load data
    lc = pd.read_table(
        input_dir + fname, header=None, names=["t", "flux", "err"], sep=sep
    )
    assert not (
        np.isnan(lc["flux"]).all()
    ), f"NaN values encountered in {input_dir+fname}. Check if the file format is correct."

    lc.sort_values(by="t", ascending=True, inplace=True)

    # Average the observations with identical times
    lc = lc.groupby("t").mean().reset_index()

    # Save .lc file (condensed)
    if savelc == True:
        lc.to_csv(output_dir + out_name, index=False)
        print(f"{out_name} written (contains {len(lc.t)} points)")

    return lc


def simerr(flux, err):
    """
    Generating a "true" signal by substracting gaussian noise from observations.
    The error err is the ABSOLUTE error in flux.
    """
    n = len(err)
    mc = flux - (err * rndnrm(n))
    return mc


def rndnrm(n):
    """
    Generating an array (size: n) of standard Gaussian deviates using the
    Box-Muller method.
    """
    n2 = int((n + 1) / 2)
    x1 = np.random.uniform(size=n2)
    x2 = np.random.uniform(size=n2)
    x2 = x2 * 2 * np.pi
    C = np.sqrt(-2 * (np.log(x1)))
    x1 = C * np.cos(x2)
    x2 = C * np.sin(x2)
    rnd = np.ndarray(n)
    rnd[:n2] = x1
    rnd[n2:] = x2[: n - n2]

    return rnd


def fishs(r, n):
    """
    Fisher's small sample approximation for s(z) (Kendall + Stuart Vol. 1 p.391)
    """
    fish = (
        1.0
        / (n - 1)
        * (
            1
            + (4 - r**2) / 2.0 / (n - 1)
            + (22 - (6 * r**2) - 3 * (r**4)) / 6 / (n - 1) ** 2
        )
    )
    try:
        fish = np.sqrt(np.maximum(np.zeros(len(fish)), fish))
    except (TypeError, ValueError):
        fish = np.sqrt(max(0, fish))

    return fish


def fishe(r, n):
    """
    Fisher's small sample approximation for E(z) (Kendall + Stuart Vol. 1 p.391)
    """
    fish = np.log((1 + r) / (1 - r)) / 2.0 + r / 2 / (n - 1) * (
        1
        + (5 + r**2) / 4 / (n - 1)
        + (11 + 2 * r**2 + 3 * r**4) / 8 / (n - 1) ** 2
    )
    return fish


def tlag_pts(a, b, verbose=True):
    """
    Generate all possible point pairs {a_i, b_j} and order them by their 
    associated time lags.
    """

    a_n = len(a.t)
    b_n = len(b.t)
    a_t = a["t"].values
    b_t = b["t"].values

    # Allocating the dynamic work areas
    # -9999 is an arbitrary number different than zero (zero could be mistaken
    # for a proper index value)
    waidx = np.full(a_n * b_n, -9999, dtype=np.int32)
    wbidx = np.full(a_n * b_n, -9999, dtype=np.int32)
    wtau = np.full(a_n * b_n, -9999, dtype=np.float32)

    first_time = True

    if first_time:
        first_time = False
        if verbose:
            print(
                f"\nBinning with minimum of {Minpts} points "
                f"per bin and resolution of {RSLUTN}\n"
            )

    # Calculate all the time lag points

    if Autocf:

        l1 = 0
        for i in range(a_n):
            tij = b_t[i:] - b_t[i]
            wb = np.arange(i, b_n)
            if (NoZeroLag == True) & (np.any(tij == 0)):
                wb = wb[tij != 0]
                tij = tij[tij != 0]
            l2 = len(tij) + l1
            wtau[l1:l2] = tij
            waidx[l1:l2] = i
            wbidx[l1:l2] = wb
            l1 = l2

    else:

        l1 = 0
        for i in range(a_n):
            tij = b_t[:] - a_t[i]
            wb = np.arange(b_n)
            if (NoZeroLag == True) & (np.any(tij == 0)):
                wb = wb[tij != 0]
                tij = tij[tij != 0]
            l2 = len(tij) + l1
            wtau[l1:l2] = tij
            waidx[l1:l2] = i
            wbidx[l1:l2] = wb
            l1 = l2

    npp = len(wtau[wtau != -9999])  # Number of all pairs

    wtau = wtau[:npp]
    waidx = waidx[:npp]
    wbidx = wbidx[:npp]

    # Sort index according to increasing time lag
    idx = np.argsort(wtau, kind="heapsort")

    return wtau, waidx, wbidx, idx, npp


def range_inclusive(start, end, incr):
    """
    Range function which includes both starting and ending values.
    It is implemented for easier translation of the code from the
    Fortran programming language.
    """
    if incr == 1:
        return range(start, end + 1, incr)
    elif incr == -1:
        return range(start, end - 1, incr)


def alcbin(a, b, sparse="auto", verbose=True):
    """
    Generate and order all possible point pairs by their associated time-lag. 
    "The ordered list is then divided, bin by bin, into bins of Minpts pairs. 
    In the process of adding pairs to a bin, a new pair whose a or b points 
    have previously appeared in that bin, is discarded. The allocation of pairs
    to bins starts at the lag where the pairs are densest. For the 
    auto-correlation function, this means that the binning proceeds from τ = 0 
    up to τmax. For the CCF, the binning proceeds from the median τ up to τmax
    and then from the median τ down to τmin." (Alexander, 2013).
    
    Sparse matrices are added as a new feature in this function. They are
    utilized when allocating work areas in case when there is > 3000 points per
    light curve to economize RAM usage.
    """

    # Calculating all time-lag points and their index
    wtau, waidx, wbidx, idx, npp = tlag_pts(a, b, verbose=verbose)

    # Calculating the tolerance level for lags to be considered the same
    tij = wtau[idx[npp - 1]] - wtau[idx[0]]
    tolrnc = tij * RSLUTN

    # The dcf arrays are allocated here with the maximal number of bins
    # possible for the requested Minpts. The actual number used is
    # subsequently determined on the fly and stored in dcf_nbins variable.
    # It is typically much smaller than mbins.

    mbins = int((npp + 1) / Minpts) + 1

    dcf_t = np.full(mbins, -9999, dtype=float)
    dcf_sigtm = np.full(mbins, -9999, dtype=float)
    dcf_sigtp = np.full(mbins, -9999, dtype=float)
    dcf_inbin = np.full(mbins, -9999, dtype=np.int32)

    # Automatically set flag for choosing data type of work area matrix
    if sparse == "auto":
        if max(a_n, b_n) > 3000:  # sparse matrix in case of long light curves
            sparse = True
        else:
            sparse = False  # otherwise, use numpy arrays

    # Optimize data type usage to save memory
    if max(a_n, b_n) > 32000:
        int_dtype = np.int32
    else:
        int_dtype = np.int16

    # Allocate work areas (used later by clcdcf) in order to economize
    # memory usage.

    if sparse == True:
        # use sparse work area matrix to save memory
        wa_i = sp.lil_matrix((maxpts, mbins), dtype=int_dtype)
        wb_i = sp.lil_matrix((maxpts, mbins), dtype=int_dtype)
        if verbose:
            print("using sparse matrices...")

    elif sparse == False:
        # use numpy arrays (more memory used, but faster indexing)
        wa_i = np.zeros((maxpts, mbins), dtype=int_dtype)
        wb_i = np.zeros((maxpts, mbins), dtype=int_dtype)

    nbins = -1

    # If binned CCF: binning from median time-lag upwards and backwards!

    if Autocf or UniSample:
        pfr = 0
        pmax = npp - 1
        incr = 1

    else:
        pfr = int((npp / 2) - 1)
        pmax = 0
        incr = -1

    # Looping on bins

    # bin_loop

    while True:
        inb = 0
        tij = wtau[idx[pfr]]

        nbins += 1

        # This shouldn't happen...
        if nbins > mbins:
            sys.exit("alcbin: nbins > mbins (this shouldn't happen...)")

        dcf_t[nbins] = 0.0

        # Initialize used point flag vectors
        a_used = np.zeros(a_n, dtype=bool)
        b_used = np.zeros(b_n, dtype=bool)

        # Collect points into bins that contain at least "ENOUGH" points,
        # but do not break points with the same lag (up to the tolerance)
        # into separate bins

        # tau_loop

        continue_bin_loop = False

        for i in range_inclusive(pfr, pmax, incr):

            p = idx[i]

            # Check whether bin is full
            if (
                (abs(wtau[p] - tij) > tolrnc)
                & ((inb >= Minpts) | UniSample)  # tij = previous lag
            ) | (i == pmax):
                # Bin is full: Calculating tau and its std
                # (before proceeding to the next bin)
                dcf_inbin[nbins] = inb
                dcf_t[nbins] = dcf_t[nbins] / inb

                if UniSample:
                    dcf_sigtm[nbins] = 0.0
                    dcf_sigtp[nbins] = 0.0
                    pfr = i
                    # If not ENOUGH points in bin, ignore it
                    if inb < Minpts:
                        nbins = nbins - 1
                    if pfr != pmax:
                        continue_bin_loop = True
                        break
                else:
                    # If the last point is alone in its bin, we ignore it
                    # (to avoid subsequent divisions by zero) and get
                    # immediately out of tau loop.
                    if inb <= 1:
                        nbins = nbins - 1
                        break

                    # Finding the 0.3413 (+/-1sig) time lags above and below the bin mean time lag
                    if incr == 1:
                        plo = pfr
                        phi = i - 1
                    else:
                        plo = i + 1
                        phi = pfr

                    # midpnt_loop
                    for p in range_inclusive(plo, phi, 1):
                        if wtau[idx[p]] >= dcf_t[nbins]:
                            j = p  # The point closest to the bin's mean
                            break

                    p = max(j - round(float(j - plo) * 0.3413 * 2), plo)
                    dcf_sigtm[nbins] = max(dcf_t[nbins] - wtau[idx[p]], 0.0)
                    p = min(j + round(float(phi - j) * 0.3413 * 2), phi)
                    dcf_sigtp[nbins] = max(wtau[idx[p]] - dcf_t[nbins], 0.0)
                    pfr = i

                    if pfr != pmax:
                        continue_bin_loop = True
                        break

                # If no more points - get out of tau loop.
                break

            # Adding another point to the bin...

            # Avoiding correlated pairs
            if (not a_used[waidx[p]]) and (not b_used[wbidx[p]]):
                inb += 1
                # This shouldn't happen...
                if inb > maxpts:
                    sys.exit("ALCBIN: inb > maxpts = {}.".format(maxpts))

                a_used[waidx[p]] = True
                b_used[wbidx[p]] = True
                tij = wtau[p]
                dcf_t[nbins] = dcf_t[nbins] + tij
                wa_i[inb - 1, nbins] = waidx[p]
                wb_i[inb - 1, nbins] = wbidx[p]
                # assert nbins <= 1000, "nbins > 1000!"

        if continue_bin_loop:
            continue

        # Binning is finished
        if not (Autocf or UniSample or (incr == 1)):
            # Now, go back and bin the other half lag axis
            pfr = int(npp / 2) + 1 - 1
            pmax = npp - 1
            incr = 1
            nnegtv = nbins
            continue
        break

    # If CCF (and NOT uniform sampling): Sort the bins into increasing
    # chronological order: The nnegtv negative bins are at the
    # beginning but at reverse order.

    if not (Autocf or UniSample):
        for i in range_inclusive(0, int(nnegtv / 2), 1):
            j = nnegtv - i
            # swapping
            dcf_inbin[i], dcf_inbin[j] = dcf_inbin[j], dcf_inbin[i]
            dcf_t[i], dcf_t[j] = dcf_t[j], dcf_t[i]
            dcf_sigtp[i], dcf_sigtp[j] = dcf_sigtp[j], dcf_sigtp[i]
            dcf_sigtm[i], dcf_sigtm[j] = dcf_sigtm[j], dcf_sigtm[i]
            for p in range(max(dcf_inbin[i], dcf_inbin[j])):
                wa_i[p, i], wa_i[p, j] = wa_i[p, j], wa_i[p, i]
                wb_i[p, i], wb_i[p, j] = wb_i[p, j], wb_i[p, i]
                # assert i <= 1000, "i > 1000!"
                # assert j <= 1000, "j > 1000!"

    dcf_nbins = nbins + 1
    if dcf_nbins == 0:
        sys.exit("alcbin: nbins = 0")

    # Num. of bins is calculated, so we can cut these arrays to save memory...
    wa_i = wa_i[:, :dcf_nbins]
    wb_i = wb_i[:, :dcf_nbins]
    dcf_t = dcf_t[:dcf_nbins]
    dcf_sigtm = dcf_sigtm[:dcf_nbins]
    dcf_sigtp = dcf_sigtp[:dcf_nbins]
    dcf_inbin = dcf_inbin[:dcf_nbins]

    if sparse:
        # convert back to numpy arrays
        wa_i = wa_i.todense()
        wb_i = wb_i.todense()

    # Pack the dcf related values and work areas into tuples for code clarity
    dcf_vars = (dcf_t, dcf_sigtm, dcf_sigtp, dcf_inbin, dcf_nbins)
    work_areas = (wa_i, wb_i)

    return dcf_vars, work_areas, mbins


def dcf_pairs(dcf_inbin, dcf_nbins):
    """
    Calculate the number of inter-dependent pairs.
    """
    dcf_used = np.sum(dcf_inbin[:dcf_nbins])

    if Autocf:
        if NoZeroLag:
            dcf_unused = (a_n * (a_n - 1) / 2) - dcf_used
        else:
            dcf_unused = ((a_n**2) / 2) - dcf_used
    else:
        dcf_unused = (a_n * b_n) - dcf_used

    return dcf_unused


def clcdcf(a, b, dcf_vars, work_areas, MC=False):
    """
    Calculating the discrete correlation function. POSITIVE lag values mean
    b lags after a.
    """
    # Unpack input variables
    dcf_t, dcf_sigtm, dcf_sigtp, dcf_inbin, dcf_nbins = dcf_vars
    wa_i, wb_i = work_areas

    # Intialize dcf r and its errors
    dcf_r = np.zeros(dcf_nbins)
    dcf_sigrm = np.zeros(dcf_nbins)
    dcf_sigrp = np.zeros(dcf_nbins)

    if MC:
        # If another Monte Carlo simulation, use the "true signal" values
        a_flux = a["MC"].values
        b_flux = b["MC"].values
    else:
        # If not a Monte Carlo simulation, use obesrvational flux values
        a_flux = a["flux"].values
        b_flux = b["flux"].values

    # After allocating pairs to bins: calculating the dcf
    for ibin in range(dcf_nbins):
        # Collecting the points of the bin
        n = dcf_inbin[ibin]  # Assumed >1
        if n < 2:
            continue  # This shouldn't happen...

        wa_x = a_flux[wa_i[:n, ibin]]
        wb_x = b_flux[wb_i[:n, ibin]]

        expa = np.sum(wa_x[:n]) / n
        expb = np.sum(wb_x[:n]) / n
        vara = np.sum(wa_x[:n] ** 2)
        varb = np.sum(wb_x[:n] ** 2)

        # Dividing by (n-1) for an unbiased estimator of the
        # correlation coefficient cf Barlow / Statistics, p. 80

        vara = (vara - n * expa**2) / (n - 1)
        varb = (varb - n * expb**2) / (n - 1)
        vnorm = vara * varb

        if vnorm <= 0:
            # Pathological case: normalization factor <= 0
            dcf_r[ibin] = 0
        else:
            expbin = np.sum((wa_x[:n] - expa) * (wb_x[:n] - expb))
            expbin = expbin / np.sqrt(vnorm) / (n - 1)
            # Making sure -1 < r < 1
            try:
                assert (expbin < 1 - EPS) & (
                    expbin > -1 + EPS
                ), "dcf_r is out of bounds of (-1,1)"
            except AssertionError:
                if expbin > 1 - EPS:
                    expbin = 1 - EPS
                elif expbin < -1 + EPS:
                    expbin = -1 + EPS

            dcf_r[ibin] = expbin

        # Calculating the +/- 1 Sigma limits from Fisher's z

        # NOTE: This error estimation is by "bootstrapping": fishe & fishs give
        # the true E(z) and S(z) when the TRUE correlation coefficient is
        # given. We are using the empirical r itself, similarily to the
        # common poissonian estimate of n +/- sqrt(n)

        # z = np.log((1+expbin)/(1-expbin))/2  # never used..
        sigz = fishs(expbin, n)
        expz = fishe(expbin, n)
        dcf_sigrm[ibin] = expbin - np.tanh(expz - sigz)
        dcf_sigrp[ibin] = np.tanh(expz + sigz) - expbin

    dcf_taus = (dcf_t, dcf_sigtm, dcf_sigtp)
    dcf_rs = (dcf_r, dcf_sigrm, dcf_sigrp)
    dcf_other = (dcf_inbin, dcf_nbins)

    return dcf_taus, dcf_rs, dcf_other


def check_user_input(usr_input):
    try:
        # Convert it into integer
        val = int(usr_input)
        return val
    except ValueError:
        try:
            # Convert it into float
            val = float(usr_input)
            return val
        except ValueError:
            print("Invalid input. Please enter a number.")


def user_input(interactive=True, verbose=True, parameters={}):
    """
    Returns input parameters for `pyzdcf` by asking the user to provide them 
    interactively or allow to set them manually.

    Parameters
    ----------
    interactive : bool, optional
        If True, prompt the user to provide parameters interactively.
        Otherwise, use parameters stored in `parameters` keyword argument.
        Defaults to True.

    verbose : bool, optional
        Print additional information while running the calculations. Defaults
        to True.

    params : dict, optional
        A dictionary containing all the required parameters (keys are parameter
        names). Defaults to a dictionary with None values for all parameters.
        
        Required parameters key names: `autocf`, `prefix`, `uniform_sampling`,
        `omit_zero_lags`, `minpts`, `num_MC`, `lc1_name`. The key `lc2_name` 
        needs to be provided only if you choose to perform cross-correlation by
        setting `autocf` to `False`. Check function returns for the meaning of
        these parameters. Used only if the `interactive` argument is set to 
        False.

    Returns
    -------
    Autocf : bool
        If True, calculate autocorrelation. If False, calculate 
        cross-correlation.
    Prefix : str
        Output file prefix.
    UniSample : bool
        If True, set flag to perform uniform sampling of the light curve.
    NoZeroLag : bool
        If True, omit zero lag points.
    Minpts : int
        Minimal number of points per bin.
    nMC : int
        Number of Monte Carlo simulations for error estimation.
    Name1 : str
        Name of the first light curve file.
    Name2 : str
        Name of the second light curve file.
    """

    if verbose:

        print("\npyZDCF begins:\n")

    if interactive:
        inpi = input("Auto-correlation or cross-correlation? (1/2): ")

        if inpi == "1":
            Autocf = True
        elif inpi == "2":
            Autocf = False

        Prefix = input("Enter output files prefix: ")

        inpc = input("Uniform sampling of light curve? (y/n): ")

        if (inpc == "y") | (inpc == "Y"):
            UniSample = True
        elif (inpc == "n") | (inpc == "N"):
            UniSample = False

        Minpts = input("Enter minimal number of points per bin (0 for default): ")
        Minpts = check_user_input(Minpts)
        if Minpts <= 0:
            Minpts = ENOUGH

        inpz = input("Omit zero-lag points? (y/n): ")
        if (inpz == "y") | (inpz == "Y"):
            NoZeroLag = True
        elif (inpz == "n") | (inpz == "N"):
            NoZeroLag = False

        nMC = input("How many Monte Carlo runs for error estimation? ")
        nMC = check_user_input(nMC)
        if nMC <= 1:
            nMC = 0

        if Autocf:
            Name1 = input("Enter name of 1st light curve file: ")
            Name2 = Name1

        else:
            Name1 = input("Enter name of 1st light curve file: ")
            Name2 = input("Enter name of 2nd light curve file: ")

    else:
        # Manual input
        input_keys = ['autocf', 'prefix', 'uniform_sampling', 'omit_zero_lags',
                'minpts', 'lc1_name']
        
        # Check if all required keys are provided by the user
        for i in input_keys:
            if i not in parameters.keys():
                raise KeyError(f"Missing key in input dictionary 'parameters': {i}")  
        
        if parameters['autocf'] == False:
            if 'lc2_name' not in parameters.keys():
                raise KeyError("Missing key in input dictionary 'parameters': lc2_name")
            
        Autocf = parameters['autocf']
        Prefix = parameters['prefix']
        UniSample = parameters['uniform_sampling']
        NoZeroLag = parameters['omit_zero_lags']
        Minpts = parameters['minpts']
        nMC = parameters['num_MC']
        Name1 = parameters['lc1_name']
        if Autocf:
            Name2 = Name1
        else:
            Name2 = parameters['lc2_name']
            
        # Check if input values are valid
        Minpts = check_user_input(Minpts)
        nMC = check_user_input(nMC)
        if Minpts <= 0:
            Minpts = ENOUGH
        if nMC <= 1:
            nMC = 0

    if verbose:
        print(30 * "=")
        print("pyZDCF PARAMETERS:\n")
        print("Autocorrelation?  ", Autocf)
        print("Uniform sampling? ", UniSample)
        print("Omit zero lags?   ", NoZeroLag)
        print("Minimal # in bin: ", Minpts)
        print("# of Monte Carlo: ", nMC)
        print("Monte Carlo seed: ", ISEED)
        print(30 * "=")

    return Autocf, Prefix, UniSample, NoZeroLag, Minpts, nMC, Name1, Name2


# Main function --> pyZDCF
def pyzdcf(
    input_dir,
    output_dir,
    intr=True,
    sep=",",
    verbose=True,
    sparse="auto",
    savelc=False,
    parameters={"autocf": False,
                "prefix": 'ccf',
                "uniform_sampling" : False,
                "omit_zero_lags" : True,
                "minpts" : 0,
                "num_MC" : 100,
                "lc1_name" : 'lc1_example',
                "lc2_name" : 'lc2_example'}):
    """
    Calculates auto-correlation or cross-correlation function of user-provided
    light curve(s) using the ZDCF (Z-transformed Discrete Correlation Function)
    method (Alexander, 1997).
    
    Parameters
    ----------
    input_dir : str
        Path to the directory containing input light curve(s).

    output_dir : str
        Path to the directory for storing the results.

    intr : bool, optional
        If True, use interactive user input. Otherwise, use a 
        dictionary with input parameters in the 'parameters' keyword argument. 
        Defaults to True.

    sep : str, optional
        Delimiter to use when reading the file(s). Defaults to ','.
    
    verbose : bool, optional
        Print additional information while running the calculations.
        Defaults to True.
        
    sparse : {'auto', True, False}, optional
        If True, use sparse matrices when allocating work areas in 
        order to save memory. If False, use classical numpy arrays (not 
        recommended for light curves containing > 3000 points when running 
        pyZDCF on a 8 GB RAM personal computer). Default value is 'auto' and it
        forces the program to use sparse matrices with light curves of 3000 or 
        more points.

    savelc : bool, optional
        If True, save the condensed light curve, otherwise skip this
        step. Defaults to False.

    parameters : dict, optional
        A dictionary containing all the required parameters (keys are parameter
        names). Defaults to a dictionary with placeholder values, it is up to the
        user to provide correct values for given keys (input parameters). 
        Required parameters key names (types) are: `autocf` (bool), `prefix` (str),
        `uniform_sampling` (bool), `omit_zero_lags` (bool), `minpts` (str/int),
        `num_MC` (str/int), `lc1_name` (str). The key `lc2_name` needs to be provided 
        only if you choose to perform cross-correlation (by setting `autocf` to `False`).
        The `params` argument is used only in manual mode (`intr` argument set to `False`).
        For description of all input parameters, consult the "How to use" section in the docs.

    Returns
    -------
    dcf_df : pandas.DataFrame
        A pandas DataFrame organized in 7 columns displaying the results for
        each time-lag bin. The columns are: time-lag, negative time-lag std, 
        positive time-lag std, zdcf, negative zdcf sampling error, positive 
        zdcf sampling error, number of points per bin. 
    """

    # User inputs are treated as global constants
    global Autocf, UniSample, NoZeroLag, Minpts, nMC

    # Interactive user input
    if intr:
        Autocf, Prefix, UniSample, NoZeroLag, Minpts, nMC, Name1, Name2 = user_input(
            interactive=True, verbose=verbose
        )

    # Manual input
    else:
        Autocf, Prefix, UniSample, NoZeroLag, Minpts, nMC, Name1, Name2 = user_input(
            interactive=False, parameters=parameters, verbose=verbose
        )

    # Constants derived from user input are also global constants
    global a_n, b_n, maxpts

    # Load data
    a = read_obs(
        Name1, input_dir, output_dir, out_name=Prefix + ".lc1", sep=sep, savelc=savelc
    )
    if Autocf:
        b = a.copy()
    else:
        b = read_obs(
            Name2,
            input_dir,
            output_dir,
            out_name=Prefix + ".lc2",
            sep=sep,
            savelc=savelc,
        )

    # Initialize global constants infered from user input
    a_n = len(a.t)
    b_n = len(b.t)
    maxpts = min(a_n, b_n)

    # Additional info: maxpts is the maximal number of points in a bin.
    # In principle maxpts = npp ( the number of all pairs: a huge number).
    # In practice, for non uniform sampling it is not supposed to be much
    # larger than Minpts (sometimes there are a few extra points per
    # bin if they are closer to the bin's boundary than the resolution
    # tolerance).  For uniform sampling, it is bound by min(a_n,b_n)
    # To be on the safe side, min(a_n,b_n) is assumed for both cases.

    # Estimating the effects of the measurement errors by Monte Carlo simulations

    np.random.seed(ISEED)

    if nMC > 1:
        # Calculate the ZDCF with Monte Carlo errors
        for i in range(nMC):
            a["MC"] = simerr(a.flux, a.err)
            b["MC"] = simerr(b.flux, b.err)

            if i == 0:  # First MC iteration
                dcf_vars, work_areas, mbins = alcbin(
                    a, b, sparse=sparse, verbose=verbose
                )
                dcf_taus, dcf_rs, dcf_other = clcdcf(
                    a, b, dcf_vars, work_areas, MC=True
                )
                dcf_avz = np.zeros(len(dcf_rs[0]))
                if verbose:
                    dcf_unused = dcf_pairs(*dcf_other)
                    print(
                        f"{dcf_other[1]} bins actually used, {dcf_unused} inter-dependent pairs discarded."
                    )
                dcf_avz = dcf_avz + np.log((1 + dcf_rs[0]) / (1 - dcf_rs[0])) / 2
                continue

            dcf_taus, dcf_rs, dcf_other = clcdcf(a, b, dcf_vars, work_areas, MC=True)
            dcf_r, dcf_sigrm, dcf_sigrp = dcf_rs
            dcf_inbin = dcf_other[0]

            # Making sure -1 < dcf < 1 in all bins
            try:
                assert (np.max(dcf_r) < 1 - EPS) & (
                    np.min(dcf_r) > -1 + EPS
                ), "dcf_r is out of bounds of (-1,1)"
            except AssertionError:
                cond1 = dcf_r < -1 + EPS
                cond2 = dcf_r > 1 - EPS
                dcf_r[cond1] = -1 + EPS
                dcf_r[cond2] = 1 - EPS

            # The summing and averaging is done in z-space.
            dcf_avz = dcf_avz + np.log((1 + dcf_r) / (1 - dcf_r)) / 2

        dcf_avz = dcf_avz / nMC
        dcf_r = np.tanh(dcf_avz)

        # Making sure -1 < dcf < 1 in all bins
        try:
            assert (np.max(dcf_r) < 1 - EPS) & (
                np.min(dcf_r) > -1 + EPS
            ), "dcf_r is out of bounds of (-1,1)"
        except AssertionError:
            cond1 = dcf_r < -1 + EPS
            cond2 = dcf_r > 1 - EPS
            dcf_r[cond1] = -1 + EPS
            dcf_r[cond2] = 1 - EPS

        dcf_avz = np.log((1 + dcf_r) / (1 - dcf_r)) / 2
        dcf_expz = fishe(dcf_r, dcf_inbin)
        dcf_sigz = fishs(dcf_r, dcf_inbin)
        dcf_sigrm = dcf_r - np.tanh(dcf_expz - dcf_sigz)
        dcf_sigrp = np.tanh(dcf_expz + dcf_sigz) - dcf_r

        dcf_t, dcf_sigtm, dcf_sigtp = dcf_taus

    else:
        # calculate the ZDCF w/o Monte Carlo errors.
        dcf_vars, work_areas, mbins = alcbin(a, b, sparse=sparse, verbose=verbose)
        dcf_taus, dcf_rs, dcf_other = clcdcf(a, b, dcf_vars, work_areas)
        dcf_t, dcf_sigtm, dcf_sigtp = dcf_taus
        dcf_r, dcf_sigrm, dcf_sigrp = dcf_rs
        dcf_inbin = dcf_other[0]
        if verbose:
            dcf_unused = dcf_pairs(*dcf_other)
            print(
                f"{dcf_other[1]} bins actually used, {int(dcf_unused)} inter-dependent pairs discarded."
            )

    # Save dcf table
    d = np.column_stack(
        (dcf_t, dcf_sigtm, dcf_sigtp, dcf_r, dcf_sigrm, dcf_sigrp, dcf_inbin)
    )
    lim = dcf_other[1]
    d = d[:lim, :]
    np.savetxt(
        output_dir + Prefix + ".dcf",
        d,
        fmt=["%11.3e", "%11.3e", "%11.3e", "%11.3e", "%11.3e", "%11.3e", "%5.1i"],
    )

    dcf_df = pd.DataFrame(
        data=d,
        columns=[
            "tau",
            "-sig(tau)",
            "+sig(tau)",
            "dcf",
            "-err(dcf)",
            "+err(dcf)",
            "#bin",
        ],
    )
    if verbose:
        print("\n")
        print(Prefix + ".dcf written...")
        print("\npyZDCF ended.")
        print("")

    return dcf_df


if __name__ == "__main__":
    # Use interactive mode when running pyzdcf as a script
    input_dir = input("Enter the path to the directory containing input data: ")
    output_dir = input("Enter the path to the directory for storing the results: ")
    dcf_df = pyzdcf(input_dir, output_dir, intr=True, sep=",")
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(dcf_df)
    # to get more pretty table you can use `tabulate` module (pip install tabulate)
    #print(tabulate(dcf_df, headers="keys", tablefmt="psql"))

    # Manual mode (change parameters to your liking and then run this script)
    # input_dir = './'
    # output_dir = './'
    # lc1 = 'name1'
    # lc2 = 'name2'
    # params = {"autocf": False, 
    #           "prefix": 'pref', 
    #           "uniform_sampling" : False, 
    #           "omit_zero_lags" : True,
    #           "minpts" : 0,
    #           "num_MC" : 100,
    #           "lc1_name" : lc1,
    #           "lc2_name" : lc2}
    # #params = (False,'ccfFT100',False,True,0,100,lc1,lc2)
    # dcf_df = pyzdcf(input_dir,
    #                 output_dir,
    #                 intr=True,
    #                 parameters = params,
    #                 sep=',',
    #                 verbose=True,
    #                 sparse='auto')
