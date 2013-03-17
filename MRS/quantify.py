#!/usr/bin/python

import numpy as np

import sklearn.linear_model as lm
import nitime as nt

import MRS.files as fi
import MRS.analysis as ana
import MRS.utils as ut

def lorentzian(l,w,x):
    """
    Calculate a bell-shaped Lorentzian function centered at w

    Parameters
    ----------

    l : float
        Width parameter.

    w : float
        Location parameter.

    x : ndarray
        Evaluate the function at these values.

    Returns
    -------
    ndarray : The Lorentzian function evaluated at x.
    
    Notes
    -----
    This is the lorentzian function define as equation 1 in _[Soher1996]

    [Soher1996] Soher, BJ, Hurd, RE, Sailasuta, N, Barker, PB (1996)
    Quantitation of automated single voxel proton MRS using cerebral water as
    an internal reference. MRM 36: 335-339
    """

    lor = l/(l**2+(x-w)**2)
    # Normalize it to sum to 1, because we are approximating the integral:
    return lor/np.sum(lor)


def reconstruct(pfile, max_ppm=4.3, min_ppm=-0.8, sampling_rate=5000):
    """
    Parameters
    ----------
    pfile : str
        name of pfile with data for reconstruction

    min_ppm : float (optional)
        Lower limit of x axis in ppm. Default: -0.8

    max_ppm : flat (optional)
        Upper limit of x axis in ppm. Default: 4.3

    sampling_rate: float (optional)
        The sampling rate (in Hz). Default: 5000.0 
    
    Returns
    -------
    out_array : recarray
    """
    # Get data from file: 
    data = fi.get_data(pfile)
    # Use the water unsuppressed data to combine over coils:
    w_data, w_supp_data = ana.coil_combine(data.squeeze())
    # Once we've done that, we only care about the water-suppressed data
    f_nonw, nonw_sig = ana.get_spectra(nt.TimeSeries(w_supp_data,
                                        sampling_rate=sampling_rate))
    # The first echo (off-resonance) is in the first output 
    echo1 = nonw_sig[0]
    # The second output is the difference between off- and on-resonance echos:
    echo2 = nonw_sig[0] - nonw_sig[1]    
    f_ppm = ut.freq_to_ppm(f_nonw)
    idx0 = np.argmin(np.abs(f_ppm - min_ppm))
    idx1 = np.argmin(np.abs(f_ppm - max_ppm))
    idx = slice(idx1, idx0)
    # Convert from Hz to ppm and extract the part you are interested in.
    f_ppm = f_ppm[idx]
    e1 = echo1[:,idx]
    e2 = echo2[:,idx]

    return f_ppm, e1, e2


class Fitter(object):
    """

    """
    def __init__(self, ppm, alpha=0.001, l1_ratio=0.5):
        """

        """
        self.ppm = ppm
        
    
    def design_matrix(self):
        """
        Generate a design matrix

        """
        design_matrix = np.array([lorentzian(0.1, w, self.ppm) for w in
                                  self.ppm]).T

        return design_matrix

    def fit(self, data):
        """
        Fit a model to the data 
        """
        design_matrix = self.design_matrix()
        EN = lm.ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
        EN.fit(design_matrix, spectrum)
        self.coef = EN.coef_

