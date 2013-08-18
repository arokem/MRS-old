import os
import tempfile

import numpy as np
import numpy.testing as npt
import matplotlib
matplotlib.use('agg')

import nitime as nt
import nibabel as nib

import MRS
import MRS.utils as ut
import MRS.analysis as ana

test_path = os.path.join(MRS.__path__[0], 'tests')
file_name = os.path.join(test_path, 'pure_gaba_P64024.nii.gz')

def test_separate_signals():
    """
    Test separation of signals for water-suppressed and other signals
    """
    data = np.transpose(nib.load(file_name).get_data(), [1,2,3,4,5,0]).squeeze()
    w_data, w_supp_data = ana.separate_signals(data)
    # Very simple sanity checks
    npt.assert_equal(w_data.shape[-1], data.shape[-1])
    npt.assert_equal(w_supp_data.shape[-1], data.shape[-1])
    npt.assert_array_equal(data[1], w_data[0])
    
def test_coil_combine():
    """
    Test combining of information from different coils
    """
    data = np.transpose(nib.load(file_name).get_data(), [1,2,3,4,5,0]).squeeze()
    w_data, w_supp_data = ana.coil_combine(data)
    # Make sure that the time-dimension is still correct: 
    npt.assert_equal(w_data.shape[-1], data.shape[-1])
    npt.assert_equal(w_supp_data.shape[-1], data.shape[-1])

    # Check that the phase for the first data point is approximately the
    # same and approximately 0 for all the water-channels:
    npt.assert_array_almost_equal(np.angle(w_data)[:,:,0],
                                  np.zeros_like(w_data[:,:,0]),
                                  decimal=1)  # We're not being awfully strict
                                              # about it.
        
def test_get_spectra():
    """
    Test the function that does the spectral analysis
    """
    data = np.transpose(nib.load(file_name).get_data(), [1,2,3,4,5,0]).squeeze()
    w_data, w_supp_data = ana.coil_combine(data)

    # XXX Just basic smoke-testing for now:
    f_nonw, nonw_sig1 = ana.get_spectra(nt.TimeSeries(w_supp_data,
                                       sampling_rate=5000))

    f_nonw, nonw_sig2 = ana.get_spectra(nt.TimeSeries(w_supp_data,
                                       sampling_rate=5000),
                                       line_broadening=5)

    f_nonw, nonw_sig3 = ana.get_spectra(nt.TimeSeries(w_supp_data,
                                       sampling_rate=5000),
                                       line_broadening=5,
                                       zerofill=1000)
    
def test_mrs_analyze():
    """
    Test the command line utility
    """
    mrs_path = MRS.__path__[0]
    out_name = tempfile.NamedTemporaryFile().name
    # Check that it runs through:
    cmd = '%s/../bin/mrs-analyze.py %s '%(mrs_path, file_name)
    npt.assert_equal(os.system(cmd),0)
    # XXX We might want to analyze the output file here...


