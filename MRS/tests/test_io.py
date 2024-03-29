import os
import numpy as np

import numpy.testing as npt
import scipy.io as sio

import nibabel as nib

import MRS

test_path = os.path.join(MRS.__path__[0], 'tests')

file_name = os.path.join(test_path, 'pure_gaba_P64024.nii.gz')

def test_get_data():
    """
    Test that the data is the same as the one read using Matlab tools.

    """ 
    matlab_data = sio.loadmat(os.path.join(test_path, 'data.mat'),
                              squeeze_me=True)['data']
    mrs_data = np.transpose(nib.load(file_name).get_data(),
                            [1,2,3,4,5,0]).squeeze()

    # There are differences in the data because the numbers are very large, so
    # we check that the differences are minscule in relative terms:
    mat_dat_e8 = matlab_data.squeeze() / 1e8
    mrs_dat_e8 = mrs_data.squeeze() / 1e8
    npt.assert_almost_equal(mrs_dat_e8, mat_dat_e8)
    
