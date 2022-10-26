import numpy as np
from pyzdcf import pyzdcf

test_data = '../test_data/test_lc2/'
results = '../test_data/test_lc2/results_py/'

# Length = 3500 days, R_BLR = 34.1 days
lc1 = 'lcx_3500'
lc2 = 'lcy_3500'

sparse = 'auto'
eps = 0.01

# Create input dictionary
def create_input(params):
    keys = ["autocf", "prefix", "uniform_sampling", "omit_zero_lags", "minpts",
            "num_MC", "lc1_name", "lc2_name"]
    return dict(zip(keys,params))

# Define tests
def checks(dcf_py, dcf_for):
    if dcf_for.ndim == 1:
        dcf_for = dcf_for.reshape(1,-1)
        assert dcf_py.shape == dcf_for.shape
    else:
        mid = int(dcf_py.shape[0]*0.5)
        assert dcf_py.shape == dcf_for.shape
        assert np.isclose(dcf_for[:20,:], dcf_py[:20,:], rtol=eps).all()
        assert np.isclose(dcf_for[-20:,:], dcf_py[-20:,:], rtol=eps).all()
        assert np.isclose(dcf_for[mid,:], dcf_py[mid,:], rtol=eps).all()
        
# Perform tests with different input parameters
def test_ccf_uniT_zeroT_MC1():
    dcf_for = np.loadtxt(test_data+'ccfTT1.dcf')
    params=(False,'ccfTT1_test',True,True,'0','1',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)

def test_ccf_uniF_zeroT_MC1():
    dcf_for = np.loadtxt(test_data+'ccfFT1.dcf')
    params=(False,'ccfFT1_test',False,True,'0','1',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    

    
def test_ccf_uniT_zeroF_MC1():
    dcf_for = np.loadtxt(test_data+'ccfTF1.dcf')
    params=(False,'ccfTF1_test',True,False,'0','1',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)


def test_ccf_uniF_zeroF_MC1():
    dcf_for = np.loadtxt(test_data+'ccfFF1.dcf')
    params=(False,'ccfFF1_test',False,False,'0','1',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
    
