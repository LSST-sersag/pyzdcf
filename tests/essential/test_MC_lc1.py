import numpy as np
from pyzdcf import pyzdcf

test_data = '../test_data/test_lc1/'
results = '../test_data/test_lc1/results_py/'

# Length = 1000 days, R_BLR = 34.1 days
lc1 = 'lcx_1000'
lc2 = 'lcy_1000'

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
        assert dcf_py.shape == dcf_for.shape
        assert np.isclose(dcf_for[:20,0:3], dcf_py[:20,0:3], rtol=eps).all()
        assert np.isclose(dcf_for[-20:,0:3], dcf_py[-20:,0:3], rtol=eps).all()
        
# Perform tests with different input parameters
def test_acf_uniT_zeroT_MC3():
    dcf_for = np.loadtxt(test_data+'acfTT3.dcf')
    params=(True,'acfTT3_test',True,True,'0','3',lc1)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
def test_acf_uniF_zeroT_MC3():
    dcf_for = np.loadtxt(test_data+'acfFT3.dcf')
    params=(True,'acfFT3_test',False,True,'0','3',lc1)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
def test_acf_uniT_zeroF_MC3():
    dcf_for = np.loadtxt(test_data+'acfTF3.dcf')
    params=(True,'acfTF3_test',True,False,'0','3',lc1)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
def test_acf_uniF_zeroF_MC3():
    dcf_for = np.loadtxt(test_data+'acfFF3.dcf')
    params=(True,'acfFF3_test',False, False,'0','3',lc1)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)


def test_ccf_uniT_zeroT_MC3():
    dcf_for = np.loadtxt(test_data+'ccfTT3.dcf')
    params=(False,'ccfTT3_test',True,True,'0','3',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
def test_ccf_uniF_zeroT_MC3():
    dcf_for = np.loadtxt(test_data+'ccfFT3.dcf')
    params=(False,'ccfFT3_test',False,True,'0','3',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
def test_ccf_uniT_zeroF_MC3():
    dcf_for = np.loadtxt(test_data+'ccfTF3.dcf')
    params=(False,'ccfTF3_test',True,False,'0','3',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
    
def test_ccf_uniF_zeroF_MC3():
    dcf_for = np.loadtxt(test_data+'ccfFF3.dcf')
    params=(False,'ccfFF3_test',False, False,'0','3',lc1,lc2)
    d = create_input(params)
    dcf_py = pyzdcf(input_dir=test_data, output_dir=results, 
                    intr=False, parameters = d, sparse=sparse)
    
    dcf_py = dcf_py.to_numpy()
    checks(dcf_py, dcf_for)
