from ..data_management.utils import  _download_data


filename_a = 'rap/examples/legacy_gui_ps'
filename_b = 'rap/examples/Coupling6_rescale2_full_23pH.s2p'
filename_c = 'rap/examples/52a_survey.mat'
url_a = 'http://cosmology.berkeley.edu/~miguel/ResearchWebPages/49a_survey.mat'
url_b = 'http://cosmology.berkeley.edu/~miguel/private/KIPS_Data/49a_survey.mat' #password protected

def test_load_legacy_sweep_gui_data(SWP):
    SWP.load_legacy_sweep_gui_data(filename_a)
    assert SWP.Sweep_Array['Frequencies_Syn'][0].sum() > 0

def test_load_touchstone(SWP):
    SWP.load_touchstone(filename_b)
    assert True

def test_download_data():
    mat, filename_or_path = _download_data(url_a)
    # print(mat)
    assert isinstance(mat, dict) and ('ScanData' in mat.keys())

def test_load_scandata(SWP):
    SWP.load_scandata(filename_c, Verbose = False)
    # print(SWP.metadata.Run)
    assert True


# def test_load_scandata_pw(SWP):
#     SWP.load_scandata(url_b, Verbose = False, username = '', password = '')
#     # print(SWP.metadata.Run)
#     assert True
