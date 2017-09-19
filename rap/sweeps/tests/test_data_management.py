from ..data_management.utils import  _download_data

# def test_read_scandata_from_file():
# 	#SWP.load_scandata('rap/sweeps/tests/49a_survey.mat')
# 	mat, filename_or_path = _read_scandata_from_file('rap/examples/52a_survey.mat')
# 	assert isinstance(mat, dict) and ('ScanData' in mat.keys())


def test_download_data():
	#SWP.load_scandata('rap/sweeps/tests/49a_survey.mat')
	mat, filename_or_path = _download_data('http://cosmology.berkeley.edu/~miguel/ResearchWebPages/49a_survey.mat') 
	print(mat)
	assert isinstance(mat, dict) and ('ScanData' in mat.keys())

def test_load_scandata(SWP):
	SWP.load_scandata('rap/examples/52a_survey.mat')
	print(SWP.metadata.Run)
	assert True