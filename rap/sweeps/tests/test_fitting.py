
from ..sweep import sweep

import pytest


filename_a = 'rap/examples/Coupling6_rescale2_full_23pH.s2p'
filename_b = 'rap/examples/S5_truegapdistance.s2p'

@pytest.fixture(scope="module", params=[filename_a, filename_b])
def SWP_Fitting(request): # create a sweep object
	
	swp = sweep()
	swp.load_touchstone(request.param)
	print('Create SWP_Fitting object')
	yield swp
	print('Teardown SWP_Fitting object')
	#return swp


def test_circle_and_phase_fit(SWP_Fitting):
	SWP_Fitting.pick_loop(0)
	SWP_Fitting.circle_fit(Show_Plot = False)
	SWP_Fitting.phase_fit(Show_Plot = False, Verbose = False)

# def test_pick_loop(SWP_Fitting):
# 	SWP_Fitting.pick_loop(0)
# 	assert True