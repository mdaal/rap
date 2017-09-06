# go to directory containing rap pckage
# issue command 'pytest -v' for verbose or -s for silent

from ..sweep import sweep

import pytest


@pytest.fixture
def SWP(request): # create a sweep object
	print('Creating a sweep object')
	swp = sweep()
	def fin():
		print('Deleting the sweep object')
	request.addfinalizer(fin)
	return swp
	


# from nose import with_setup # optional


# from ..sweep import sweep

# def setup_package():
# 	print("") # this is to get a newline after the dots
# 	print("setup_module before anything in this file")
# 	swp = sweep()
# 	return swp
	

 
# def teardown_package():
# 	print("teardown_module after everything in this file")
# 	#del(swp)