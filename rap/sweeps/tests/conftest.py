# go to directory containing rap pckage
# issue command 'pytest -v' for verbose or -s for silent or -vv -s for showing all print() statements

from ..sweep import sweep

import pytest


@pytest.fixture(scope="session") # the returned fixture value will be shared for all tests needing it
def SWP(request): # create a sweep object
	print('Creating a global sweep object')
	swp = sweep()
	def fin():
		print('Deleting the global sweep object')
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