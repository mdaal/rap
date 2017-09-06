import numpy.ma as ma
import numpy as np

def _points_removed(initial, final):
	''' Compute and return the number of point removed from inital due to a cut. 
	return this number and the number of points in final'''
	try:
		initial_number = initial.size - initial.mask.sum()
	except:
		initial_number = initial.size

	try: 	
		final_number = final.size - final.mask.sum()
	except:
		final_number = final.size
	return (initial_number - final_number), final_number


def _angle(z, deg = 0, return_offset = False):
	''' 
	IF Z IS A VECTOR, THEN ANGLE IS SHIFTED WRT FIRST ELEMENT!!!!

	If z is a masked array. angle(z) returns the angle of the elements of z
	within the branch [0,360] instead of [-180, 180], which is the branch used
	in np.angle(). The mask of angle(z) is set to be the mask of the input, z.

	If z is not a masked array, then angle(z) is the same as np.angle except 
	that range is [0,360] instead of [-180, 180]

	If z is a vector, then an angle shift is added to z  so the z[0] is 0 degrees
	If z is a number, then dont shift angle'''
	a = np.angle(z, deg = deg)
	
	try:
		offset = a[0] #if a is not a vector, then a[0] will throw an error
		a = a - offset  
	except:
		pass
	p = np.where(a<=0,1,0)
	n = 2
	units = n*np.pi if deg == 0 else n*180
	try:
		a = ma.array(a + p*units, mask =z.mask) 
	except:
		a = a + p*units #if z is not a masked array do this
	
	if return_offset:
		return a, offset
	else:
		return a

