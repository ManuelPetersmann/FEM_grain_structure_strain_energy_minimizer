''' It also holds the parameters for the semianalytical calculation of
all IEs constituting a part of the energy barrier in the transformation '''

import numpy as np # used for a faster function of tensor rotation using outer products
import math	   
#my modules
import mathutils 

def eigenstrains():
''' initializes the matrices representing the transformation strains'''
	stress_driving_forces = []
	# transformation strain components
	[ v, u, w ] = [ 0.0381, 0.0234 0.1186 ]
	# e1 corresponds to laminate 1, e2 to laminate 2 etc.
	e1 = [ u, -v,  u, 0,-w, 0 ]
	e2 = [ u,  u, -v, w,  0, 0 ]
	e3 = [ u, -v,  u, 0,  w, 0 ]
	e4 = [ u,  u,  -v,-w, 0, 0 ]
	e5 = [ v,  u,  u, 0,  0, -w ]
	e6 = [ -v, u,  u, 0,  0, w ]
	e = [ e1, e2, e3, e4, e5, e6]
	transforming_strains =  []
	for i in e:
		transforming_strains.append( mathutils.fillmatrix( i ) )
	return transforming_strains


def calc_draggingForces( volume_grain ):
	''' Given the grainvolume of a certain grain in [nm^3] the specific dragging forces for
	the grain are calculated. Therefore the grain is assumed to be a sphere with equal volume
	as the grain. '''
	#
	# ----- < P A R A M E T E R S > -----
	sigma_twin = 0.014E-9 # [nJ/nm^2] specific twin interfaceenergy
	poissons_ratio_austenite = 0.4 
	alpha_twin = 0.856 * (1+poissons_ratio_austenite) / (8*math.pi) 
	shear_distortion = 0.1677 # 2*math.sqrt( e12**2 + e23**2 ) 
	youngs_modulus_austenite = 72.E-9  
	shear_modulus_austenite = youngs_modulus_austenite / (2*(1+poissons_ratio_austenite)) 
	inclusion_energy = 2.42E-10 # 1.177E-10 
	fit_C =  inclusion_energy / (alpha_twin*2*1.5*shear_distortion**2 * shear_modulus_austenite) 
	# 1.5nm = d, 2d = D, the thickness of the 1nm model is not written explicitly
	Fc = 5.8E-12 # [nJ/nm^3] work of friction, dissipated energy
	delta_surfaceEnergy = 0.1E-9 # [nJ/nm^2]; values from 0.1 ... 0.4 are reasonable
	#
	# At first the equal sphere's radius and surface is calculated
	diameter, surf = equal_lengths(volume_grain)
	# next the optimal twin width minimizing the dragging force is calculated
	d_opt = math.sqrt( sigma_twin / (12 * fit_C * alpha_twin * shear_modulus_austenite * \
	                   shear_distortion**2) ) * math.sqrt(diameter)
	#
	# The dragging forces follow as:
	total_interfaceEnergy = ( diameter / (6*d_opt) ) * sigma_twin 
	twin_surfaceEnergy = fit_C*alpha_twin*shear_modulus_austenite*shear_distortion**2 *2 *d_opt
	dragging_forces = (total_interfaceEnergy + delta_surfaceEnergy + twin_surfaceEnergy) \
	                   * (surf / volume_grain) + Fc
	return dragging_forces


def equal_lengths(volume):
	''' given the volume of a grain calculates the diameter
	and surface of a sphere with equal volume '''
	diameter = ( (6*volume) / math.pi ) ** (1./3.)
	surf = diameter**2 * math.pi
	return diameter, surf


def selfconsistent_matrix( oris, graindata, Vinner ):
	'''This function averages anisotropic elastic constants ( refering to local coordinate 
	systems respectively) to global isotropic elastic constants considering phase fractions.
	The nearly isotropic elastic constants are used as the matrix material property. In 
	micromechanics this is commonly called "self consistence scheme" '''
	#
	# define isotropic elastic constants for austenite
	E_aust = 70e-9	
	poissons_ratio_aust = 0.4
	prefactor_austenite = E_aust / (( 1 + poissons_ratio_aust)*( 1 - 2*poissons_ratio_aust))
	# Assemble isotropic elastic constants of austenite as a fourth order tensor
	# note that the the following lines hold the 21 independent entries of the elastic tensor
	# respectively. The order is the same as in the abaqus inputfile. 
	# first line in inputfile
	A1111 = prefactor_austenite * (1 - poissons_ratio_aust)
	A1122 = A2211 = prefactor_austenite * poissons_ratio_aust
	A2222 = prefactor_austenite * (1 - poissons_ratio_aust)
	A1133 = A3311 = prefactor_austenite * poissons_ratio_aust
	A2233 = A3322 = prefactor_austenite * poissons_ratio_aust
	A3333 = prefactor_austenite * (1 - poissons_ratio_aust)
	A1112 = A1211 = A1121 = A2111 = 0.
	A2212 = A1222 = A2221 = A2122 = 0. 
	# second line in inputfile
	A3312 = A1233 = A3321 = A2133 =  0.
	A1212 = A2112 = A1221 = A2121 = prefactor_austenite * ((1 - 2*poissons_ratio_aust) / 2)
	A1113 = A1311 = A1131 = A3111 = 0.
	A2213 = A1322 = A2231 = A3122 = 0.
	A3313 = A1333 = A3331 = A3133 = 0.
	A1213 = A1312 = A2113 = A1231 = A3112 = A1321 = A2131 = A3121 = 0. 
	A1313 = A3113 = A1331 = A3131 = prefactor_austenite * ((1 - 2*poissons_ratio_aust) / 2)
	A1123 = A2311 = A1132 = A3211 = 0. 
	# thirA line in inputfile
	A2223 = A2322 = A2232 = A3222 = 0.
	A3323 = A2333 = A3332 = A3233 = 0.
	A1223 = A2312 = A2123 = A1232 = A3212 = A2321 = A2132 = A3221 = 0.
	A1323 = A2313 = A3123 = A1332 = A3213 = A2331 = A3132 = A3231 = 0.
	A2323 = A3223 = A2332 = A3232 = prefactor_austenite * ((1 - 2*poissons_ratio_aust) / 2)
	# anisotropic elastic constants of martensite given in the basis of the tetragonal unit cell 
	M1111 = 2.54e-07
	M1122 = M2211 = 1.04e-07
	M2222 = 1.8e-07
	M1133 = M3311 = 1.36e-07
	M2233 = M3322 = 1.51e-07
	M3333 = 2.48e-07
	M1112 = M1211 = M1121 = M2111 = 0. 
	M2212 = M1222 = M2221 = M2122 = 0. 
	M3312 = M1233 = M3321 = M2133 = 0.
	M1212 = M2112 = M1221 = M2121 = 9.1e-8
	M1113 = M1311 = M1131 = M3111 = 0.
	M2213 = M1322 = M2231 = M3122 = 0.
	M3313 = M1333 = M3331 = M3133 = 0.
	M1213 = M1312 = M2113 = M1231 = M3112 = M1321 = M2131 = M3121 = -3.e-9
	M1313 = M3113 = M1331 = M3131 = 9.3e-08
	M1123 = M2311 = M1132 = M3211 = 2.1e-08
	M2223 = M2322 = M2232 = M3222 = 0.
	M3323 = M2333 = M3332 = M3233 = -6.e-09
	M1223 = M2312 = M2123 = M1232 = M3212 = M2321 = M2132 = M3221 = 0.
	M1323 = M2313 = M3123 = M1332 = M3213 = M2331 = M3132 = M3231 = 0.
	M2323 = M3223 = M2332 = M3232 = 5.e-09
	#
	Ca =[ [ [ [A1111, A1112, A1113], [A1121, A1122, A1123], [A1131, A1132, A1133] ],
	         [ [A1211, A1212, A1213], [A1221, A1222, A1223], [A1231, A1232, A1233] ],
	         [ [A1311, A1312, A1313], [A1321, A1322, A1323], [A1331, A1332, A1333] ] ],
	       [ [ [A2111, A2112, A2113], [A2121, A2122, A2123], [A2131, A2132, A2133] ],
	         [ [A2211, A2212, A2213], [A2221, A2222, A2223], [A2231, A2232, A2233] ],
	         [ [A2311, A2312, A2313], [A2321, A2322, A2323], [A2331, A2332, A2333] ] ],
	       [ [ [A3111, A3112, A3113], [A3121, A3122, A3123], [A3131, A3132, A3133] ],
	         [ [A3211, A3212, A3213], [A3221, A3222, A3223], [A3231, A3232, A3233] ],
	         [ [A3311, A3312, A3313], [A3321, A3322, A3323], [A3331, A3332, A3333] ] ] ]
	#
	Cm =[ [ [ [M1111, M1112, M1113], [M1121, M1122, M1123], [M1131, M1132, M1133] ],
	         [ [M1211, M1212, M1213], [M1221, M1222, M1223], [M1231, M1232, M1233] ],
	         [ [M1311, M1312, M1313], [M1321, M1322, M1323], [M1331, M1332, M1333] ] ],
	       [ [ [M2111, M2112, M2113], [M2121, M2122, M2123], [M2131, M2132, M2133] ],
	         [ [M2211, M2212, M2213], [M2221, M2222, M2223], [M2231, M2232, M2233] ],
	         [ [M2311, M2312, M2313], [M2321, M2322, M2323], [M2331, M2332, M2333] ] ],
	       [ [ [M3111, M3112, M3113], [M3121, M3122, M3123], [M3131, M3132, M3133] ],
	         [ [M3211, M3212, M3213], [M3221, M3222, M3223], [M3231, M3232, M3233] ],
	         [ [M3311, M3312, M3313], [M3321, M3322, M3323], [M3331, M3332, M3333] ] ] ]
	# create numpy arrays
	Ca = np.array(Ca)
	Cm = np.array(Cm)
	C_ave = np.zeros((3,3,3,3))
	for i in range(len(graindata)):
		# get rotationmatrix between local and global coordinate system. 
		rot_euler = np.array( mathutils.calc_rotmatrix_euler( oris[i][0], oris[i][1] ) )
		# distinguish between austenite and martensite
		if graindata[i][2] == 'AUSTENITE': C = Ca
		else: C = Cm
		Vi = graindata[i][1]
		# calculate  C_averaged = sum_x[ (Vx/Vinner) * ( Rim Rjn Rkp Rlq Cmnpq ) ]
		# where Cmnpq can be C_a or C_m	
		C_ave = np.add( C_ave, np.multiply( (Vi/Vinner), rotateElasticTensor(C,rot_euler) ) )
	C_selfconsistent_voigt = voigt_notation( C_ave )	
	# convert float entries to strings with five decimals in order to write to inputfile
	for i in range( len(C_selfconsistent_voigt) ):
		C_selfconsistent_voigt[i] = '{0:.5e}'.format(float( C_selfconsistent_voigt[i]) ) 
		# positional argument 0 in python 2.x required.
	return C_selfconsistent_voigt
