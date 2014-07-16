''' This module contains all mathematical operations needed for the simulation '''

import math
import izip
import numpy as np # installed along with abaqus. Available after invoking abaqus python

# -----< used matrix operations >---------------------------------------------------------#

def doubledot_product(A, B):
	'''returns the doubledot product of two 3x3 matrices'''
	result = 0	
	for i in range(3):	
		for j in range(3):		
			result = result + A[i][j]*B[i][j]			
	return result

def fillMatrix(A):
	''' creates a symmetric 3x3 matrix from its elements
	A = [11,22,33,12,13,23] '''
	result = [ [A[0],A[3],A[4]],[A[3],A[1],A[5]],[A[4],A[5],A[2]] ]
	return result

# -----< used vector operations >---------------------------------------------------------#

def vecDotProduct(vec1,vec2):
	''' calculates the dot or scalar product of two arbitrary vectors'''
	return sum(x1 * x2 for x1, x2 in izip(vec1, vec2))

def vecNorm(vec):
	''' calculates the Euclidean norm of a vector '''
	return vecDotProduct(vec,vec)**0.5

def angle_between_vectors(vec1,vec2):
	'''calculates the angle between two vectors in radians'''
	return math.acos( vecDotProduct(vec1,vec2) / ( vecNorm(vec1)*vecNorm(vec2) ) )

def vecCrossProduct(vec1, vec2):
	''' calculates the cross product between two vectors '''
	return( [ vec1[1] * vec2[2] - vec1[2] * vec2[1],
	          vec1[2] * vec2[0] - vec1[0] * vec2[2],
	          vec1[0] * vec2[1] - vec1[1] * vec2[0])

# -----[ calculation of rotation matrix relating two Cartesion coordinate systems and ]---#
#      [ rotations of fourth order tensors from one coordinate system to the other ] 

def voigt_notation( C )
	''' this function takes the elastic fourth order tensor and returns its Voigt notation,
	    which is a 6 x 6 matrix '''
	return C_voigt = [ C[0][0][0][0],C[0][0][1][1],C[1][1][1][1],
	                     C[0][0][2][2],C[1][1][2][2],C[2][2][2][2],
	                     C[0][0][0][1],C[1][1][0][1],C[2][2][0][1],
	                     C[0][1][0][1],C[0][0][0][2],C[1][1][0][2],
	                     C[2][2][0][2],C[0][1][0][2],C[0][2][0][2],
	                     C[0][0][1][2],C[1][1][1][2],C[2][2][1][2],
	                     C[0][1][1][2],C[0][2][1][2],C[1][2][1][2] ]



def calc_rotmatrix_euler(a,b):
	'''Calculates the rotation matrix R that takes the Cartesion coordinate system
	defined by the two normal unit vectors a and b to the reference system given by 
	x[1,0,0] y[0,1,0] z[0,0,1]. R is defined using three Eulerian angles whereby
	the x-convention (z,y',z'') is used. '''
	#
	# the third vector of the rotated system is calculated
	c = vecCrossProduct(a,b)
	# the Eulerian angles are calculated 
	alpha, beta, gamma = eulerian_angles(a,b,c)
	# the Rotationsmatrix entrys are given as:
	r11 = - math.sin(alpha)*math.sin(gamma) + math.cos(alpha)*math.cos(beta)*math.cos(gamma)
	r12 = math.cos(alpha)*math.sin(gamma) + math.sin(alpha)*math.cos(beta)*math.cos(gamma)
	r13 = - math.sin(beta)*math.cos(gamma)
	r21 = - math.sin(alpha)*math.cos(gamma) - math.cos(alpha)*math.cos(beta)*math.sin(gamma)
	r22 = math.cos(alpha)*math.cos(gamma) - math.sin(alpha)*math.cos(beta)*math.sin(gamma)
	r23 = math.sin(beta)*math.sin(gamma)
	r31 = math.cos(alpha)*math.sin(beta)
	r32 = math.sin(alpha)*math.sin(beta)
	r33 = math.cos(beta)
	return ((r11, r12, r13),(r21, r22, r23),(r31, r32, r33))

def eulerian_angles(a,b,c):
	''' calculates the eulerangles between two Cartesian coordinate systems with 
	the same origin. The 3 Eulerian angles are dependent from each other, which means
	the order of plane rotations carried out to transform the coordinate system is 
	definite. The angle beta is simply the angle between the z axes of both coordinate
	systems (Z and c here). The angle alpha is the angle between the X axis of the 
	reference coordinate system and the projection of c into the X,Y plane. Finally, 
	gamma is the angle between the b - axis and Y'
	'''
	beta = angle_between_vectors( (0.,0.,1.),c )
	alpha = angle_between_vectors( (1.,0.,0.), (c[0],c[1],0) )
	# for gamma Y has to be rotated first
	Y_rotated = [0,0,0]
	Y_rotated[0] = - math.sin(alpha) * 1
	Y_rotated[1] = math.cos(alpha) * 1
	gamma = angle_between_vectors(Y_rotated, b)
	return alpha, beta, gamma

def rotate_indicial( C, R ):
	''' calculates the rotation of a fourth order tensor C by the rotation defined by R, 
	    where both are given as standard python lists '''
  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          for m in range(3):
            for n in range(3):
              for o in range(3):
                for p in range(3):
                  Crot[i][j][k][l]=R[i][m]*R[j][n]*R[k][o]*R[l][p]*C[m][n][o][p]+Crot[i][j][k][l]
	return Crot

def rotateElasticTensor(C, R):
    ''' rotates a forth order tensor C by a rotation given by the matrix R. 
	see: numpy reference-linear-algebra p653 '''
    RR = np.outer(R, R)
    RRRR = np.outer(RR, RR).reshape(4 * R.shape)
    axes = ((0, 2, 4, 6), (0, 1, 2, 3))
    return np.tensordot(RRRR, C, axes)
