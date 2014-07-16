''' This module automates the calculation of all generated inputfiles as well as the 
evaluation of the Outputdatabases. Also preselection parameters are specified.'''

import os
import time
import sys # to exit the subprocess after the abaqus standard job has finished
from multiprocessing import Process # manages subprocesses calculating the abaqus jobs
import psutil # library for retrieving information on running processes 
# my modules
import mathutils


def preselect( totalGrainAmount, martensiteAmount ):
	''' This function reduces the number of calculations carried out in a increment
	by preselecting more likely states known from a previous increment. The here 
	defined parameters reduce the number of total calculations for this simulation
	by a factor of 6 while the results are the same. Note that especially, low 
	fractions at early increments lower the number of total calculations because the
	number of not transformed grains and all possibilites are related multiplicatively.'''
	
	# here the fraction of all possible transformations is declared
	if martensiteAmount < int(totalGrainAmount*.7): 
		calc_fraction = (1./7)	
	if martensiteAmount>=int(totalGrainAmount*.7) and martensiteAmount<int(totalGrainAmount*.9):
		calc_fraction = (1./5)
	if martensiteAmount>=int(totalGrainAmount*.9) and martensiteAmount<int(totalGrainAmount*.94):
		calc_fraction = (1./2)
	if martensiteAmount >= int(totalGrainAmount*.95): 
		calc_fraction = 1
	deltas = []
	# variant_preselection is taken from the last state in which 
	# all variant permutations are known
	allstates = 'saves/allruns_' + str( martensiteAmount - 1 ) 
	with open(allstates,'r') as allstates:
		for index, line in enumerate( allstates ):
			line = line.rstrip('\n')
			if index == 0: continue # ignore the headerline
			else:
				# data = [ delta_ener, grain_nr, laminate_nr ]
				data = [ float(line.split()[5]), line.split()[0], line.split()[1] ]
				deltas.append( data )
				del data
	deltas.sort() # sort ascending
	deltas.reverse() # reverse to get descending sort
	amount = int( len(deltas)* calc_fraction ) 
	deltas = deltas[ 0 : amount ] 
	preselection = []
	for i in deltas:
		preselection.append( [ int(i[1]), int(i[2]) ] ) # [grainNr, laminate] 
	return preselection
	
def submitjobs( austeniteGrains, martensiteAmount, laminateVariants, preselection, \
                selected_variants, timeout ):
	''' handles automatic submission of all inputfiles, created in an increment '''
	#
	wait_array = []
	jobNr = 0
	#
	for austeniteGrain in austeniteGrains:
		for laminate in laminateVariants:
			if ( [austeniteGrain[0], laminate] in selected_variants) or preselection == False:
				#
				jobNr = jobNr + 1				
				inputname = 'Inputfile_'+ str(martensiteAmount)+'_'+ \
				str(austeniteGrain[0])+ '_'+ str(laminate) +'.inp'
				outputname = 'Outputfile_'+ str(martensiteAmount)+ '_'+ \
				str( austeniteGrain[0] )+'_'+ str(laminate)
				#
				pid = os.fork() # start (=fork) subprocesses
				if pid == 0: # 0 is the child process
					# in the child process invoke the standard solver 
		    			os.system( '/opt/abaqus/Commands/abq6123 job=' + outputname + \
						' interactive cpus=2 scratch=/dev/shm input=' + \
						inputname +' mp_mode=threads standard_parallel=all')
					# /dev/shm tmpfs directory
					sys.exit() # close subprocess after the job has finished
				else:	# in the parent process the child process are managed
					process = psutil.Process( pid )
					wait_array.append( process)
				# 					
				if jobNr % 6 == 0 or (preselection == True and jobNr == len(selected_variants) ):		
					for iprocess in wait_array:
						try:
							iprocess.wait( timeout )
							# timeout is the time in seconds the script waits for completion of the job
						except psutil.TimeoutExpired:
							print 'process running after timeout, killing proces...'
							os.system('killall -9 standard.exe')
							time.sleep(30)
					wait_array = [] # empty wait array for new processes
	del inputname, outputname, wait_array, jobNr


def findMinimumEnergy( austeniteGrains, martensiteAmount, laminateVariants, preselection, \ 
                       selected_variants = [], total_strainEner_cell_before = 0, \
                       chemical_drivingForce = 0 ):
	''' Reads totalstrainergy from outputfiles and evaluates the transformation that minimizes
	    the total strain energy density '''
	#
	evaluationData = [] # Define list for calculation results
	#
	for austeniteGrain in austeniteGrains:
		for laminate in laminateVariants:
			if ( [austeniteGrain[0], laminate] in selected_variants) or preselection == False:
				#
				odbname = 'Outputfile_'+ str(martensiteAmount)+'_'+ str( austeniteGrain[0] )+ \
						    '_'+ str(laminate)+'.odb'
				#
				# if the calculation was terminated and a .lck file exist ignore that .odb
				if os.path.isfile( odbname[0:len(odbname)-3] +'lck' ): continue 
				# read the totalstrainenergy of the whole model via the odb file
				# alternatively it could be extracted from the .dat file
				else: total_strainEner_cell = evaluate_odb.(odbname, var=1)
				#
				# Calculate difference of free energy density to previous increment:
				# first calculate specific strain energy of transformed grain
				delta_totalStrain = total_strainEner_cell - total_strainEner_cell_before
				delta_totalStrain_spec = delta_totalStrain / austeniteGrain[1]
				# next calculate specific interface energy barrier of transformed grain
				dragEner_spec = forces.calc_draggingForces( austeniteGrain[1] )
				#
				# In the first run the chemical_drivingForce is determined as the sum 
				# of dragging energies, thus it a negative value
				delta_G = chemical_drivingForce - ( dragEner_spec + delta_totalStrain_spec )
				# CAUTION! : line continuation with - \ - gives +! 1--1 = 2!
				#
				evaluationData.append( [delta_G, austeniteGrain[0], laminate, austeniteGrain[1], \
							dragEner_spec, delta_totalStrain_spec, total_strainEner_cell ] )	
				#
				del total_strainEner_cell, delta_totalStrain, delta_totalStrain_spec, \
					dragEner_spec, delta_G
	return evaluationData 


def evaluate_odb(odbname, var = 0):
	''' this function has three different return values (var = 0, 1, 2)
	per default (0) weighted strain energy densities are returned. var = 1: only the total 
	strain energy of the model is returned. var = 2: The transformation criterion for 
	the LTC is returned based on the specific IE energy barrier and the double dot product
	of the averaged strain 	tensor in a grain and its possible eigenstrains'''
	# create odb singular object
	odb = openOdb(path=odbname)
	# go through the object model
	trans_step = odb.steps['Transformation']
	last_frame = trans_step.frames[-1] # [-1] gives last frame
	
	if var == 1:
		histreg = trans_step.historyRegions['Assembly Assembly-1'] 
		# Assembly Assembly-1 is the default repository key that is generated 
		all_total_energies = histreg.historyOutputs['ALLIE'].data 
		# .data is the key in the dictionary for the Allenergies (total energies) array.
		odb.close()
		return all_total_energies[1][1] # [1][1] is the totalstrainenergy
	
	# --- integration point variables ---
	stresses = last_frame.fieldOutputs['S']	# stress tensor components of integration points 
	seners = last_frame.fieldOutputs['SENER'] # strain energy densities of integration points
	ivols = last_frame.fieldOutputs['IVOL']   # integration point volumes
	'''# ---  alternatively wole element variables can be read ---
        last_frame.fieldOutputs['ESEDEN']    # equivalent to SENER for whole elements
	last_frame.fieldOutputs['ELSE']	     # strain energy for all whole Elements
	last_frame.fieldOutputs['EVOL']      # equivalent to IVOL for whole elements '''
	#
	# initialize variables to be weighted from the integration point level
	total_strainEner = total_aveSener = ivol_total = ivol_aust = ivol_mart = 0
	tot_strainEner_aust = tot_strainEner_mart = 0
	#
	for igrain in odb.sections.values(): 
	# or: for i in odb.rootAssembly.elementSets.keys() -> i = set_name from inputfile
		#
		material = igrain.material
		set_name = igrain.name.lstrip('Section-')
		#
		# 'PART-1-1' is the default name of the first created part if none is specified
		grain = odb.rootAssembly.instances['PART-1-1'].elementSets[set_name]
		# integration point variable subsets
		set_stresses = stresses.getSubset(region = grain, position=INTEGRATION_POINT)
		set_ivol = ivols.getSubset(region = grain, position=INTEGRATION_POINT)
		set_sener = seners.getSubset(region = grain, position=INTEGRATION_POINT)
		# for whole element variables the same can be done without the parameter position=...
		#
		# evaluate the elementset related to the section
		grain_sener_sum = grain_ges_ivol = 0
		for i in range(len(set_stresses.values)): # number of all integration points
			# in the set SENER must be weighted with the element volume since 
			# not elements are of the same size. Same as ELSE 
		    	grain_sener_sum=grain_sener_sum+set_sener.values[i].data*set_ivol.values[i].data
			grain_ges_ivol = grain_ges_ivol + set_ivol.values[i].data
		grain_ave_sener = grain_sener_sum / grain_ges_ivol
		total_strainEner = total_strainEner + grain_sener_sum
		# the volume of the matrix must not be considered for the random RVE cell !
		if 'TRANSIG' in set_name: ivol_total = ivol_total + grain_ges_ivol
		#		
		# additionally calculate the strain energy developement in each phase respectively
		if material == 'AUSTENITE':
			ivol_aust = ivol_aust + grain_ges_ivol
			tot_strainEner_aust = tot_strainEner_aust + grain_ave_sener * grain_ges_ivol
		#
		if 'LAMINATE' in material:
			ivol_mart = ivol_mart + grain_ges_ivol
			tot_strainEner_mart = tot_strainEner_mart + grain_ave_sener * grain_ges_ivol

		if var == 2 and material == 'AUSTENITE':	
			# initialize transformation strains
			transforming_strains = eigenstrains()
			# Define average "effective" stress tensor components
			sig_eff_11 = sig_eff_22 = sig_eff_33 = sig_eff_12 = sig_eff_13 = sig_eff_23 = 0
			#
			for i in range(len(set_stresses.values)):
				grain_ges_ivol = grain_ges_ivol + set_ivol.values[i].data
				#
				sig_eff_11 = sig_eff_11 + set_stresses.values[i].data[0]*set_ivol.values[i].data
	    			sig_eff_22 = sig_eff_22 + set_stresses.values[i].data[1]*set_ivol.values[i].data
	    			sig_eff_33 = sig_eff_33 + set_stresses.values[i].data[2]*set_ivol.values[i].data
	    			sig_eff_12 = sig_eff_12 + set_stresses.values[i].data[3]*set_ivol.values[i].data
	    			sig_eff_13 = sig_eff_13 + set_stresses.values[i].data[4]*set_ivol.values[i].data
	    			sig_eff_23 = sig_eff_23 + set_stresses.values[i].data[5]*set_ivol.values[i].data
	    			#
			sig_eff_11 = sig_eff_11 / grain_ges_ivol 
			sig_eff_22 = sig_eff_22 / grain_ges_ivol
			sig_eff_33 = sig_eff_33 / grain_ges_ivol
			sig_eff_12 = sig_eff_12 / grain_ges_ivol
			sig_eff_13 = sig_eff_13 / grain_ges_ivol
			sig_eff_23 = sig_eff_23 / grain_ges_ivol
			#
			sigma = [sig_eff_11, sig_eff_22, sig_eff_33, sig_eff_12, sig_eff_13, sig_eff_23]
			sigma = mathutils.fillmatrix(sigma)
			#
			spec_grain_drag = forces.calc_draggingForces( grain_ges_ivol )
			for laminate in range(1, 6+1): # = [1,2,3,4,5,6]
				spec_driving_stress=mathutils.doubledotproduct(transforming_strains[laminate-1],sigma)
				check_term = spec_driving_stress - spec_grain_drag
				if check_term > driving_force
					driving_force = check_term
					found_grain = [set_name, laminate, driving_force]
	odb.close()
	#
	# calculate total model values
	total_ivol = ivol_aust + ivol_mart 
	# note that for the random RVE this is only the graincluster without the matrix
	total_aveSener = total_strainEner /  ivol_total
	total_aveSener_mart = total_strainEner_mart / ivol_mart
	if ivol_aust != 0: total_aveSener_aust = total_strainEner_aust / ivol_aust

	if var == 2: 
		return found_grain, total_strainEner, total_aveSener, ivol_austenite, ivol_transformed, \
		       tot_aveSener_aust, tot_aveSener_mart 
	
	return total_strainEner, total_aveSener, ivol_austenite, ivol_transformed, \
	       tot_aveSener_aust, tot_aveSener_mart 






