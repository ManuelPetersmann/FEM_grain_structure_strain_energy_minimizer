"""This script evaluates the grain-laminate pair which minimizes the total free energy
density of the specified RVE upon transformation. If previous increments are fully 
calculated the script automatically starts from the last state that is saved in file 
'continuing_data'. The script needs N = (1 + maxGrainNr)* 6 * (maxGrainNr / 2) calculations
to find the energy minimizing state of a full transformation. A preselection of more likely
states based on previous results can be done reducing N significantly. At the beginning 
the script-parameters and the used textfiles have to be specified."""

# python modules
import cPickle as pickle  # Phython module to save intermediate results conveniently
import shutil  # high level file operations like copying
import glob  # Unix style pathname pattern expansion
import os  # miscellaneous operating system interfaces
import time  # module for time access
# my modules
import write
import automate
import material


# -----< SPECIFY SCRIPT-PARAMETERS >-------------------------------------------------------#
pbc = True  # use the specified periodic boundary equations or a self-consistent matrix
preselection = True  # define if presection is used
# define at which increments (number of transformed grains) all possibillities for the 
# preselection are calculated 
selected_steps = [1, 15, 30, 45, 60, 75, 85, 95, 101, 107, 113, 117, 121, 124]
# set the timeout after which a single calculation is killed if it has not finished
timeout = 1200
# initialize array of numbers that define the transformed material behavior (here laminates)		
laminate_variants = [1, 2, 3, 4, 5, 6]
# choose between periodic boundary conditions for the regular tesselation or the self 
# consistent matrix for the random microstrocture and specify the required files
if pbc == True:
    material_jobData_filename = ' path/to/file'  # holds the equations for the PBCs
    geometry_filename = ' path/to/file'  # fixed mesh and orientations for the PBC model
    total_grain_amount = 128  # total number of equally sized octahedra in the RVE
    grain_volume = 268000.0  # Volume of a sphere of 80nm diameter
else:
    material_jobData_filename = ' path/to/file'  # static inputfile section of random RVE
    geometry_filename = ' path/to/file'  # fixed mesh for the ESCM
    # here the orientations are written explicitly since they are also used for the
    # averaging of the material properties in each increment
    orientation_filename = ' path/to/file'


#-----< DETERMINE LAST STATE of simulation (if) or  PREPARE CALCULATION (else) >----------#
if os.path.isfile('saves/continuing_data'):  # check if there are calculation results
    # open file and load results in the same order they were written
    cont = open('saves/continuing_data', 'rb')
    martensite_amount = pickle.load(cont)
    oris = pickle.load(cont) if pbc == False else    0  # ternary operator assignment
    martensite_grains = pickle.load(cont)
    austenite_grains = pickle.load(cont)
    chemical_drivingForce = pickle.load(cont)
    total_strain_energy_cell_before = pickle.load(cont)
    selected_variants = pickle.load(cont)
    cont.close()

    # get amount of grains for the randomly generated microstructure
    total_grain_amount = len(austenite_grains) + len(martensite_grains)
    # Define name of .odb file containing latest evaluated result
    odbname = 'saves/Outputfile_' + str(martensite_amount - 1) + '_' + \
              str(martensite_grains[-1][0]) + '_' + str(martensite_grains[-1][1]) + '.odb'
#
else:  # create save directory and files
    martensite_grains = selected_variants = []
    martensite_amount = 1
    os.system('mkdir saves')  # create directory where results are saved
    with open('saves/save_grain', 'w') as save_grain:
        save_grain.write('grainNr\tlaminateNr\tgrainVol\tdragEner_spec\t\t' + \
                         'delta_totStrainEner_spec\tdelta_allEnergies\ttransformingEnergy\n')
    with open('saves/save_model', 'w') as save_model:
        save_model.write('tot_strainEner\t\tdelta_tot_strainEner\ttot_aveSener' + \
                         '\t\tivol_aust\t\tivol_mart\t\tave_sener_aust\t\tave_sener_mart\n')
    #
    if pbc == False:
        odbname = exodb_filename
        austenite_grains = automate.get_volumes_and_laminates(exodb_filename)
        # austeniteGrains [grainNr, grainvolume, grainmaterial]
        total_grain_amount = len(austenite_grains)  # get Nr of grains in the ESCM
    else:
        austenite_grains = []
        for i in range(total_grain_amount):
            austenite_grains.append([i + 1, grain_volume, 0])
        # recall that the grain volume is equal for all octahedra


#-----< CALCULATE SELF CONSISTENT MATERIAL PROPERTIES and PRESELECT more likely states >--#
if pbc == False:
    # list [grainnumber,  grainvolume,  grainmaterial]
    # necessary for averaging anisotropic data
    graindata, Vinner = automate.get_volumes_and_laminates.get(odbname)
    #
    # calculate the averaged material properties from the last energy-minimizing state
    C_ave = material.selfconsistent_matrix(oris, graindata, Vinner)
    del graindata, Vinner
#
# Limit number of calculations by preselecting more likely states
if preselection == True:
    # calculate all possible states only in every selected stepwidth
    if martensite_amount in selected_steps:
        preselection = False
    if (martensite_amount - 1) in selected_steps:
        selected_variants = automate.preselect(total_grain_amount, martensite_amount)


    #-----< INPUTFILE CREATION >--------------------------------------------------------------#
# All possible or preselected states of one more transformed grain are evaluated
for austenite_grain in austenite_grains:
    # Every not transformed grain can transform in multiple ways
    for laminate in laminate_variants:
        if ( [austenite_grain[0], laminate] in selected_variants) or preselection == False:
            write.writeInputfile(martensite_amount, austenite_grain, \
                                 austenite_grains, martensite_grains, laminate, pbc, \
                                 geometry_filename, material_jobData_filename, C_ave)



#-----< JOB SUBMISSION of all Jobs that were created >------------------------------------#
automate.submitjobs(austenite_grains, martensite_amount, laminate_variants, preselection, \
                    selected_variants, timeout)
# delay to finish operations on the .odb files so that no *lck files are created
time.sleep(60)


#-----< EVALUATE ALL jobs and SET PARAMETERS for the transformation of the next grain >---#   
if martensite_amount == 1:
    evaluationData = automate.find_minimum_energy(austenite_grains, martensite_amount, \
                                                  laminate_variants, preselection)
    # evaluationData =  [0-delta_G, 1-GrainNr, 2-GrainLaminate, 3-GrainVol,
    #             4-dragEner_spec, 5-delta_totalStrain_spec,   6-total_strainEner_cell]
    chemical_drivingForce = max(evaluationData)[0]  # note that this is a negative value
else:
    evaluationData = automate.find_minimum_energy(austenite_grains, martensite_amount, \
                                                  laminate_variants, preselection, selected_variants, \
                                                  total_strain_energy_cell_before, chemical_drivingForce)
# # the optimum grain is that with minimum delta_G
found_grain = min(abs(evaluationData))
# the max function acts on the first entry which is 'delta_G'
#
# if delta_G reaches a new negative maximum the chemical driving force
# has to be increased for further transformations
if found_grain[0] < 0:  # if delta_G < 0
    chemical_drivingForce = found_grain[0]
#
total_strain_energy_cell_before = found_grain[6]


#-----< WRITE DATA of all runs and energy-minimizing configuration to files >-------------#
write.writeSaves(martensite_amount, evaluationData, chemical_drivingForce, found_grain)
del evaluationData


#-----< MOVE FOUNDGRAIN from austeniteGrains to martensiteGrains >------------------------#
for index, iGrain in enumerate(austenite_grains):
    # remove found grain from austeniteGrains
    if iGrain[0] == found_grain[1]:
        austenite_grains.pop(index)
# add found grain - material pair to martensiteGrains
martensite_grains.append([found_grain[1], found_grain[2]])
#
# remove foundgrain from selected_variants if preselection is used
if preselection == True:
    for index, iGrain in enumerate(selected_variants):
        if iGrain[0] == found_grain[1]:
            selected_variants.pop(index)


#-----< SAVE (PICKLE) EVALUATED NECESSARY VARIABLES for the next increment >--------------#
# it is crucial that the values are 'loaded' in the same order they are 'dumped'
cont = open('saves/continuing_data', 'wb')
pickle.dump(martensite_amount + 1, cont)
if pbc == False: pickle.dump(oris, cont)
pickle.dump(martensite_grains, cont)
pickle.dump(austenite_grains, cont)
pickle.dump(chemical_drivingForce, cont)
pickle.dump(total_strain_energy_cell_before, cont)
pickle.dump(selected_variants, cont)
cont.close()


#-----< SAVE FILES OF FOUNDGRAIND AND DELETE THE REST >-----------------------------------#
savefilenames = glob.glob('*_' + str(martensite_amount) + '_' + str(found_grain[1]) + \
                          '_' + str(found_grain[2]) + '*')  # example '*_17_44_5*'
for i in savefilenames:
    shutil.move(i, 'saves')  # generally: src --> destination, here: i --> saves
#
# delete all other files
os.system('rm *.*')


#-----< CREATE STOPPINGFILE >-------------------------------------------------------------#
# Obviously Abaqus only allows a maximum number of around 1000 jobs in one interactive
# Python # session. As a workaround, the script is written for one increment and is 
# recalled from an external bash script. When the transformation is finished, the file 
# finish_loop file is created # causing the bash-script not to recall the script anymore
if martensite_amount == total_grain_amount:
    with open('saves/finish_loop', 'w') as f:
        pass





