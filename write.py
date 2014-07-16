"""This module creates the specified inputfiles and writes simulation results to files """

import shutil
from automate import evaluate_odb


def writeInputfile(martensite_amount, austenite_grain, austenite_grains, martensite_grains, \
                   laminate, pbc, geometry_filename, material_jobdata_filename, C_ave=0):
    """ creates an inputfile according to the specified parameters """

    #
    # specify the name of the created inputfile
    inputFile_name = 'Inputfile_' + str(martensite_amount) + '_' + str(austenite_grain[0]) + '_' + \
                     str(laminate) + '.inp'
    # the first section of the created inputfile is the used mesh from the specified
    # external file. Copy this external file and rename it to the specified inputfile name
    shutil.copy2(geometry_filename, inputFile_name)
    # append section definitions and assignments to inputfile according to previous results
    with open(inputFile_name, 'a') as ifile:
        # ----- write sections -----
        for iGrain in austenite_grains:
            # only write currently transforming grain once
            if austenite_grain == iGrain: continue
            string = '*Solid Section, elset=transig_' + str(iGrain[0]) + \
                     ', orientation=Ori_' + str(iGrain[0]) + ', material=austenite\n'
            ifile.write(string)
        # write sections for already transformed grains
        for iGrain in martensite_grains:
            # iGrain = [ grainNr, laminate ]
            string = '*Solid Section, elset=transig_' + str(iGrain[0]) + \
                     ', orientation=Ori_' + str(iGrain[0]) + ', material=laminate' + \
                     str(iGrain[1]) + '\n'
            ifile.write(string)
        # write section of additional grain transforming in this increment
        string = '*Solid Section, elset=transig_' + str(austenite_grain[0]) + \
                 ', orientation=Ori_' + str(austenite_grain[0]) + ', material=laminate' + \
                 str(laminate) + '\n'
        ifile.write(string)
        # ---- write material and jobdata -----
        if pbc == False:
            # write section for self consistent matrix, an orientation is
            # needed because self consistent isotropic properties are given as
            # averaged anisotropic tensor
            string = '*Solid Section, elset=matrix, orientation=Ori_1,' + \
                     'material=selfconsistentIsotropic\n'
            ifile.write(string)
            ifile.write('*Material, name=selfconsistentIsotropic\n*Elastic, type=ANISOTROPIC\n')
            # the specification due to abaqus is first and second line 8
            # and third line 4 entries, see keyword *elastic, type=anisotropic
            entry = 0
            for i in C_ave:
                entry = entry + 1
                ifile.write(i + '\t,')
                if entry == 8:
                    ifile.write('\n')
                    entry = 0
            ifile.write('\n')
        # finally write laminate and job information
        with open(material_jobdata_filename, 'r') as material_jobData:
            for line in material_jobData:
                ifile.write(line)


def writeSaves(martensiteAmount, evaluationData, chem_drivForce, foundGrain):
    """write data of Energy minimizing-configuration in the two files 'save_gain'
    containing grain specific data and 'save_model' containing model specific data.
    Also save the information of all other transformations in an increment as a reference.
    foundGrain = [0 - delta_allEner, 1 - GrainNr, 2 - GrainLaminate, 3 - GrainVol,
                  4 - dragEner_spec, 5 - delta_totalStrain_spec,
                  6 - stress_drivingForce, 7 - total_strainEner_cell ] """
    #
    # write data from all runs of the actual increment
    with open('saves/allruns_' + str(martensiteAmount), 'w') as save_all:
        save_all.write('grainNr\tlaminateNr\tgrainVol\tdragEner_spec\t\t' +
                       'delta_totStrainEner_spec\tdelta_allEnergies\n')
        for dat in evaluationData:
            save_all.write(str(dat[1]) + '\t' + str(dat[2]) + '\t\t' + str(dat[3]) + '\t' +
                           str(dat[4]) + '\t' + str(dat[5]) + '\t\t' + str(dat[0]) + '\n')
            # write grain data:
    with open('saves/save_grain', 'a') as save_grain:
        save_grain.write(str(fg[1]) + '\t' + str(fg[2]) + '\t\t' + str(fg[3]) + '\t' + str(fg[4]) +
                         '\t' + str(fg[5]) + '\t\t' + str(fg[0]) + '\t' + str(- hD + chem_drivForce) + '\n')
    # gather and write model data:
    odbname = 'Outputfile_' + str(martensiteAmount) + '_' + str(fg[1]) + \
              '_' + str(fg[2]) + '.odb'
    md = automate.evaluate_odb(odbname, var=0)  # md ... modeldata
    md_6 = fg[3] * fg[5]  # this is delta_total_strainEner

    # md = [0 - tot_strainEner, 1 - tot_aveSener, 2 - ivol_aust, 3 - ivol_mart,
    # 4 - aveSener_aust,  5 - aveSener_mart ]
    with open('saves/save_model', 'a') as save_model:
        save_model.write(str(md[0]) + '\t' + str(md_6) + '\t' + str(md[1]) + '\t' + str(md[2]) + '\t' +
                         str(md[3]) + '\t' + str(md[4]) + '\t' + str(md[5]) + '\n')








