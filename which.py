
#
# will provide the contains of the h5 files associated to heckle code :
# - for each files (fields.h5, species.h5, time.h5 a restarts.h5)
# - the names of each groups & sub-groups
# - the name (& possibly the content) of each attributes
# - the name of each datasets
#
# to make it work, just provide the path of the run directory :
# python which.py ~/my/run/path/
#


import sys
import os
import h5py
import numpy as np

path = sys.argv[1]
fieldFile = os.path.join(path, "fields.h5")
specieFile = os.path.join(path, "species.h5")
timeFile = os.path.join(path, "time.h5")
restartFile = os.path.join(path, "restarts.h5")

if os.path.isdir(path) is True :

    if os.path.isfile(fieldFile) is True :
        # _____ fields.h5 _____
        f = h5py.File(fieldFile)
        print("\n### file : {0} ###".format(fieldFile))

        # attribute
        print("\n  # __ list of attributes :\n")
        attributes = f.attrs.items()
        for attribute in attributes :
            print("  attribute [ {0:>20} ] : {1}".format(attribute[0], attribute[1]))

        # groups for each times
        print("\n  # __ list of groups (time) in fields.h5 :\n")
        times = list(f.keys())
        for time in times :
            print("  {0}".format(time))

        # datasets (they are all the sames in a given group,
        # so let's see the content of the first one
        print("\n  # __ list of associated dataset (for each group) :\n")
        firstTimeGrp = '/'+times[0]
        inTheGroup = f.get(firstTimeGrp)
        for dataSet in inTheGroup.items() :
           print("  {0}".format(dataSet[0]))

        f.close()
    else :
        print("\n### file fields.h5 does not exist")


    if os.path.isfile(specieFile) is True :
        # _____ species.h5 _____
        f = h5py.File(specieFile)
        print("\n### file : {0} ###".format(specieFile))

        # attribute
        print("\n  # __ list of attributes :\n")
        attributes = f.attrs.items()
        for attribute in attributes :
            print("  attribute [ {0:>20} ] : {1}".format(attribute[0], attribute[1]))

        # groups for each times
        print("\n  # __ list of groups (time) in species.h5 :\n")
        times = list(f.keys())
        for time in times :
            print("  {0}".format(time))

        # subgroups (for species) are the same for each time group
        # let's see the content of the first one
        print("\n  # __ list of sub-groups (specie) in each time group :\n")
        time0 = times[0]
        species = list(f.get(time0).items())
        for specie in species :
            print("  {0}".format(specie[0]))

        # datasets (all the sames in a given group/subgroup
        print("\n  # __ list of associated dataset (for each sub-group) :\n")
        specie0 = species[0][0]
        firstTimeGrpSpecieGrp = '/'+time0+'/'+specie0
        inTheGroup = f.get(firstTimeGrpSpecieGrp)
        for dataSet in inTheGroup.items() :
           print("  {0}".format(dataSet[0]))

        f.close()

    else :
        print("\n### file species.h5 does not exist")


    if os.path.isfile(timeFile) is True :
        # _____ time.h5 _____
        f = h5py.File(timeFile)
        print("\n### file : {0} ###".format(timeFile))

        # datasets
        print("\n  # __ list of associated dataset (for each sub-group) :\n")
        inRoot = f.get('/')
        for dataSet in inRoot.items() :
           print("  {0}".format(dataSet[0]))

        f.close()
    else :
        print(time)
        print("\n### file time.h5 does not exist")


    if os.path.isfile(restartFile) is True :
        # _____ restarts.h5 _____
        f = h5py.File(restartFile)
        print("\n### file : {0} ###".format(restartFile))

        # attribute
        print("\n  # __ list of attributes :\n")
        attributes = f.attrs.items()
        for atribute in attributes :
            print("  attribute [ {0} ] : {1}".format(atribute[0], atribute[1]))

        # groups for each times
        print("\n  # __ list of groups (time) in species.h5 :\n")
        times = list(f.keys())
        for time in times :
            print("  {0}".format(time))

        # attribute
        inTheGroup = f.get(times[0])
        print("\n  # __ list of attributes in each time group :\n")
        attributes = inTheGroup.attrs.items()
        for atribute in attributes :
            print("  attribute [ {0} ] : {1}".format(atribute[0], atribute[1]))

        # subgroups (for species) are the same for each time group
        # let's see the content of the first one
        print("\n  # __ list of sub-groups (ranks) in each time group :\n")
        ranks = list(inTheGroup.items())
        for rank in ranks :
            print("  {0}".format(rank[0]))

        # attribute
        rank0 = ranks[0][0]
        firstTimeGrpRankGrp = '/'+times[0]+'/'+rank0
        inTheGroup = f.get(firstTimeGrpRankGrp)
        attributes = inTheGroup.attrs.items()
        print("\n  # __ list of attributes in each rank sub-group :\n")
        for attribute in attributes :
            #print(attrib[0], attrib[1])
            print("  attribute [ {0} ] : {1}".format(attribute[0], attribute[1]))

        # datasets
        print("\n  # __ list of associated dataset (for each sub-group) :\n")
        for dataSet in inTheGroup.items() :
           print("  {0}".format(dataSet[0]))

        f.close()
    else :
        print("\n### file restarts.h5 does not exist")

else :
    print("\n### the dir {0} does not exist".format(path))
