
import h5py


filePath = '/home/smets/shErpA/blAckDog/run/b3/'
initialFileName = 'fields.h5'
extractFileName = 'extract.h5'

fi = h5py.File(filePath+initialFileName)
fe = h5py.File(filePath+extractFileName)

times = list(fi.keys())
#print("groups : ", times)

attributes = list(fi.attrs.keys())

for attribute in attributes:
    fe.attrs[attribute] = fi.attrs[attribute]

times  = [24.0, 36.0]

for time in times:
    initialGrp =  'time : {0:.6f}'.format(time)
    print("extracted groups : ", initialGrp)
    grpPath = fi[initialGrp].parent.name
    extractGrp = fe.require_group(grpPath)
    fi.copy(initialGrp, extractGrp, name=initialGrp)

fi.close()
fe.close()

