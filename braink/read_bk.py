# -*- coding: utf-8 -*-

import numpy as np

##--------------------------------------------------------------------------------------
def lfm(lfmfilename, nE):
    nC = nE + 1 
    fd = open(lfmfilename, 'rb')
    Kori = np.fromfile(file=fd, dtype=np.dtype('>d'))
    fd.close()
 
    nV =  Kori.shape[0] / nE
    Kori = Kori.reshape(nE, nV)
 
    H = np.eye((nC))- np.ones((nC, nC))/nC
    Kori2 = np.concatenate((Kori, np.zeros((1, nV)))) 
    K = np.dot(H, Kori2) 
    return K 

##--------------------------------------------------------------------------------------
def dsa(dsafilename):
    fd = open(dsafilename, 'rb')
    A_dsa = np.fromfile(file=fd, dtype=np.dtype('i4'))
    fd.close()

    nV = int(np.sqrt(A_dsa.shape[0]))
    A_dsa = A_dsa.reshape(nV, nV)

    return A_dsa 


##---------------------------------------------------------------------------------------
def bkd(fname):
    #fname = '/Users/jsong/Python2.7/108/108_HeadModel_MRI/108_01_2000.Dipoles_1_MRI.bkd' 

    fid = open(fname,'r') 
    #fullname = fid.name  #nbytes = os.path.getsize(fullname)
    n_dipoles = np.fromfile(file=fid, count=1, dtype=np.dtype('<u4'), sep="")[0]

    #% load dipole location indices 
    all_dli = np.fromfile(file=fid, count= 3 * n_dipoles,dtype=np.dtype('<i4'),sep="")
    dipole_location_index_x=all_dli[0::3]
    dipole_location_index_y=all_dli[1::3] 
    dipole_location_index_z=all_dli[2::3]
  
    # load dipole location float values
    n_triples_with_fv = np.fromfile(file=fid, count=1, dtype=np.dtype('<u4'),sep="")[0]
    if not(n_triples_with_fv==0 or n_triples_with_fv == n_dipoles):
        print('Invalid n_triples_with_fv \n') 
    elif (n_triples_with_fv>0):
        all_dlf = np.fromfile(file=fid, count=3*n_triples_with_fv, dtype=np.dtype('<f4'),sep="")
        dipole_location_x=all_dlf[0::3] 
        dipole_location_y=all_dlf[1::3] 
        dipole_location_z=all_dlf[2::3]
    else:
        dipole_location_x=[] 
        dipole_location_y=[] 
        dipole_location_z=[] 
  
    # load oriented normals
    n_oriented_normals = np.fromfile(file=fid, count=1, dtype=np.dtype('<u4'),sep="")[0] #% number of oriented normals
    if  not(n_oriented_normals==0 or n_oriented_normals == n_dipoles):
        print('invalid n_oriented_normals\n') 
    elif (n_oriented_normals>0):
        all_dipole_orientations = np.fromfile(file=fid, count=3*n_oriented_normals, dtype=np.dtype('<f4'),sep="") 
        dipole_orientation_x=all_dipole_orientations[0::3] 
        dipole_orientation_y=all_dipole_orientations[1::3]  
        dipole_orientation_z=all_dipole_orientations[2::3] 
    else:
        dipole_orientation_x=[] 
        dipole_orientation_y=[] 
        dipole_orientation_z=[] 
    
    #% load dipole distances
    n_dipole_dc = np.fromfile(file=fid, count=1, dtype=np.dtype('<u4'),sep="")[0] #number of dipoles with calculated distances
    if not(n_dipole_dc == 0 or n_dipole_dc == n_dipoles):
        print('invalid n_dipole_dc\n')
    elif (n_dipole_dc > 0):
        all_dipole_distances = np.fromfile(file=fid, count=n_dipole_dc**2, dtype=np.dtype('<f4'),sep="") 
        dipole_distances = all_dipole_distances.reshape(n_dipole_dc, n_dipole_dc) 
    else:
        dipole_distances = [] 

    fid.close()

    R = {"ndipoles": n_dipoles, 
    "dipole_location_index":{"x":dipole_location_index_x, "y":dipole_location_index_y, "z":dipole_location_index_z}, 
    "dipole_location":{"x":dipole_location_x, "y":dipole_location_y, "z":dipole_location_z},
    "dipole_orientation":{"x":dipole_orientation_x, "y":dipole_orientation_y, "z":dipole_orientation_z},
    "dipole_distances":dipole_distances} 

    return R

##------------------------------------------------------------------------------

def bkm(fname):
    
# fname = '/Users/jsong/Python2.7/108/108_HeadModel_MRI/108_01_2000.Cortical_surface_1_MRI.bkm' 
# This file typically contains the cortical surface mesh
#
# Inputs :
# fname - filename of bkm file (must be on matlab search path)
#
# Outputs :
# R - structure containing fields
#  R.endinanness - 'l' or 'b' indicating endianness of file
#  R.npoints - number of points in mesh
#  R.x, R.y, R.z - npoints size arrays with x,y,z coordinates of mesh points
#  R.ntriangles - number of triangles in mesh
#  R.tri - ntriangles x 3 array, kth row gives indices of kth triangle
#   e.g, R.x(R.tri(k,1)) is x coordinate of first vertex of kth triangle
#  This format is compatible with matlab's trisurf
#
#  R.nindices - number of indices (if voxels are given) : should be same as
#  npoints, or 0
#  R.voxindices - nindices size array giving voxel of each index
#
# Author: David Hammond, NeuroInformatics Center, University of Oregon
# Date : May, 2010
    
    fid = open(fname,'r') 
    # endianness
    tmp = np.fromfile(file=fid, count=1, dtype=np.dtype('<u4'),sep="")[0] 
    if tmp == 0:
        endianness='l' 
    else:
        endianness='b'  

    # npoints
    npoints = np.fromfile(file=fid, count=1, dtype=np.dtype('<u4'),sep="")[0]  
    
    # point coordinates
    # now read x(1),y(1),z(1),x(2),y(2),z(2)...
    all_pc = np.fromfile(file=fid, count=3*npoints, dtype=np.dtype('<f4'),sep="")
    x=all_pc[0::3]
    y=all_pc[1::3]
    z=all_pc[2::3]

    # ntriangles
    ntriangles = np.fromfile(file=fid, count=1, dtype=np.dtype('uint32'),sep="")[0]   

    # triangle point id's
    all_tpi = np.fromfile(file=fid, count=3*ntriangles, dtype=np.dtype('uint32'))# ,sep="") 
    # careful to pack triangles in correct order to give ntriangles x 3 matrix
    # add 1 to shift to 1-based indexing for matlab
    all_tpi =all_tpi.reshape(ntriangles,3) 
    tri = all_tpi + 1

    #% number of indices
    nindices = np.fromfile(file=fid, count=1, dtype=np.dtype('uint32'),sep="")[0]  

    if nindices > 0:
        voxindices = np.fromfile(file=fid, count=nindices, dtype=np.dtype('uint32'), sep="")  
 
    fid.close()

    R = {"endianness":endianness, "npoints":npoints, "x":x, "y":y, "z":z, 
        "ntriangles":ntriangles, "tri":tri, "nindices":nindices,"voxindices":voxindices} 

    return R
    
##------------------------------------------------------

#fname = '/Users/jsong/Python2.7/108/108_HeadModel_MRI/108_01_2000.Chaco_assignment_1_MRI.bka' 
def bka(fname):

    matrix = []
    for i in open(fname, 'r'):
        matrix.append( map(int, i.split()) )

    myarray = np.asarray(matrix) 
    numtriangles = myarray[0][0]
    dipolenumbers = myarray[1:, 0]
    dipolenumbers = dipolenumbers+1 #; % convert to 1 based indexing

    R = {"numtriangles":numtriangles, "dipolenumbers":dipolenumbers}
        
    return R

##------------------------------------------------------
#fname = '/Users/jsong/Python2.7/108/108_HeadModel_MRI/108_01_2000.Sensors_GPS_MRI.bks' 
def bks(fname):
 
    elname = []
    elpos_x = []
    elpos_y = []
    elpos_z = []

    alist = [line.rstrip() for line in open(fname)]

    for i in range(len(alist)):
        list_of_items_in_line = alist[i].split(" ")
        if list_of_items_in_line[0][0] == 'E': 
            elname.append(list_of_items_in_line[0])
            elpos_x.append(list_of_items_in_line[1])
            elpos_y.append(list_of_items_in_line[2])
            elpos_z.append(list_of_items_in_line[3])

        elif list_of_items_in_line[0] == 'FidNz':
            FidNz = [list_of_items_in_line[1], list_of_items_in_line[2],list_of_items_in_line[3]]
            FidNz = map(float, FidNz) 

        elif list_of_items_in_line[0] == 'FidT9':
            FidT9 = [list_of_items_in_line[1], list_of_items_in_line[2],list_of_items_in_line[3]]
            FidT9 = map(float, FidT9) 

        elif list_of_items_in_line[0] == 'FidT10':
            FidT10 = [list_of_items_in_line[1], list_of_items_in_line[2],list_of_items_in_line[3]]
            FidT10 = map(float, FidT10) 

        elif list_of_items_in_line[0] == 'Cz':
            Cz = [list_of_items_in_line[1], list_of_items_in_line[2],list_of_items_in_line[3]]
            Cz = map(float, Cz) 
 
    elpos_x = map(float, elpos_x) 
    elpos_y = map(float, elpos_y) 
    elpos_z = map(float, elpos_z) 
    xx = np.matrix(elpos_x) 
    yy = np.matrix(elpos_y)
    zz = np.matrix(elpos_z)
    elpos = np.concatenate((xx, yy,zz ), axis=0)

    R = {'elpos':elpos, 'elname':elname, 'FidNz':FidNz, 'FidT9':FidT9, 'FidT10':FidT10, 'Cz':Cz}
    
    return R
    
##------------------------------------------------------
