#% read_bkm : read braink .bkm file
#%
#% This file typically contains the cortical surface mesh
#%
#% Inputs :
#% fname - filename of bkm file (must be on matlab search path)
#%
#% Outputs :
#% R - structure containing fields
#%  R.endinanness - 'l' or 'b' indicating endianness of file
#%  R.npoints - number of points in mesh
#%  R.x, R.y, R.z - npoints size arrays with x,y,z coordinates of mesh points
#%  R.ntriangles - number of triangles in mesh
#%  R.tri - ntriangles x 3 array, kth row gives indices of kth triangle
#%   e.g, R.x(R.tri(k,1)) is x coordinate of first vertex of kth triangle
#%  This format is compatible with matlab's trisurf
#%
#%  R.nindices - number of indices (if voxels are given) : should be same as
#%  npoints, or 0
#%  R.voxindices - nindices size array giving voxel of each index
#%
#% Author: David Hammond, NeuroInformatics Center, University of Oregon
#% Date : May, 2010
#
#function R=read_bkm(fname)
#  fid=fopen(fname,'r');
#  if fid==-1
#    error('unable to open file %s',fname);
#  end
#  fullname=which(fname);
#  if isempty(fullname)
#    % happens if fname not on path ... in which case it is a fullname
#    fullname=fname;
#  end
#  nbytes=getfield(dir(fullname),'bytes');
#
#  % endianness
#  tmp=fread(fid,1,'uint32',0,'l');
#  if tmp==0
#    R.endianness='l';
#  else
#    R.endianness='b';
#  end
#  mf=R.endianness;
#
#  % npoints
#  R.npoints=fread(fid,1,'uint32',0,mf);
#
#  % point coordinates
#  % now read x(1),y(1),z(1),x(2),y(2),z(2)...
#  all_pc = fread(fid,3*R.npoints,'float32',0,mf);
#  R.x=all_pc(1:3:end);
#  R.y=all_pc(2:3:end);
#  R.z=all_pc(3:3:end);
#
#  % ntriangles
#  R.ntriangles=fread(fid,1,'uint32',0,mf);
#
#  % triangle point id's
#  all_tpi=fread(fid,3*R.ntriangles,'uint32');
#  % careful to pack triangles in in correct order to give ntriangles x 3 matrix
#  % add 1 to shift to 1-based indexing for matlab
#  R.tri=reshape(all_tpi,[3 R.ntriangles])' + 1;
#  % number of indices
#  R.nindices=fread(fid,1,'uint32',0,mf);
#
#  if R.nindices>0
#    R.voxindices=fread(fid,R.nindices,'uint32',0,mf);
#  end
#
#  % we should be at end of file. Read one more character to check
#  fread(fid,1,'char');
#  if ~feof(fid)
#    warning('Expected end of file ... perhaps file was not interpreted correctly?');
#  end
#
#  fclose(fid);
