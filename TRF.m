
 
 for m=0:1:3
   
 r=30;    
 
 ncidv0 = netcdf.open(['ANLEC_r' num2str(r) 'h4v1_TRF_phase.' num2str(r,'%04d') '-v' num2str(m) '.nc'], 'NC_NOWRITE') ; 
 varnamev0=netcdf.inqVar(ncidv0,0);
 varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
 datav0=netcdf.getVar(ncidv0,varidv0);
 datav1=double(datav0);
 
 nx=size(datav0,1);
 ny=size(datav0,2);
 nz=size(datav0,3);
 
 parpool(24)
 parfor z=1:nz
     for y=1:ny
         for x=1:nx

             datav1(x,y,z)=datav0(nx+1-x,ny+1-y,z);
 
         end
     end
 end
  
 delete(gcp('nocreate'))
%  
% 
%  k=2;   i='0110'; 


 
 % recovery of classified results - data of MFs after renormalization using GMM.
filename = ['ANLEC_r' num2str(r) 'h4v1_TRF_phase.' num2str(r,'%04d') '-v' num2str(m) '_Transformed180.nc'];
voxelsize = 1.0;

 
origin = [0 0 0];
valid_range=single([min(min(min(datav1))) max(max(max(datav1)))]);
%valid_range=single([min(min(min(datav1))) max(max(max(datav1)))]);
  %morphy_netCDF(a, voxelsize, filename)

 
ncid = netcdf.create(filename,'CLOBBER');
tomo_xdimID = netcdf.defDim(ncid,'tomo_xdim',nx);
tomo_ydimID = netcdf.defDim(ncid,'tomo_ydim',ny);
tomo_zdimID = netcdf.defDim(ncid,'tomo_zdim',nz);
attGlob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,attGlob,'number_of_files',1);
netcdf.putAtt(ncid,attGlob,'voxel_size_xyz',single([voxelsize voxelsize voxelsize]));
%netcdf.putAtt(ncid,attGlob,'voxel_size_xyz',single([voxelsize voxelsize voxelsize]));
netcdf.putAtt(ncid,attGlob,'voxel_unit','um');
netcdf.putAtt(ncid,attGlob,'zdim_range',int32([0 nz-1]));
netcdf.putAtt(ncid,attGlob,'zdim_total',int32(nz));
netcdf.putAtt(ncid,attGlob,'coordinate_origin_xyz',origin);
netcdf.putAtt(ncid,attGlob,'history_gen','Matlab write');
varid = netcdf.defVar(ncid,'tomo','NC_FLOAT',[tomo_xdimID,tomo_ydimID,tomo_zdimID]);
netcdf.putAtt(ncid,varid,'data_description','drop');
netcdf.putAtt(ncid,varid,'valid_range',valid_range);
netcdf.putAtt(ncid,varid,'_FillValue',single(-127));
%netcdf.putAtt(ncid,varid,'_FillValue',int8(-127));
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,datav1);
netcdf.close(ncid);

clear

 end