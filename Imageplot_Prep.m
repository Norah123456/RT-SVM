 %%GMM_classified

 clear

 r=30; nx=900;ny=400;


 ncidv0 = netcdf.open(['classification_' num2str(r) '_innerc_4_RT2_renormalized_unsupervised_resampled_test.nc']); 
 varnamev0=netcdf.inqVar(ncidv0,0);
 varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
 datav0=netcdf.getVar(ncidv0,varidv0);
 datav0=double(datav0);
 
  %%Image transform 180 degrees.
 
 temp=flip(datav0,1);
 datav1=flip(temp,2);
 
 %Image subset
 datav2=datav1(1:nx,1:ny,:);
nz=size(datav2,3);
 
 % recovery of classified results - data of MFs after renormalization using GMM.
filename = ['segmentedclassification_' num2str(r) '_innerc_4_RT2_renormalized_unsupervised_resampled_test_TRF_nx' num2str(nx) 'ny' num2str(ny) 'nz' num2str(nz) '.nc'];
voxelsize = 1.0;

origin = [0 0 0];
valid_range=single([min(min(min(datav2))) max(max(max(datav2)))]);
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
%netcdf.putAtt(ncid,varid,'_FillValue',int16(-127));
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,datav2);
netcdf.close(ncid);


clear datav* nc* t* origin f* a* v*

%%SVM_classified

ncidv0 = netcdf.open(['classification_' num2str(r) '_RT2_supervised_resampled_SVMrbf_withboundary_allcentroid_Transformed180.nc']) ; 
 varnamev0=netcdf.inqVar(ncidv0,0);
 varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
 datav0=netcdf.getVar(ncidv0,varidv0);
 datav0=double(datav0);
 
  %%Image transform 180 degrees.
 
 temp=flip(datav0,1);
 datav1=flip(temp,2);
 
 %Image subset
 datav2=datav1(1:nx,1:ny,:);
nz=size(datav2,3);
 
 % recovery of classified results - data of MFs after renormalization using GMM.
filename = ['segmentedclassification_' num2str(r) '_RT2_supervised_resampled_SVMrbf_withboundary_allcentroid_Transformed180_TRF_nx' num2str(nx) 'ny' num2str(ny) 'nz' num2str(nz) '.nc'];
voxelsize = 1.0;

origin = [0 0 0];
valid_range=single([min(min(min(datav2))) max(max(max(datav2)))]);
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
netcdf.putVar(ncid,varid,datav2);
netcdf.close(ncid);

