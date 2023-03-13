% format shortg
% c1=clock

clear


diary('DiaryGMM_unsupervised_Transformed180.txt');
diary on;

r=30;
% rock types to classify
k=2;
%data transformation from nc file.
% 
 ncidv0 = netcdf.open(['ANLEC_r' num2str(r) 'h4v1_TRF_phase.' num2str(r, '%04d') '-v0_Transformed180.nc'], 'NC_NOWRITE') ; 
 varnamev0=netcdf.inqVar(ncidv0,0);
 varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
 datav0_raw=netcdf.getVar(ncidv0,varidv0);
 datav0_raw=double(datav0_raw);
 
 %data transformation from nc file.

 ncidv1 = netcdf.open(['ANLEC_r' num2str(r) 'h4v1_TRF_phase.' num2str(r, '%04d') '-v1_Transformed180.nc'], 'NC_NOWRITE') ; 
 varnamev1=netcdf.inqVar(ncidv1,0);
 varidv1 = netcdf.inqVarID(ncidv1,varnamev1);
 datav1_raw=netcdf.getVar(ncidv1,varidv1);
 datav1_raw=double(datav1_raw);
 %data transformation from nc file.

 ncidv2 = netcdf.open(['ANLEC_r' num2str(r) 'h4v1_TRF_phase.' num2str(r, '%04d') '-v2_Transformed180.nc'], 'NC_NOWRITE') ; 
 varnamev2=netcdf.inqVar(ncidv2,0);
 varidv2 = netcdf.inqVarID(ncidv2,varnamev2);
 datav2_raw=netcdf.getVar(ncidv2,varidv2);
 datav2_raw=double(datav2_raw); %data transformation from nc file.

 ncidv3 = netcdf.open(['ANLEC_r' num2str(r) 'h4v1_TRF_phase.' num2str(r, '%04d') '-v3_Transformed180.nc'], 'NC_NOWRITE');  
varnamev3=netcdf.inqVar(ncidv3,0);
 varidv3 = netcdf.inqVarID(ncidv3,varnamev3);
 datav3_raw=netcdf.getVar(ncidv3,varidv3);
 datav3_raw=double(datav3_raw);

  % adjust the sidelength of z dimension to the same as r110.
 dist=110-r;
 nz_raw=length(datav0_raw(1,1,:));
 datav0=datav0_raw(:,:,(dist+1):(nz_raw-dist));
 datav1=datav1_raw(:,:,(dist+1):(nz_raw-dist));
 datav2=datav2_raw(:,:,(dist+1):(nz_raw-dist));
 datav3=datav3_raw(:,:,(dist+1):(nz_raw-dist));
 
 
nx=size(datav0,1);
ny=size(datav0,2);
nz=size(datav0,3);


datav0_innerc=datav0(2:2:end,2:2:end,2:2:end);
datav1_innerc=datav1(2:2:end,2:2:end,2:2:end);
datav2_innerc=datav2(2:2:end,2:2:end,2:2:end);
datav3_innerc=datav3(2:2:end,2:2:end,2:2:end);
% datav4_innerc=datav4_temp(2:2:end,2:2:end,2:2:end);



nx=length(datav0_innerc(:,1,1));
ny=length(datav0_innerc(1,:,1));
nz=length(datav0_innerc(1,1,:));


 

 % FIRST CLASSIFICATION
 % data preparation for Gaussian mixture model calculation with morphological, electrical and elastic properties.
 X_DEM(:,1)=datav0_innerc(:);
 X_DEM(:,2)=datav1_innerc(:);
 X_DEM(:,3)=datav2_innerc(:);
 X_DEM(:,4)=datav3_innerc(:);
%  X_DEM(:,5)=DCInterpdata_DEM(:);
%  X_DEM(:,6)=DCInterpdata_ELAairDEM_K(:);
%  X_DEM(:,7)=DCInterpdata_ELAairDEM_G(:);
%  X_DEM(:,5)=sa_innerc(:);
 
 % data preparation for Gaussian mixture model calculation with morphological properties.
 
 X_DEM_1(:,1:4)=X_DEM(:,1:4);
 
  
 %data renormalization
 
 mu=mean(datav0_innerc(:));
 sig=std(datav0_innerc(:));
 a=1;

  parfor m=1:nx
     for n=1:ny
         for j=1:nz
             datav0_innerc_1(m,n,j)=(datav0_innerc(m,n,j)-mu)/sig;
             a=a+1;
         end
     end
 end
 delete(gcp('nocreate'))
 
 %renormalize the data of surface volume.
 
 mu=mean(datav1_innerc(:));
 sig=std(datav1_innerc(:));
 a=1;
  parfor m=1:nx
     for n=1:ny
         for j=1:nz
             datav1_innerc_1(m,n,j)=(datav1(m,n,j)-mu)/sig;
             a=a+1;
         end
     end
 end
delete(gcp('nocreate'))

 %renormalize the data of mean curvature.
 
 mu=mean(datav2_innerc(:));
 sig=std(datav2_innerc(:));
 a=1;
  parfor m=1:nx
     for n=1:ny
         for j=1:nz
             datav2_innerc_1(m,n,j)=(datav2_innerc(m,n,j)-mu)/sig;
             a=a+1;
         end
     end
 end
 delete(gcp('nocreate'))

 %renormalize the data of total curvature.
 
 mu=mean(datav3_innerc(:));
 sig=std(datav3_innerc(:));
 a=1;
  parfor m=1:nx
     for n=1:ny
         for j=1:nz
            datav3_innerc_1(m,n,j)=(datav3_innerc(m,n,j)-mu)/sig;
             a=a+1;
         end
     end
 end
delete(gcp('nocreate')) 

        
 
  % data preparation for Gaussian mixture model calculation with morphological, electrical and elastic properties with renormalization.
 
 X_DEM_2(:,1)=datav0_innerc_1(:);
 X_DEM_2(:,2)=datav1_innerc_1(:);
 X_DEM_2(:,3)=datav2_innerc_1(:);
 X_DEM_2(:,4)=datav3_innerc_1(:);
%  X_DEM_2(:,5)=DCInterpdata_DEM_1(:);
%  X_DEM_2(:,6)=DCInterpdata_ELAairDEM_K_1(:);
%  X_DEM_2(:,7)=DCInterpdata_ELAairDEM_G_1(:);
%  X_DEM_2(:,5)=sa_innerc_1(:);  
 % data preparation for Gaussian mixture model calculation with morphological properties with renormalization.
 
 X_DEM_3(:,1:4)=X_DEM_2(:,1:4);


 
clear datav*
format shortg
c1=clock

% classification using MFs and Sw.

parpool('local',24,'IdleTimeout',120)  % parallel computing


% % Classification of 4 attributes using data without renormalization.
% options_1=statset('MaxIter', 40000);
% obj_1=fitgmdist(X_DEM_1,k,'regularize', 0.0001, 'options', options_1);
% idx_1=cluster(obj_1,X_DEM_1);
% cluster1_1=X_DEM_1(idx_1==1,:);
% cluster2_1=X_DEM_1(idx_1==2,:);
% cluster3_1=X_DEM_1(idx_1==3,:);
% Classification of 4 attributes using data with renormalization.

options_3=statset('MaxIter', 40000);
obj_3=fitgmdist(X_DEM_3,k,'regularize', 0.0001, 'options', options_3);
idx_3=cluster(obj_3,X_DEM_3);
cluster1_3=X_DEM_3(idx_3==1,:);
cluster2_3=X_DEM_3(idx_3==2,:); 
cluster3_3=X_DEM_3(idx_3==3,:);
% % Classification of 5 attributes(MFs and Sw) using data without renormalization.
% 
% 
% options_8=statset('MaxIter', 40000);
% obj_8=fitgmdist(X_DEM_8,k,'regularize', 0.0001, 'options', options_8);
% idx_8=cluster(obj_8,X_DEM_8);
% cluster1_8=X_DEM_8(idx_8==1,:);
% cluster2_8=X_DEM_8(idx_8==2,:);
% 
% % Classification of 5 attributes(MFs and Sw) using data with renormalization.
% 
% options_9=statset('MaxIter', 40000);
% obj_9=fitgmdist(X_DEM_9,k,'regularize', 0.0001, 'options', options_9);
% idx_9=cluster(obj_9,X_DEM_9);
% cluster1_9=X_DEM_9(idx_9==1,:);
% cluster2_9=X_DEM_9(idx_9==2,:);


delete(gcp('nocreate')) % shut down parallel computing.

format shortg
c2=clock
% clear X_DEM*  ncidv* varnamev* varidv* datav*
% c2=clock
% save(['classification_' num2str(i) '_RT' num2str(k) '_innerc.mat'],'-v7.3')

rock_type(1,1)=length(cluster1_3(:,1));
rock_type(1,2)=length(cluster2_3(:,1));
rock_type(1,3)=rock_type(1,1)/(rock_type(1,1)+rock_type(1,2));
rock_type(1,4)=rock_type(1,2)/(rock_type(1,1)+rock_type(1,2));


data=rock_type;
[m,n]=size(data);
data_cell=mat2cell(data,ones(m,1),ones(n,1));
% title={'MF1','MF2','MF3','MF1_re','MF2_re','MF3_re'};
title={'MF1_re','MF2_re','Prop_RT1','Prop_RT2'};
result=[title;data_cell];
writecell(result,['rock_type_' num2str(r) '_RT' num2str(k) '_unsupervised_Transformed180.xls']);
 




% recovery of classified results - data of MFs after renormalization using GMM.
filename = ['classification_' num2str(r) '_innerc_4_RT' num2str(k) '_renormalized_unsupervised_resampled_test.nc'];
voxelsize = 1.0;

a=zeros(nx, ny, nz);

b=1;
 for m=1:nz
     for n=1:ny
         for j=1:nx
             a(j,n,m)=idx_3(b);
             b=b+1;
         end
     end
 end

% recover to original size (2*2*2) change a point into a cube of the same
% values.
a_resample=zeros(2*nx,2*ny,2*nz);
for z=1:nz
    for y=1:ny
        for x=1:nx
            x1=(x-1)*2+1;
            y1=(y-1)*2+1;
            z1=(z-1)*2+1;
            x2=2*x;
            y2=2*y;
            z2=2*z;
          a_resample(x1:x2,y1:y2,z1:z2)=a(x,y,z);  
            
        end
    end
end 

 nx=2*nx;
 ny=2*ny;
 nz=2*nz;
 
 
origin = int32([0 0 0]);
valid_range=int32([min(min(min(a))) max(max(max(a)))]);
  %morphy_netCDF(a, voxelsize, filename)

 
ncid = netcdf.create(filename,'CLOBBER');
seg_xdimID = netcdf.defDim(ncid,'segmented_xdim',nx);
seg_ydimID = netcdf.defDim(ncid,'segmented_ydim',ny);
seg_zdimID = netcdf.defDim(ncid,'segmented_zdim',nz);
attGlob = netcdf.getConstant('NC_GLOBAL');
netcdf.putAtt(ncid,attGlob,'number_of_files',1);
netcdf.putAtt(ncid,attGlob,'voxel_size_xyz',single([voxelsize voxelsize voxelsize]));
netcdf.putAtt(ncid,attGlob,'voxel_unit','um');
netcdf.putAtt(ncid,attGlob,'zdim_range',int32([0 nz-1]));
netcdf.putAtt(ncid,attGlob,'zdim_total',int32(nz));
netcdf.putAtt(ncid,attGlob,'coordinate_origin_xyz',origin);
netcdf.putAtt(ncid,attGlob,'history_gen','Matlab write');
varid = netcdf.defVar(ncid,'segmented','byte',[seg_xdimID,seg_ydimID,seg_zdimID]);
netcdf.putAtt(ncid,varid,'data_description','drop');
netcdf.putAtt(ncid,varid,'valid_range',valid_range);
netcdf.putAtt(ncid,varid,'_FillValue',int8(-127));
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,a_resample);
netcdf.close(ncid);


nx=size(a,1);
ny=size(a,2);
nz=size(a,3);

% verification of classification results

ncidv0 = netcdf.open('classification_0110_innerc_4_RT2_renormalized_unsupervised_resampled_Reference2020.nc'); 
varnamev0=netcdf.inqVar(ncidv0,0);% beware of the integer here, 1 or 0, for the calling of intensities;
varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
datav0=netcdf.getVar(ncidv0,varidv0);
datav0=double(datav0);



 datav1=a_resample;
 datav1_comparison=datav1==datav0;
 datav1_comparison=double(datav1_comparison);
 datav1_comparison_true=find(datav1_comparison==1);
 MFs_supervised=length(datav1_comparison_true(:))/length(datav1_comparison(:));
 MFs_supervised

 diary off;


 

 
 
