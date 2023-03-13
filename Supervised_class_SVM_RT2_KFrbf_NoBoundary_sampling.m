
clear

% format shortg
% c1=clock
%data transformation from nc file.
diary('DiarySVMSamplingTest1_1000_NoBoundary.txt');
diary on

 i='30'; r=30;

 
ncidv0 = netcdf.open(['ANLEC_r' num2str(i) 'h4v1_TRF_phase.' num2str(r,'%04d') '-v0_Transformed180.nc'], 'NC_NOWRITE') ;
 varnamev0=netcdf.inqVar(ncidv0,0);
 varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
 datav0_raw=netcdf.getVar(ncidv0,varidv0);
 datav0_raw=double(datav0_raw);
 %data transformation from nc file.

 ncidv1 = netcdf.open(['ANLEC_r' num2str(i) 'h4v1_TRF_phase.' num2str(r,'%04d') '-v1_Transformed180.nc'], 'NC_NOWRITE') ; 
 varnamev1=netcdf.inqVar(ncidv1,0);
 varidv1 = netcdf.inqVarID(ncidv1,varnamev1);
 datav1_raw=netcdf.getVar(ncidv1,varidv1);
 datav1_raw=double(datav1_raw);
 %data transformation from nc file.

 ncidv2 = netcdf.open(['ANLEC_r' num2str(i) 'h4v1_TRF_phase.' num2str(r,'%04d') '-v2_Transformed180.nc'], 'NC_NOWRITE') ; 
 varnamev2=netcdf.inqVar(ncidv2,0);
 varidv2 = netcdf.inqVarID(ncidv2,varnamev2);
 datav2_raw=netcdf.getVar(ncidv2,varidv2);
 datav2_raw=double(datav2_raw); %data transformation from nc file.

 ncidv3 = netcdf.open(['ANLEC_r' num2str(i) 'h4v1_TRF_phase.' num2str(r,'%04d') '-v3_Transformed180.nc'], 'NC_NOWRITE') ; 
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
 
 clear datav*_raw
 
  
 % pick every other data point for a smaller dataset (faster computation with little information loss).
 datav0_innerc=datav0(2:2:end,2:2:end,2:2:end);
 datav1_innerc=datav1(2:2:end,2:2:end,2:2:end);
 datav2_innerc=datav2(2:2:end,2:2:end,2:2:end);
 datav3_innerc=datav3(2:2:end,2:2:end,2:2:end);
 
 clear datav0 datav1 datav2 datav3

 
 % the size of the shrinked dataset. 
 nx=size(datav0_innerc,1);
 ny=size(datav0_innerc,2);
 nz=size(datav0_innerc,3);

 % FIRST CLASSIFICATION
 % data preparation for SVM with morphological properties(MFs).
 X_DEM_1(:,1)=datav0_innerc(:);
 X_DEM_1(:,2)=datav1_innerc(:);
 X_DEM_1(:,3)=datav2_innerc(:);
 X_DEM_1(:,4)=datav3_innerc(:);

 
%Import supervised information from original Gaussian Random Field.

ncidv1 = netcdf.open('classification_0110_innerc_4_RT2_renormalized_unsupervised_resampled_Reference2020.nc'); 
varnamev1=netcdf.inqVar(ncidv1,0);% beware of the integer here, 1 or 0, for the calling of intensities;
varidv1 = netcdf.inqVarID(ncidv1,varnamev1);
data_class1=netcdf.getVar(ncidv1,varidv1);
data_class1=double(data_class1);

% keep up with the original dataset.
data_class1_innerc=data_class1(2:2:end,2:2:end,2:2:end);

asp=1;  rock_types=2;

sidelen=[nx,ny,nz];
data_class1_innerc_fft=fftn(data_class1_innerc);

% find out the biggest window available for the SVM.
for sl=200:-10:10
    sl
    sl_window=[sl,sl,sl]; 
    side_start=[sl,sl,sl]; 

    data_conv=zeros(nx,ny,nz);
    data_conv(1:sl,1:sl,1:sl)=1; 
    data_conv=circshift(data_conv,-sl/2,1);
    data_conv=circshift(data_conv,-sl/2,2);
    data_conv=circshift(data_conv,-sl/2,3);
    data_conv_fftx=fftn(data_conv);
    
    C=data_class1_innerc_fft.*data_conv_fftx;
    data_window=real(ifftn(C));
    % Get the average value of the current window
    data_window=data_window/(sl_window(1,1)*sl_window(1,2)*sl_window(1,3));
    
    %Get rid of the boundary caused by shifting the window.
    %considering data_window(50,50,50)=mean(data_class1_innerc(1:100,1:100,1:100),'all')
    data_window1=data_window((sl_window(1)/2):(nx-sl_window(1)/2-1),(sl_window(2)/2):(ny-sl_window(2)/2-1),(sl_window(3)/2):(nz-sl_window(3)/2)-1);
    
    if (sum(data_window1(:)==1.00000)>0) && (sum(data_window1(:)==2.00000)>0)
        fprintf('the sidelength of the cubic window is: %6d \n',sl_window);
        break       
    end

end    

clear data_conv_fftx data_class1_innerc_fft

[x1,y1,z1] = ind2sub(size(data_window1),find(data_window1 == 1.00000));
[x2,y2,z2] = ind2sub(size(data_window1),find(data_window1 == 2.00000));


% The position of the centroids of each RT (we use the first spot respectively)
wi1=round(length(x1)/2); % we are going to deal with the first of the proper window spots.
wi2=round(length(x2)/2);
Cposition1=[x1(wi1),y1(wi1),z1(wi1)];
Cposition2=[x2(wi2),y2(wi2),z2(wi2)];
% check if the window is correctly determined.
fprintf('mean value of the chosen window of RT1: %14.6e \n',mean(data_class1_innerc(x1(wi1):(x1(wi1)+sl-1),y1(wi1):(y1(wi1)+sl-1),z1(wi1):(z1(wi1)+sl-1)),"all"));
fprintf('mean value of the chosen window of RT2: %14.6e \n',mean(data_class1_innerc(x2(wi2):(x2(wi2)+sl-1),y2(wi2):(y2(wi2)+sl-1),z2(wi2):(z2(wi2)+sl-1)),"all"));
%-------------------------------------------------------------------------%
% % case2: using XZ layers for supervised classification (with the boundary information between the RTs) %
% 
% ny_start=ny*3/5;
% sample1(:,1)=reshape(datav0_innerc(:,ny_start:5:(ny_start+5*10-1),:),[],1);
% sample1(:,2)=reshape(datav1_innerc(:,ny_start:5:(ny_start+5*10-1),:),[],1);
% sample1(:,3)=reshape(datav2_innerc(:,ny_start:5:(ny_start+5*10-1),:),[],1);
% sample1(:,4)=reshape(datav3_innerc(:,ny_start:5:(ny_start+5*10-1),:),[],1);
% 
% label_SVMclassifier1=reshape(data_class1_innerc(:,ny_start:5:(ny_start+5*10-1),:),[],1);


sample1_1(:,1)=reshape(datav0_innerc(x1(wi1):(x1(wi1)+sl-1),y1(wi1):(y1(wi1)+sl-1),z1(wi1):(z1(wi1)+sl-1)),[],1);
sample1_1(:,2)=reshape(datav1_innerc(x1(wi1):(x1(wi1)+sl-1),y1(wi1):(y1(wi1)+sl-1),z1(wi1):(z1(wi1)+sl-1)),[],1);
sample1_1(:,3)=reshape(datav2_innerc(x1(wi1):(x1(wi1)+sl-1),y1(wi1):(y1(wi1)+sl-1),z1(wi1):(z1(wi1)+sl-1)),[],1);
sample1_1(:,4)=reshape(datav3_innerc(x1(wi1):(x1(wi1)+sl-1),y1(wi1):(y1(wi1)+sl-1),z1(wi1):(z1(wi1)+sl-1)),[],1);

label_SVMclassifier1_1=reshape(data_class1_innerc(x1(wi1):(x1(wi1)+sl-1),y1(wi1):(y1(wi1)+sl-1),z1(wi1):(z1(wi1)+sl-1)),[],1);

sample1_2(:,1)=reshape(datav0_innerc(x2(wi2):(x2(wi2)+sl-1),y2(wi2):(y2(wi2)+sl-1),z2(wi2):(z2(wi2)+sl-1)),[],1);
sample1_2(:,2)=reshape(datav1_innerc(x2(wi2):(x2(wi2)+sl-1),y2(wi2):(y2(wi2)+sl-1),z2(wi2):(z2(wi2)+sl-1)),[],1);
sample1_2(:,3)=reshape(datav2_innerc(x2(wi2):(x2(wi2)+sl-1),y2(wi2):(y2(wi2)+sl-1),z2(wi2):(z2(wi2)+sl-1)),[],1);
sample1_2(:,4)=reshape(datav3_innerc(x2(wi2):(x2(wi2)+sl-1),y2(wi2):(y2(wi2)+sl-1),z2(wi2):(z2(wi2)+sl-1)),[],1);

label_SVMclassifier1_2=reshape(data_class1_innerc(x2(wi2):(x2(wi2)+sl-1),y2(wi2):(y2(wi2)+sl-1),z2(wi2):(z2(wi2)+sl-1)),[],1);

sample1=cat(1,sample1_1,sample1_2);
label_SVMclassifier1=cat(1,label_SVMclassifier1_1,label_SVMclassifier1_2);



clear datav*_innerc sample1_1 sample1_2 label_SVMclassifier1_1 label_SVMclassifier1_2

%--------Temporary test for the time usage of fitcsvm here.
sample1_temp=sample1(1:1000:end,:);
label_SVMclassifier1_temp=label_SVMclassifier1(1:1000:end,:);

format shortg
c1=clock

SVMmodel_temp=fitcsvm(sample1_temp,label_SVMclassifier1_temp,'KernelFunction','rbf','Standardize',true);

format shortg
c2=clock

%--------Temporary test for the time usage of fitcsvm here.



X_DEM_1_temp=X_DEM_1(2:2:end,:);

%parpool('local',24,'IdleTimeout',480) 

[label_X_DEM1,score_X_DEM1] = predict(SVMmodel_temp,X_DEM_1_temp);

%delete(gcp('nocreate'))

format shortg
c3=clock





% recover the size of training data
label_X_DEM1_double=kron(label_X_DEM1,[1;1]);

a=zeros(nx, ny, nz);
Eff=find(a==0);
a(Eff)=label_X_DEM1_double(:);


cluster1_3=X_DEM_1(label_X_DEM1_double==1,:);
cluster2_3=X_DEM_1(label_X_DEM1_double==2,:); 

rock_type(1,1)=length(cluster1_3(:,1));
rock_type(1,2)=length(cluster2_3(:,1));
rock_type(1,3)=rock_type(1,1)/(rock_type(1,1)+rock_type(1,2));
rock_type(1,4)=rock_type(1,2)/(rock_type(1,1)+rock_type(1,2));

k=2;
data=rock_type;
[m,n]=size(data);
data_cell=mat2cell(data,ones(m,1),ones(n,1));
% title={'MF1','MF2','MF3','MF1_re','MF2_re','MF3_re'};
title={'MF1_re','MF2_re','Prop_RT1','Prop_RT2'};
result=[title;data_cell];
writecell(result,['rock_type_' num2str(r) '_RT' num2str(k) '_SVMrbf_Transformed180_NoBoundary.xls']);


 
% recovery of classified results - data of MFs using SVM.
filename = ['classification_' num2str(i) '_RT2_supervised_resampled_SVMrbf_withboundary_allcentroid_Transformed180_NoBoundary.nc'];
voxelsize = 1.0;
 % recover to original size (2*2*2) change a point into a cube of the same
% values.
a_resample=zeros(2*nx,2*ny,2*nz);
for z=1:nz
    for y=1:ny
        for x=1:nx
            rx1=(x-1)*2+1;
            ry1=(y-1)*2+1;
            rz1=(z-1)*2+1;
            rx2=2*x;
            ry2=2*y;
            rz2=2*z;
          a_resample(rx1:rx2,ry1:ry2,rz1:rz2)=a(x,y,z);  
            
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
ncidv0 = netcdf.open('classification_0110_innerc_4_RT2_renormalized_unsupervised_resampled_Reference2020.nc') ; 
varnamev0=netcdf.inqVar(ncidv0,0);
varidv0 = netcdf.inqVarID(ncidv0,varnamev0);
datav0=netcdf.getVar(ncidv0,varidv0);
datav0=double(datav0);

datav1=a_resample;
datav1_comparison=datav1==datav0;
datav1_comparison=double(datav1_comparison);
datav1_comparison_true=find(datav1_comparison==1);
MFs_supervised=length(datav1_comparison_true(:))/length(datav1_comparison(:));

MFs_supervised



%-------------------------------------------------------------------------%
% end of case1



%-------------------------------------------------------------------------%
% case2: using slice-shaped GRF data for supervised classification (with the boundary information between the RTs) %


% % % data documentation
% % clear X_DEM*  ncidv* varnamev* varidv* datav*
% % save(['classification_' num2str(i) '_RT2_innerc_SVMrbf_withBoundary.mat'],'-v7.3')

diary off;


