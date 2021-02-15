% Calculate T1 maps from variable flip angle 3D SPGR data
% Using a lnear model with fixed TR and TE
% Daniel Bulte, University of Oxford, November 2016
% Edited by E Bluemke 2019/20

% Note: you will need to download: DICOM series reader https://uk.mathworks.com/matlabcentral/fileexchange/71807-dicom-series-reader

%%
clear all
close all

%%%%%%%%%%%%%%%%%%%
disp('Select SPGR Folder set')
dirName = uigetdir(); 
options = struct('recursive', true, 'verbose', true, 'loadCache', false);
[partitions, meta] = readDicomSeries(dirName, options);
 % Return values:
%   imagePartitions: Array of structs containing all partitions found
%   metaFilenames: Cell array of dicom filenames that contain no images

% Read image by partition index
% readDicomSeriesImage reads a dicom image (and optional dicominfo) from a
% dicom series partition that was found with the readDicomSeries function.
%
% This function can be used in two ways:
% [image, info] = readDicomSeriesImage(directory, partition):
% Reads a partition specified by the caller. partition should be one
% element of the imagePartitions structure returned by readDicomSeries.
% directory should be the same directory that was provided to
% readDicomSeries.
%
% The image return value will contain only the frames specified in the
% partition, typically in a 3D matrix. The type is the same as returned
% from dicomread (usually int16).
%
% The info return value is either a dicominfo structure in case of an
% enhanced dicom file, or a cell array containing dicominfo structures in
% case of a series of classic dicom files.
[image1, info1] = readDicomSeriesImage(dirName, partitions(1));


nbrow = size(image1,1);
nbcol = size(image1,2);
nbslice = size(image1, 3);
nbseries = length(partitions);

data = zeros(nbrow,nbcol,nbslice,nbseries); % Complex data


for k = 1:nbseries
    [image, info] = readDicomSeriesImage(dirName, partitions(k));
	dataTmp = image;
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice 
		data(:,:,ss,k) = dataTmp(:,:,ss); 
    end
end 
size(data)

vfa_nii_1=load_nii('ZZ_TEST_OEMR_VFA2deg_20210203103845_32_reg.nii.gz');
vfa_nii_1_img=fliplr(rot90(vfa_nii_1.img,1));
data(:,:,:,1)=vfa_nii_1_img

vfa_nii_2=load_nii('ZZ_TEST_OEMR_VFA5deg_20210203103845_31_reg.nii.gz');
vfa_nii_2_img=fliplr(rot90(vfa_nii_2.img,1));
data(:,:,:,2)=vfa_nii_2_img

vfa_nii_3=load_nii('ZZ_TEST_OEMR_VFA10deg_20210203103845_30_reg.nii.gz');
vfa_nii_3_img=fliplr(rot90(vfa_nii_3.img,1));
data(:,:,:,3)=vfa_nii_3_img

vfa_nii_4=load_nii('ZZ_TEST_OEMR_VFA15degAX_20210203103845_28_reg.nii.gz');
vfa_nii_4_img=fliplr(rot90(vfa_nii_4.img,1));
data(:,:,:,4)=vfa_nii_4_img

vfa_nii_5=load_nii('ZZ_TEST_OEMR_VFA20deg_AX_20210203103845_29.nii.gz');
vfa_nii_5_img=fliplr(rot90(vfa_nii_5.img,1));
data(:,:,:,5)=vfa_nii_5_img

% data should now be 256x256x14x4

figure
imagesc(data(:,:,nbslice/2));
name='Raw loaded SPGR'
caxis([0,3500]);
axis off
title(sprintf('%s.png',name))
colorbar
daspect([1 1 1])
saveas(gcf,sprintf('%s.png',name))

% flip_ang = [2,5,10,15, 20]; degrees
flip_ang = [0.0349066;0.0872665;0.174533;0.261799;0.349066]; % radians
%flip_ang = [0.349066;0.261799;0.174533;0.0872665;0.0349066]; % radians flipped
sinfa = sin(flip_ang);
tanfa = tan(flip_ang);

sfa1 = sinfa(1);
sfa2 = sinfa(2);
sfa3 = sinfa(3);
sfa4 = sinfa(4);
sfa5 = sinfa(5);

tfa1 = tanfa(1);
tfa2 = tanfa(2);
tfa3 = tanfa(3);
tfa4 = tanfa(4);
tfa5 = tanfa(5);

tr = 5;

% create a mask to speed up calc
% mask1 = dataTmp(:,:,:) - min(min(min(dataTmp(:,:,:))));
% mask = mask1./max(max(max(mask1)));
% mask(le(mask,0.08))=0;                     
% mask(ge(mask,0.08))=1;

% initialise matrices 
t1map = zeros(nbrow,nbcol,nbslice);
line = zeros(5,nbrow*nbcol*nbslice);

% fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[3000,10000],'StartPoint',[10,1]);

% myfittype = fittype('Eone * Si + Mo','dependent',{'y'},...
%     'independent',{'Si'},'coefficients',{'Eone','Mo'},'options',fo); % Eone * Si + Mo * (1 - Eone)


%% Calculate T1
for z=1:nbslice
    for y=1:nbcol
        for x=1:nbrow
        %if (mask(x,y,z)==1)
            line(:,x*y*z) = data(x,y,z,:);
            why = [line(1,x*y*z)./sfa1;line(2,x*y*z)./sfa2;line(3,x*y*z)./sfa3;line(4,x*y*z)./sfa4];
            echs = [line(1,x*y*z)./tfa1;line(2,x*y*z)./tfa2;line(3,x*y*z)./tfa3;line(4,x*y*z)./tfa4];
%             f = fit(echs,why,myfittype);
            f = polyfit(echs,why,1);
%             coeffvals = coeffvalues(f);
            Tone = -tr/log(f(1));
            t1map(x,y,z)= real(Tone);
            
            if (isnan(t1map(x,y,z)) || t1map(x,y,z)<0 || isinf(t1map(x,y,z)))
                t1map(x,y,z)=0;
                
            end
        %end
        end
    end
    z % counter to show how far through calc is
end

% go from slope to t1... somehow, problem is variable TR


for i=1:nbslice

figure
imagesc(t1map(:,:,i));
title(sprintf('%s-%d.png',dirName,i))
caxis([0,3500]);
colorbar
axis off
daspect([1 1 1])
saveas(gcf,sprintf('%s_%d.png',dirName,i))

end

dicomt1map = uint16(reshape(t1map,[nbrow nbcol 1 nbslice]));
% reshape undoes the squeeze, which removed the colour dimension

%try without copy and createmode
% try with lowercase copy
metadata=info1(1);
dicomwrite(dicomt1map,sprintf('%s_T1map.dcm',dirName), metadata{1,1},'CreateMode','Copy');
%niftiwrite(dicomt1map,sprintf('%s_T1map.nii',dirName),metadata{1,1});

beep % will alert you when it's done running :-)

%% end
