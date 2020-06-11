% Calculate T1 maps from variable flip angle 3D SPGR data
% Using a lnear model with fixed TR and TE
% Daniel Bulte, University of Oxford, November 2016
% Edited by E Bluemke 2019/20

%%
clear all
% close all]

currpath = pwd;
% Where to find the DICOM data 
loadpath = uigetdir(currpath,'Select DICOM directory');

metapath = uigetdir(loadpath,'Select a DICOM VFA sub dir');
cd(metapath)

metafile = uigetfile('.dcm','Select a file');
metadata=dicominfo(metafile);

cd(loadpath)

d = dicomLoadAllSeries('.');
cd(currpath)

nbrow = size(d(1).imData,1);
nbcol = size(d(1).imData,2);
nbslice = size(d(1).imData, 3);
nbseries = length(d);

data = zeros(nbrow,nbcol,nbslice,nbseries); % Complex data


for k = 1:nbseries
	dataTmp = d(k).imData;
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice %#ok<ALIGN>
		data(:,:,ss,k) = dataTmp(:,:,ss); 
    end
end 

imagesc(data(:,:,nbslice/2));

% flip_ang = [2,5,10,15]; degrees
% flip_ang = [0.0349066;0.0872665;0.174533;0.261799]; % radians
flip_ang = [0.261799;0.174533;0.0872665;0.0349066]; % radians
sinfa = sin(flip_ang);
tanfa = tan(flip_ang);

sfa1 = sinfa(1);
sfa2 = sinfa(2);
sfa3 = sinfa(3);
sfa4 = sinfa(4);

tfa1 = tanfa(1);
tfa2 = tanfa(2);
tfa3 = tanfa(3);
tfa4 = tanfa(4);

tr = 4;

% create a mask to speed up calc
mask1 = dataTmp(:,:,:) - min(min(min(dataTmp(:,:,:))));
mask = mask1./max(max(max(mask1)));
mask(le(mask,0.08))=0;                     
mask(ge(mask,0.08))=1;

% initialise matrices 
t1map = zeros(nbrow,nbcol,nbslice);
line = zeros(4,nbrow*nbcol*nbslice);

% fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[3000,10000],'StartPoint',[10,1]);

% myfittype = fittype('Eone * Si + Mo','dependent',{'y'},...
%     'independent',{'Si'},'coefficients',{'Eone','Mo'},'options',fo); % Eone * Si + Mo * (1 - Eone)




%% Calculate T1
for z=1:nbslice
    for y=1:nbcol
        for x=1:nbrow
        if (mask(x,y,z)==1)
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
        end
        end
    end
    z % counter to show how far through calc is
end

% go from slope to t1... somehow, problem is variable TR

imagesc(rot90(t1map(:,:,nbslice/2)));
caxis([0,3500]);

dicomt1map = uint16(reshape(t1map,[nbrow nbcol 1 nbslice]));

%try without copy and createmode
% try with lowercase copy

dicomwrite(dicomt1map,sprintf('%s_T1map.dcm',loadpath), metadata,'CreateMode','copy');
%niftiwrite(dicomt1map,sprintf('%s_T1map.nii',loadpath),metadata);

% need to get metadata from another dicom file
% reshape undoes the squeeze, which removed the colour dimension

%% end