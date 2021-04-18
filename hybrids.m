
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('k2mask2.nii');
%PhysicianMask =imrotate(PhysicianMask ,90);
C = niftiread('k2.nii');
CT = niftiread('k2CT.nii');
%C=imrotate(C,90);
%masked region mask2_Region-2 image2masktotal ->requires 180 rotation

%%%iteration
slicef =180+1;
slicei =180+1; %
N = slicef-slicei;

Dice = zeros(N,1);

VolumePhysician = zeros(N,1);
jaccardscore = zeros(N,1);

%active contours
epsilon=0.001;
num_it = 150;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;

%% parameters
cluster=12; % the number of clustering centers
se=1; % the parameter of structuing element used for morphological reconstruction
w_size=3; % the size of fitlering window

slice = slicei;

i=1;

for slice = slicei:slicef;
    
P = PhysicianMask(:,:,slice);
P = imresize(P, [512 512]);
J = P; %imrotate(P,180);

[yCoordinates, xCoordinates] = find(P);

Y=max(yCoordinates);
y=min(yCoordinates);

add=10;

y1= y-add;
y2=Y+add;

X=max(xCoordinates);
x=min(xCoordinates);

x1= x-add;
x2=X+add;

%s=[40 40 100 60];
s=[x1 y1 x2-x1 y2-y1];

% I = C(:,:,slice);
% grayImage = imcrop(I,s);
% figure
% subplot(2, 1, 1);
% imshow(grayImage,[]);
% title('Original Image')

ICT = CT(:,:,slice);
grayCTImage = imcrop(ICT,s);
grayCTImage=im2double(grayCTImage);

I = C(:,:,slice);
I = imresize(I, [512 512]);
grayOGImage = imcrop(I,s);
grayOGImage=im2double(grayOGImage);
grayImage = 0.3*grayOGImage + 0.7*grayCTImage;
figure;
imshow(grayOGImage,[],'InitialMagnification', 1000);
figure;
imshow(grayCTImage,[],'InitialMagnification', 1000);
figure;
imshow(grayImage,[],'InitialMagnification', 1000);



% 
% % grayImage = wfusimg(grayOGImage,grayCTImage,'sym4',5,'max','max');
% grayImage = 0.3*grayOGImage + 0.7*grayCTImage;
% grayImage=im2double(grayImage);
% figure
% imshow(grayImage,[],'InitialMagnification', 1000);
% % title('Original Image')
% 
% %%%SUV calculation%%%
% mass = 53;
% height =1.59; % meters
% DoseCalibrationFactor = 45015;
% startT = [14, 26, 00]; %RadiopharmaceuticalStartDateTime
% scanT = [15, 40, 03]; %AcquisitionTime
% RadionuclideTotalDose = 221630000;
% h_life = 6586.20019531250;
% 
% decay_time = 73*60+3; %in second
% injected_dose = RadionuclideTotalDose; % in Bq
% decay_dose = injected_dose*(2.^(-decay_time/h_life));
% weight = mass; % in Kg
% SUVbw = weight*1000/decay_dose;
% 
% %%%%%%
% 
% [pixelCount, grayLevels] = imhist(grayImage);
% Id = im2double(grayImage);
% Id = Id;%*SUVbw;
% Imax = max(Id(:));
% Imean = mean(Id(:));
% Imin = min(Id(:));
% Imode = mode(Id(:));
% BS = Imode/Imax; %background to source ratio
% T = (1-(0.618*BS+0.316))*(min(Id(:)) + max(Id(:))); %%%%%%%%%%%%%%%%%
% deltaT = 0.05;
% done = false;
% counter = 1;
% while ~done
%   savedThresholds(counter) = T;
%   
%   g = Id >= T; %ITM lesion 
%   Z = g.*Id;
%   numberOfPixels = numel(g);
%   numberOfwhitePixels = sum(g(:));
%   n=numberOfwhitePixels;
%   background=sum(Z(:))/numberOfwhitePixels;
%   source=max(Z(:));
%   backtosour=background/source;
%   V1= n*4.0728*4.0728*3; %volume of the voxel in mm^3 %pixel size = 4.0728*4.0728*3
%   V=V1*0.001; % volume in mL
%   Tnext = (1-(0.078/V+0.618*BS+0.316))*(max(Id(:))); %%%%%%%%%%%%%%%%%%%%%%%
%   percent=Tnext/max(Id(:));
%   done = counter ==10; %abs(T - Tnext) < deltaT*T %counter ==12; %
%   T = Tnext;
%   
%   % Display g
%   subplot(2, 1, 2);
%   imshow(g);
%   title('Binary Image');
%  
%   counter = counter + 1;
% end
% 
% FinalT= percent;
% nF=n;
% Threshold(i) = FinalT;
% 
% volumeS=V;
% volumeSimulation(i)=volumeS;
% 
% %%%physcian mask
% grayPImage = imcrop(J,s);
% brightPixels = grayPImage >= 0.5;
% nBP = sum(brightPixels(:));
% volumeP = nBP*4.0728*4.0728*3*0.001;
% VolumePhysician(i) = volumeP;
% subplot(2, 2, 1);
% imshow(grayPImage,[]);
% title('Physician Mask')
% 
% pbi=im2bw(grayPImage,0.5); %should be binary Physican mask as binary 
% 
% Q=round(nBP/10);
% 
% if Q < 5
%     Q=6;
% else
%     Q=Q;
% end
% 
% newMask = imdilate(g, true(11)); %dilate the physican mask!!!
% 
% mask_init  = zeros(size(Id));
% mask_init(find(newMask)) = 1; %creation of inital contour
% seg = local_AC_MS(Id,mask_init,rad,alpha,num_it,epsilon);
% 
% % subplot(2, 2, 2);
% % imshowpair(seg, pbi);
% % andImage = seg & grayPImage;
% % orImage  = seg | grayPImage;
% % similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
% % Dice(i) = similarity;
% % jaccardscore(i) = jaccard(pbi,seg); %similarity./(2-similarity); returns scalar only if images are binary
% % title(['Dice Index = ' num2str(similarity)])
%  
% fn=medfilt2(grayImage);
% %% segment an image corrupted by noise
% 
% [center1,U1,~,t1]=FRFCM(double(fn),cluster,se,w_size);
% f_seg=fcm_image(fn,U1,center1);
% f_seg=im2bw(f_seg,0.5);
% f_seg = imcomplement(f_seg);
% 
% newseg = imdilate(f_seg, true(11));
% 
% mask_initf = zeros(size(fn));
% mask_initf(find(newseg)) = 1; %creation of inital contour
% segF = local_AC_MS(fn,mask_initf,rad,alpha,num_it,epsilon);
% 
% subplot(2, 2, 2);
% imshowpair(segF, pbi);
% andImage = segF & grayPImage;
% orImage  = segF | grayPImage;
% similarity =  dice(pbi,segF); %sum(andImage) / sum(orImage); returns scalar only if images are binary
% Dice(i) = similarity;
% jaccardscore(i) = jaccard(pbi,segF); %similarity./(2-similarity); returns scalar only if images are binary
% title(['Dice Index = ' num2str(similarity)])
% 
% figure;
% imshow(grayImage,[],'InitialMagnification', 2000);
% hold on
% [B,L] = bwboundaries(g,'noholes');
% % [Y,Z] = bwboundaries(segF,'noholes');
% % [S,R] = bwboundaries(pbi,'noholes');
% for k = 1:length(B)
%    boundary_res1 = B{k};
%    plot(boundary_res1(:,2), boundary_res1(:,1), 'r', 'LineWidth', 2)
%    hold on;
% end
% % for o = 1:length(Y)
% %    boundary_res2 = Y{o};
% %    plot(boundary_res2(:,2), boundary_res2(:,1), 'y', 'LineWidth', 2)
% %    hold on;
% % end
% % for p = 1:length(S)
% %    boundary_res3 = S{p};
% %    plot(boundary_res3(:,2), boundary_res3(:,1), 'g', 'LineWidth', 2)
% %    hold on
% % end
% title(['Hybrid Comparison' ])
% % legend('ITM+AC','Physician Mask')
% % legend('1-ITM+AC','FRFCM+AC','Physician Mask')
% 
% 
% slice = slice + 1;
% 
% i= i+1;
% 
end 
% 
% DiceScores = Dice 
% JaccardScores = jaccardscore
% Physcianvolume = VolumePhysician
% 
% 
