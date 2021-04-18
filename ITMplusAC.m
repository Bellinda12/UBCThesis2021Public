
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
%CT = niftiread('k2CT.nii');
%C=imrotate(C,90);
%masked region mask2_Region-2 image2masktotal ->requires 180 rotation

%%%iteration
slicef =159+1;
slicei =150+1; %
N = slicef-slicei;

Dice = zeros(N,1);

VolumePhysician = zeros(N,1);

% Threshold = zeros(N,1);
% 
% volumeSimulation = zeros(N,1);
% 
jaccardscore = zeros(N,1);

%active contours
epsilon=0.001;
num_it = 150;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;



slice = slicei;

i=1;

for slice = slicei:slicef;
    
P = PhysicianMask(:,:,slice);
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

%s=[20 60 80 60];
%s=[67  46  36  34];
s=[x1 y1 x2-x1 y2-y1]

I = C(:,:,slice);
grayOGImage = imcrop(I,s);
grayOGImage=im2double(grayOGImage);
% ICT = CT(:,:,slice);
% grayCTImage = imcrop(ICT,s);
% grayCTImage=im2double(grayCTImage);
grayImage = grayOGImage;%grayOGImage + grayCTImage;
% % figure
% % subplot(2, 1, 1);
% % imshow(grayImage,[]);
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
[pixelCount, grayLevels] = imhist(grayImage);
Id = im2double(grayImage);
Id = Id;%*SUVbw;
Imax = max(Id(:));
Imean = mean(Id(:));
Imin = min(Id(:));
Imode = mode(Id(:));
BS = Imode/Imax; %background to source ratio
T = ((0.618*BS+0.316))*(min(Id(:)) + max(Id(:))); %%%%%%%%%%%%%%%%%
deltaT = 0.05;
done = false;
counter = 1;
while ~done
  savedThresholds(counter) = T;
  
  g = Id >= T; %ITM lesion 
  Z = g.*Id;
  numberOfPixels = numel(g);
  numberOfwhitePixels = sum(g(:));
  n=numberOfwhitePixels;
  background=sum(Z(:))/numberOfwhitePixels;
  source=max(Z(:));
  backtosour=background/source;
  V1= n*4.0728*4.0728*3; %volume of the voxel in mm^3 %pixel size = 4.0728*4.0728*3
  V=V1*0.001; % volume in mL
  Tnext = ((0.078/V+0.618*BS+0.316))*(max(Id(:))); %%%%%%%%%%%%%%%%%%%%%%%
  percent=Tnext/max(Id(:));
  done = counter ==10; %abs(T - Tnext) < deltaT*T %counter ==12; %
  T = Tnext;
  
%   % Display g
%   subplot(2, 1, 2);
%   imshow(g);
%   title('Binary Image');
 
  counter = counter + 1;
end

FinalT= percent;
nF=n;
Threshold(i) = FinalT;

volumeS=V;
volumeSimulation(i)=volumeS;

%%%physcian mask
grayPImage = imcrop(J,s);
brightPixels = grayPImage >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhysician(i) = volumeP;
% subplot(2, 2, 1);
% imshow(grayPImage,[]);
% title('Physician Mask')

pbi=im2bw(grayPImage,0.5); %should be binary Physican mask as binary 

Q=round(nBP/10);

if Q < 5
    Q=6;
else
    Q=Q;
end

newMask = imdilate(g, true(11)); %dilate the physican mask!!!

mask_init  = zeros(size(Id));
mask_init(find(newMask)) = 1; %creation of inital contour
seg = local_AC_MS(Id,mask_init,rad,alpha,num_it,epsilon);

% subplot(2, 2, 2);
% imshowpair(seg, pbi);
andImage = seg & grayPImage;
orImage  = seg | grayPImage;
similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
Dice(i) = similarity;
jaccardscore(i) = jaccard(pbi,seg); %similarity./(2-similarity); returns scalar only if images are binary
% title(['Dice Index = ' num2str(similarity)])

figure;
imshow(grayImage,[],'InitialMagnification', 1500);
hold on
[B,L] = bwboundaries(seg,'noholes');
[S,R] = bwboundaries(pbi,'noholes');
for k = 1:length(B)
   boundary_res = B{k};
   plot(boundary_res(:,2), boundary_res(:,1), 'r', 'LineWidth', 2)
   hold on;
end
for k = 1:length(S)
   boundary_res = S{k};
   plot(boundary_res(:,2), boundary_res(:,1), 'g', 'LineWidth', 2)
end
title('Physician Mask')
title(['Dice Index = ' num2str(similarity)])
legend('ITM+AC','Physician Mask')


slice = slice + 1;

i= i+1;

end 

DiceScores = Dice 
JaccardScores = jaccardscore
Simulationvolume = volumeSimulation
Physcianvolume = VolumePhysician
Thresholdlist = Threshold

