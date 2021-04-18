
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('k1mask2.nii');
C = niftiread('k1.nii');
CT = niftiread('k1CT.nii');
%masked region mask2_Region-2 image2masktotal ->requires 180 rotation

slicef =136+1;
slicei =132+1; %167 final slice
N = slicef-slicei;

Dice = zeros(N,1);

VolumePhyscian = zeros(N,1);
% 
Threshold = zeros(N,1);

volumeSimulation = zeros(N,1);

jaccardscore = zeros(N,1);


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

%size=[x 40 x2-x1 y2-y1];
size=[x1 y1 x2-x1 y2-y1];
%size= [45 65 69 65]

I = C(:,:,slice);
grayOGImage = imcrop(I,size);
grayOGImage=im2double(grayOGImage);
ICT = CT(:,:,slice);
grayCTImage = imcrop(ICT,size);
grayCTImage=im2double(grayCTImage);

% grayImage = wfusimg(grayOGImage,grayCTImage,'sym4',5,'max','max');
grayImage = 0.3*grayOGImage + 0.7*grayCTImage;
grayImage=im2double(grayImage);
figure
imshow(grayImage,[],'InitialMagnification', 1000);
title('Original Image')

%%%SUV calculation%%%
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

%%%%%%

% title('Original Image', 'FontSize', fontSize);
% Let's compute and display the histogram.
[pixelCount, grayLevels] = imhist(grayImage);
Id = im2double(grayImage);
Id = Id;%*SUVbw;
Imax = max(Id(:));
Imean = mean(Id(:));
Imin = min(Id(:));
Imode = mode(Id(:));
BS = Imode/Imax;
T = (1-(0.618*BS+0.316))*(min(Id(:)) + max(Id(:))); %%%%%%%%%%%%%%%%%
deltaT = 0.05;
done = false;
counter = 1;
while ~done
  savedThresholds(counter) = T;
  %savedVolumes(counter) = V;
  
  g = Id >= T;
  Z = g.*Id;
  numberOfPixels = numel(g);
  numberOfwhitePixels = sum(g(:));
  n=numberOfwhitePixels;
  background=sum(Z(:))/numberOfwhitePixels;
  source=max(Z(:));
  backtosour=background/source;
  V1= n*4.0728*4.0728*3; %volume of the voxel in mm^3 %pixel size = 4.0728*4.0728*3
  V=V1*0.001; % volume in mL
  Tnext = (1-(0.078/V+0.618*BS+0.316))*(max(Id(:))); %%%%%%%%%%%%%%%%%%%%%%%
  percent=Tnext/max(Id(:));
  done = counter ==1; %abs(T - Tnext) < deltaT*T %counter ==12; %
  T = Tnext;
  
  % Display g
  subplot(2, 1, 2);
  imshow(g);
  title('Binary Image');
 
  counter = counter + 1;
end

FinalT= percent;
nF=n;
Threshold(i) = FinalT;

volumeS=V;
volumeSimulation(i)=volumeS;

grayPImage = imcrop(J,size);
brightPixels = grayPImage >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhyscian(i) = volumeP;
subplot(2, 2, 1);
imshow(grayPImage,[]);
title('Physician Mask')

%gbi=im2bw(g,0.5); %should be binary
pbi=im2bw(grayPImage,0.5); %should be binary 

subplot(2, 2, 2);
imshowpair(g, pbi);
andImage = g & grayPImage;
orImage  = g | grayPImage;
similarity =  dice(pbi,g); %sum(andImage) / sum(orImage); returns scalar only if images are binary
Dice(i) = similarity;
jaccardscore(i) = jaccard(pbi,g); %similarity./(2-similarity); returns scalar only if images are binary
title(['Dice Index = ' num2str(similarity)])

slice = slice + 1;

i= i+1;

end 

DiceScores = Dice 
JaccardScores = jaccardscore
Simulationvolume = volumeSimulation
Physcianvolume = VolumePhyscian
Thresholdlist = Threshold

