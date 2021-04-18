
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('k8mask1.nii');
%PhysicianMask =imrotate(PhysicianMask ,90);
C = niftiread('k8mask140.nii');
%C=imrotate(C,90);
%masked region mask2_Region-2 image2masktotal ->requires 180 rotation

%%%iteration
slicef =129+1;
slicei =115+1; %
N = slicef-slicei;

Dice = zeros(N,1);

VolumePhyscian = zeros(N,1);

% Threshold = zeros(N,1);
% 
% volumeSimulation = zeros(N,1);
% 
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

%s=[40 40 100 60];
s=[x1 y1 x2-x1 y2-y1];


I = C(:,:,slice);
grayImage = imcrop(I,s);
Im=uint8(grayImage);
% figure
% imshow(Im,[]);
% title('Original Image')


%%%physcian mask
grayPImage = imcrop(J,s);
brightPixels = grayPImage >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhyscian(i) = volumeP;
% figure
% subplot(2, 2, 1);
% imshow(grayPImage,[]);
% title('Physician Mask')

seg = Im >= 0.5;

pbi=im2bw(grayPImage,0.5); %should be binary Physican mask as binary 

% subplot(2, 2, 2);
% imshowpair(seg, pbi);
andImage = seg & grayPImage;
orImage  = seg | grayPImage;
similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
Dice(i) = similarity;
jaccardscore(i) = jaccard(pbi,seg); %similarity./(2-similarity); returns scalar only if images are binary
% title(['Dice Index = ' num2str(similarity)])

slice = slice + 1;

i= i+1;

end 

DiceScores = Dice 
JaccardScores = jaccardscore
% Simulationvolume = volumeSimulation
% Physcianvolume = VolumePhyscian
% Thresholdlist = Threshold

