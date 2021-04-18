% This Matlab file calls two localized active contour methods proposed by 
% Shawn Lankton's IEEE TIP 2008 paper 'Localizing Region-Based Active
% Contours'.
% function local_AC_UM is the localized Chan-Vese's method see TIP 2001.
% function local_AC_MS is the localized Antony's Mean Separation method see ICCV 1999  
%
% The image used in this script are from Chunming Li (http://www.engr.uconn.edu/~cmli/)

clc;
clear all;close all;

epsilon=0.001;
num_it = 150;
rad=9; %2, 9 and 19 are also OK. 2,9,19,29
alpha=0.01;

I = niftiread('k3.nii');

PhysicianMask = niftiread('k3mask1.nii'); 
newMask = imdilate(PhysicianMask, true(11)); %dilate the physican mask!!!

%%%iteration
slicef =170+1;
slicei =144+1; %167 final slice

N = slicef-slicei;

Dice = zeros(N,1);

slice = slicei;

i=1;

for slice = slicei:slicef;

Img = I(:,:,slice);
P = newMask(:,:,slice);

%%%SUV calculation%%%
mass = 53;
height =1.59; % meters
DoseCalibrationFactor = 45015;
startT = [14, 26, 00]; %RadiopharmaceuticalStartDateTime
scanT = [15, 40, 03]; %AcquisitionTime
RadionuclideTotalDose = 221630000;
h_life = 6586.20019531250;

decay_time = 73*60+3; %in second
injected_dose = RadionuclideTotalDose; % in Bq
decay_dose = injected_dose*(2.^(-decay_time/h_life));
weight = mass; % in Kg
SUVbw = weight*1000/decay_dose;

%%%

Id = im2double(Img);
Id = Id;%*SUVbw;

[yCoordinates, xCoordinates] = find(P);
mask_init  = zeros(size(Id));
mask_init(find(P)) = 1; %creation of inital contour
seg = local_AC_MS(Id,mask_init,rad,alpha,num_it,epsilon);

pbi=im2bw(P,0.5);

Dice(i) =  dice(pbi,seg);
% 
% figure;
% imshowpair(seg, pbi);
% title(['Dice Index = ' num2str(Dice(i))]);

slice = slice + 1;

i= i+1;

end 

DiceScores = Dice 









