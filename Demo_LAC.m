% This Matlab file calls two localized active contour methods proposed by 
% Shawn Lankton's IEEE TIP 2008 paper 'Localizing Region-Based Active
% Contours'.
% function local_AC_UM is the localized Chan-Vese's method see TIP 2001.
% function local_AC_MS is the localized Antony's Mean Separation method see ICCV 1999  
%
% The image used in this script are from Chunming Li (http://www.engr.uconn.edu/~cmli/)

clc;
clear all;close all;

%masked region mask2_Region-2 image2masktotal ->requires 180 rotation
%J = imrotate(P,180); only for image2masktotal.nii

slicef =246+1;
slicei =245+1; 
N = slicef-slicei;

size=[40 40 100 100];

Dice = zeros(N,1);

epsilon=0.001;
num_it = 300;
rad=29; %2, 9 and 19 are also OK.
alpha=0.01;

slice = slicei;

i=1;

for slice = slicei:slicef;

Img = niftiread('image.nii');
PhysicianMask = niftiread('mask1_Region-1.nii'); 

Img = Img(:,:,slice);

P = PhysicianMask(:,:,slice);

[yCoordinates, xCoordinates] = find(P);
Y=max(yCoordinates);
y=min(yCoordinates);

c=10;
d=10;

y1= y-c;
y2=Y+c;

X=max(xCoordinates);
x=min(xCoordinates);

x1= x-d;
x2=X+d;

mask_init  = zeros(168,168);
mask_init(x1:x2,y1:y2) = 1;
seg = local_AC_MS(Img,mask_init,rad,alpha,num_it,epsilon);

pbi=im2bw(P,0.5);

similarity =  dice(pbi,seg);
Dice(i) = similarity;

figure;
imshowpair(seg, pbi);
title(['Dice Index = ' num2str(similarity)]);

slice = slice + 1;

i= i+1;

end

DiceScores = Dice 







