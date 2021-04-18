
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

slice = 180;

C = niftiread('image2.nii');
I = C(:,:,slice);
PhysicianMask = niftiread('image2masktotal.nii'); %masked region mask2_Region-2 image2masktotal ->requires 180 rotation
P = PhysicianMask(:,:,slice);

grayImage = I;
figure
subplot(2, 1, 1);
imshow(grayImage,[]);
title('Original Image')

% title('Original Image', 'FontSize', fontSize);
% Let's compute and display the histogram.
[pixelCount, grayLevels] = imhist(grayImage);
Id = im2double(grayImage);
Imax = max(Id(:));
Imean = mean(Id(:));
Imin = min(Id(:));
Imode = mode(Id(:));
BS = Imode/Imax;
T = ((0.618*BS+0.316))*(min(Id(:)) + max(Id(:))); %%%%%%%%%%%%%%%%%
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
  Tnext = ((0.078/V+0.618*BS+0.316))*(max(Id(:))); %%%%%%%%%%%%%%%%%%%%%%%
  percent=Tnext/max(Id(:));
  done = counter ==10; %abs(T - Tnext) < deltaT*T %counter ==12; %
  T = Tnext;
  
  % Display g
  subplot(2, 1, 2);
  imshow(g);
  title('Binary Image');
 
  counter = counter + 1;
end

FinalT= percent;
nF=n;
Threshold= FinalT;

volumeS=V;
volumeSimulation=volumeS;

J = imrotate(P,180);
grayPImage = J;
brightPixels = grayPImage >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhyscian = volumeP;
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
Dice= similarity;
jaccardscore = jaccard(pbi,g); %similarity./(2-similarity); returns scalar only if images are binary
title(['Dice Index = ' num2str(similarity)])

epsilon=0.001;
num_it = 10;
rad=29; %2, 9 and 19 are also OK.
alpha=0.01;

Img = niftiread('image2.nii');

PhysicianMask = niftiread('mask2_Region-2.nii'); 
pbi = imdilate(g, true(11)); %dilate the physican mask!!!

Img = Img(:,:,slice);

mask_init  = zeros(size(Img));
mask_init(find(g)) = 1; %creation of inital contour
seg = local_AC_MS(Img,mask_init,rad,alpha,num_it,epsilon);

pbi=im2bw(pbi,0.5);

similarity =  dice(pbi,seg)

figure;
imshowpair(seg, pbi);
title(['Dice Index = ' num2str(similarity)]);