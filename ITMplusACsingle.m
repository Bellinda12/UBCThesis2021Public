
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = imread('mask_1248.jpg');
%PhysicianMask =imrotate(PhysicianMask ,90);
C = imread('image_1248.jpg');
%C=imrotate(C,90);
%masked region mask2_Region-2 image2masktotal ->requires 180 rotation

%%%iteration

%active contours
epsilon=0.001;
num_it = 20;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;

    
P = PhysicianMask;
J = P; %imrotate(P,180);

%0.690608 0.522099 0.055249 0.060773
% 1283: 0.348066 0.546961 0.066298 0.077348
centerx = 0.690608*181;
centery = 0.522099*181;
wid=0.055249*181;
hei=0.060773*181;

% centerx = 0.348066*181;
% centery = 0.546961*181;
% wid=0.066298*181;
% hei=0.077348*181;

x1=centerx-wid/2;
y1=centery-hei/2;


s=[x1 y1 wid hei];

I = C(:,:,2);
grayImage = imcrop(I,s);
grayImage=im2double(grayImage);
% ICT = CT;
% grayCTImage = imcrop(ICT,s);
% grayCTImage=im2double(grayCTImage);
% grayImage = grayOGImage + grayCTImage;
% figure
% figure;
% imshow(grayImage,[]);
% title('Original Image')

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
  Tnext = ((0.078/V+0.618*BS+0.316))*(max(Id(:))) %%%%%%%%%%%%%%%%%%%%%%%
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

volumeS=V;
volumeSimulation=volumeS;

%%%physcian mask
grayPImage = imcrop(J,s);
brightPixels = grayPImage >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhysician = volumeP;
% figure;
% imshow(grayPImage,[]);
% title('Physician Mask')

pbi=im2bw(grayPImage,0.5); %should be binary Physican mask as binary 

newMask = imdilate(g, true(3)); %dilate the physican mask!!!

mask_init  = zeros(size(Id));
mask_init(find(newMask)) = 1; %creation of inital contour
seg = local_AC_MS(Id,mask_init,rad,alpha,num_it,epsilon);

% figure;
% imshowpair(seg, pbi);
andImage = seg & grayPImage;
orImage  = seg | grayPImage;
similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
Dice = similarity
jaccardscore= jaccard(pbi,seg) %similarity./(2-similarity); returns scalar only if images are binary
% title(['Dice Index = ' num2str(similarity)])


% figure;
% imshowpair(g, pbi);
% andImage = g & grayPImage;
% orImage  = g | grayPImage;
% similarity =  dice(pbi,g); %sum(andImage) / sum(orImage); returns scalar only if images are binary
% Dice = similarity
% jaccardscore= jaccard(pbi,g) %similarity./(2-similarity); returns scalar only if images are binary
% title(['Dice Index = ' num2str(similarity)])

figure;
imshow(grayImage,[],'InitialMagnification', 1000);
hold on
[B,L] = bwboundaries(g,'noholes');
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
title(['Dice Index = ' num2str(similarity)])
legend('ITM+AC','Physician Mask')

