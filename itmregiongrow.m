
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('k2mask2.nii');
PET = niftiread('k2.nii');

%%%iteration
slicef =184+1;
slicei =184+1; 
N = slicef-slicei;

Dice = zeros(N,1);
Hausdorff = zeros(N,1);
volumeSimulation = zeros(N,1);
jaccardscore = zeros(N,1);

%active contours
epsilon=0.001;
num_it = 100;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;

slice = slicei;

i=1;

for slice = slicei:slicef;
    
P = PhysicianMask(:,:,slice);

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

s=[x1 y1 x2-x1 y2-y1];

PET = PET(:,:,slice);
grayPET=im2double(PET);
grayPET = imcrop(grayPET,s);


[pixelCount, grayLevels] = imhist(grayPET);
IPET = im2double(grayPET);
Imax = max(IPET(:));
Imean = mean(IPET(:));
Imin = min(IPET(:));
Imode = mode(IPET(:));
BS = Imode/Imax; %background to source ratio
T = (1-(0.618*BS+0.316))*(min(IPET(:)) + max(IPET(:))); %%%%%%%%%%%%%%%%%
deltaT = 0.05;
done = false;
counter = 1;
while ~done
  savedThresholds(counter) = T;
  
  g = IPET >= T; %ITM lesion 
  Z = g.*IPET;
  numberOfPixels = numel(g);
  numberOfwhitePixels = sum(g(:));
  n=numberOfwhitePixels;
  background=sum(Z(:))/numberOfwhitePixels;
  source=max(Z(:));
  backtosour=background/source;
  V1= n*4.0728*4.0728*3; %volume of the voxel in mm^3 %pixel size = 4.0728*4.0728*3
  V=V1*0.001; % volume in mL
  Tnext = (1-(0.078/V+0.618*BS+0.316))*(max(IPET(:))); %%%%%%%%%%%%%%%%%%%%%%%
  percent=Tnext/max(IPET(:));
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
%%%physcian mask
grayPImage = imcrop(P,s);
brightPixels = grayPImage >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
volumeSimulation(i)=volumeP;

pbi=im2bw(grayPImage,0.5); %should be binary Physican mask  

[yCo, xCo] = find(g)

x=round((max(xCo)+min(xCo))./2)
y=round((max(yCo)+min(yCo))./2)

[P, J] = regionGrowing(grayPET);

seg = J;

andImage = seg & grayPImage;
orImage  = seg | grayPImage;
similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
Dice(i) = similarity;
jaccardscore(i) = jaccard(pbi,seg); %similarity./(2-similarity); returns scalar only if images are binary
Hausdorff(i) = HausdorffDist(pbi,seg);
% title(['Dice Index = ' num2str(similarity)])
figure;
imshow(grayPET,[],'InitialMagnification', 2000);
hold on
[B,L] = bwboundaries(seg,'noholes');
[Z,R] = bwboundaries(pbi,'noholes');
for k = 1:length(B)
   boundary_res = B{k};
   plot(boundary_res(:,2), boundary_res(:,1), 'r', 'LineWidth', 2)
   hold on;
end
for k = 1:length(Z)
   boundary_res = Z{k};
   plot(boundary_res(:,2), boundary_res(:,1), 'g', 'LineWidth', 2)
end
title(['Dice Index = ' num2str(similarity)])
legend('hybrid','Physician Mask')

slice = slice + 1;

i= i+1;

end 

DiceScores = Dice 
JaccardScores = jaccardscore
Simulationvolume = volumeSimulation
Hausdorff = Hausdorff

