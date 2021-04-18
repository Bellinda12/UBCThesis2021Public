%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Absolute Fuzzy Connectedness segmentation example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('9702mask2.nii');
PET = niftiread('9702PET.nii');

Dice = []; 
JaccardScores = []; 
Simulationvolume = []; 
Hausdorff = []; 
Slice = [];
Structuralsimilarity = [];
i=1;

%active contours
epsilon=0.001;
num_it = 100;
rad=5; %2, 29 and 19 are also OK.
alpha=0.01;

leng= size(PhysicianMask);
set=leng(3);

cluster = 3;

for i = 1:set 
    
P = PhysicianMask(:,:,i);
N=size(P);
P = im2double(P);
    if P == zeros(N)
        continue 
    else
        [yCoordinates, xCoordinates] = find(P);
        Slice = [Slice; i];

         Y=max(yCoordinates);
         y=min(yCoordinates);
         add=4;
         y1= y-add;
         y2=Y+add;

        X=max(xCoordinates);
        x=min(xCoordinates);

        x1= x-add;
        x2=X+add;

        s=[x1 y1 x2-x1 y2-y1];

        grayPET = PET(:,:,i);
        grayPET=im2double(grayPET);
        grayPET = imcrop(grayPET,s);
        
%         filtering and enhancement 
        grayPET=imgaussfilt(grayPET,0.15); 
%         remove background
        se = strel('disk',4);
        background = imopen(grayPET,se);
        grayPET=imsubtract (grayPET,background);
        
        

%Read image and scale to [0,1]
I=grayPET;
fprintf('Image size: %d x %d\n',size(I,1),size(I,2));

%Compute adjacency
n=3;
k1=0.1;
A=adjacency(size(I),n,k1);

%Compute affinity
k2=2;
K=affinity(I,A,k2);

maxGrayLevel = max(grayPET(:));
% Find where it occurs
[x, y] = find(grayPET == maxGrayLevel);

%Seed region
S=zeros(size(I));
S([x, y])=1;

%Show seeds overlayed on image
I_rgb=repmat(I,[1,1,3]); %make rgb image (required by imoverlay)

figure(1)
image(I_rgb)
imoverlay(S,S>0);
title('Seed region');

%Compute FC
fprintf('Computing Absolute Fuzzy Connectedness map\n');
FC=afc(S,K); %Absolute FC

%Show resulting FC-map
figure(2)
imagesc(FC,[0,1])
colorbar
title('Fuzzy connectedness map');

%Threshold value
thresh=0.82;

%Show the 0.82-connected component overlayed on original image
figure(3)
image(I_rgb)
imoverlay(FC,FC>thresh);
title(sprintf('Fuzzy connected component at level %.2f',thresh));

        %%%physcian mask
        grayPImage = imcrop(P,s);
        brightPixels = grayPImage >= 0.5;
        nBP = sum(brightPixels(:));
        volumeP = nBP*4.0728*4.0728*3*0.001;
        Simulationvolume=[Simulationvolume;volumeP];

        pbi=im2bw(grayPImage,0.5); %should be binary Physican mask 

        seg=im2bw(FC);

        similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
        Dice = [Dice; similarity];
        JaccardScores = [JaccardScores;jaccard(pbi,seg)]; %similarity./(2-similarity); returns scalar only if images are binary
        Hausdorff = [Hausdorff;HausdorffDist(pbi,seg)];
        Structuralsimilarity=[Structuralsimilarity;ssim(double(pbi),double(seg))];
        
        figure;
        imshow(grayPET,[],'InitialMagnification', 2000);
        hold on
        [B,L1]=bwboundaries(seg,'noholes');
        [Z,R] = bwboundaries(pbi,'noholes');
        for k = 1:length(B)
         boundary_res1 = B{k};
        plot(boundary_res1(:,2), boundary_res1(:,1), 'r', 'LineWidth', 2)
        hold on;
        end
        for k = 1:length(Z)
         boundary_res2 = Z{k};
        plot(boundary_res2(:,2), boundary_res2(:,1), 'g', 'LineWidth', 2)
        end
        title(['Dice Index = ' num2str(similarity)])
        legend('hybrid','ITM+AC','Region Growing','Physician Mask')
        end 
i= i+1;
end

Slice = Slice;
DiceScores = Dice;
JaccardScores = JaccardScores;
Simulationvolume = Simulationvolume;
Hausdorff = Hausdorff;
Structuralsimilarity=Structuralsimilarity;
mltable = table(Slice,DiceScores,JaccardScores,Hausdorff,Structuralsimilarity,Simulationvolume)
a=max(DiceScores);
b=min(DiceScores);
c=mean(DiceScores);
d=max(JaccardScores);
e=min(JaccardScores);
f=mean(JaccardScores);
S= [a,b,c;d,e,f]
