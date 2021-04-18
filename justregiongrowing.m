
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

leng= size(PhysicianMask);
set=leng(3);

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
         add=10;
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
        
                % No filtering is best 
        
        %filtering and enhancement 
        grayPET=imgaussfilt(grayPET,0.2); 
        %remove background
        se = strel('disk',5);
        background = imopen(grayPET,se);
        grayPET=imsubtract (grayPET,background);
        grayPImage = imcrop(P,s);
        
        brightPixels = grayPImage >= 0.5;
        nBP = sum(brightPixels(:));
        volumeP = nBP*4.0728*4.0728*3*0.001;

        pbi=im2bw(grayPImage,0.5); %should be binary Physican mask  

        maxGrayLevel = max(grayPET(:));
        % Find where it occurs
        [x, y] = find(grayPET == maxGrayLevel);

        filter=medfilt2(grayPET);
%         [P, J] = regionGrowing(filter,[x,y]);
% 
%         seg = J;

        [P, J] = regiongrowchris2(filter,[x,y]);

        seg = J;
        
        segn = sum(seg(:));
        volumeseg = segn*4.0728*4.0728*3*0.001;
        Simulationvolume=[Simulationvolume;volumeseg];

        similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
        Dice = [Dice; similarity];
        JaccardScores = [JaccardScores;jaccard(pbi,seg)]; %similarity./(2-similarity); returns scalar only if images are binary
        Hausdorff = [Hausdorff;HausdorffDist(pbi,seg)];
        Structuralsimilarity=[Structuralsimilarity;ssim(double(pbi),double(seg))];
        figure;
        imshow(grayPET,[],'InitialMagnification', 2000);
        hold on
        [B,L] = bwboundaries(seg,'noholes');
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
        legend('hybrid','Physician Mask')
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
meanDice=mean(DiceScores)
