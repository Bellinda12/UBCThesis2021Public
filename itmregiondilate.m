
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('5014mask2.nii');
PET = niftiread('5014PET.nii');

Dice = []; 
JaccardScores = []; 
Simulationvolume = []; 
Hausdorff = []; 
Slice = [];
number = [];
i=1;

%active contours
epsilon=0.001;
num_it = 100;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;

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
        
        %filtering and enhancement 
        grayPET=imgaussfilt(grayPET,0.2); 
        %remove background
        se = strel('disk',5);
        background = imopen(grayPET,se);
        grayPET=imsubtract (grayPET,background);


        [pixelCount, grayLevels] = imhist(grayPET);
        IPET = im2double(grayPET);
        Imax = max(IPET(:));
        Imean = mean(IPET(:));
        Imin = min(IPET(:));
        Imode = mode(IPET(:));
        BS = Imode/Imax; %background to source ratio
        T = ((0.618*BS+0.316))*(min(IPET(:)) + max(IPET(:))); %%%%%%%%%%%%%%%%%
        deltaC = 0.05;
        over = false;
        counter = 1;
        while ~over
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
        Tnext = ((0.078/V+0.618*BS+0.316))*(max(IPET(:))); %%%%%%%%%%%%%%%%%%%%%%%
        percent=Tnext/max(IPET(:));
        over = counter ==10; %abs(T - Tnext) < deltaT*T %counter ==12; %
        T = Tnext;
 
        counter = counter + 1;
        end

        FinalT= percent;
        nF=n;
%%%physcian mask
        grayPImage = imcrop(P,s);
        brightPixels = grayPImage >= 0.5;
        nBP = sum(brightPixels(:));
        volumeP = nBP*4.0728*4.0728*3*0.001;
        Simulationvolume=[Simulationvolume;volumeP];

        pbi=im2bw(grayPImage,0.5); %should be binary Physican mask  

        maxGrayLevel = max(grayPET(:));
        filter=medfilt2(grayPET);

        [labeledImage, numberOfBlobs] = bwlabel(g);% find roughly how many lesions in image
        j=1;
     
        init=zeros(size(grayPET));
%         init(find(g))=1;
        seg=init;
        if g == zeros(size(grayPET))
            [x, y] = find(grayPET == maxGrayLevel);
            thresVal = double((max(grayPET(:)) - min(grayPET(:)))) * 0.2;
            [P, J] = regionGrowing(filter,[x,y],thresVal);
            seg = J; 
        else
            for j = 1:numberOfBlobs
            extreme = regionprops(g,grayPET,'MaxIntensity');
            extreme=extreme(j);
            value=extreme.MaxIntensity;
            [x, y] = find(grayPET == value);
            if value <= 0.5*maxGrayLevel
                thresVal = double((max(grayPET(:)) - min(grayPET(:)))) * 0.05;
            else
               thresVal = double((max(grayPET(:)) - min(grayPET(:)))) * 0.0;
            end
            [P, J] = regionGrowing(filter,[x,y],thresVal);
            seg = seg+J;
            j=j+1;
            end
        end
        seg=im2bw(seg,0.5);
%         seg = imdilate(seg, true(2)); %dilate the physican mask!!!

        similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
        Dice = [Dice; similarity];
        JaccardScores = [JaccardScores;jaccard(pbi,seg)]; %similarity./(2-similarity); returns scalar only if images are binary
        Hausdorff = [Hausdorff;HausdorffDist(pbi,seg)];
        
        figure;
        imshow(grayPET,[],'InitialMagnification', 2000);
        hold on
        [B,L]=bwboundaries(seg,'noholes');
        [B1,L1]=bwboundaries(g,'noholes');
        [Z,R] = bwboundaries(pbi,'noholes');
        for k = 1:length(B)
         boundary_res1 = B{k};
        plot(boundary_res1(:,2), boundary_res1(:,1), 'r', 'LineWidth', 2)
        hold on;
        end
        for k = 1:length(B1)
         boundary_res3 = B1{k};
        plot(boundary_res3(:,2), boundary_res3(:,1), 'b', 'LineWidth', 2)
        hold on;
        end
        for k = 1:length(Z)
         boundary_res2 = Z{k};
        plot(boundary_res2(:,2), boundary_res2(:,1), 'g', 'LineWidth', 2)
        hold on;
        end
        title(['Dice Index = ' num2str(similarity)])
%         legend('Region Growing','Physician Mask')
        end 
i= i+1;
end

Slice = Slice;
DiceScores = Dice;
JaccardScores = JaccardScores;
Simulationvolume = Simulationvolume;
Hausdorff = Hausdorff;
mltable = table(Slice,DiceScores,JaccardScores,Hausdorff,Simulationvolume)
a=max(DiceScores);
b=min(DiceScores);
c=mean(DiceScores);
d=max(JaccardScores);
e=min(JaccardScores);
f=mean(JaccardScores);
S= [a,b,c;d,e,f]
% meannumb=mean(number)
% maxnumb=max(number)
% minnumb=min(number)
