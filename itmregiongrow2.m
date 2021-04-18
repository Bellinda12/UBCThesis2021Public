
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('1586maskT.nii');
PET = niftiread('1586PET.nii');

Dice = []; 
JaccardScores = []; 
Simulationvolume = []; 
Hausdorff = []; 
Slice = [];
Structuralsimilarity = [];
i=1;

%active contours
epsilon=0.001; %was 0.001
num_it = 100;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;%was 0.01

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
        
        % No filtering is best 
        
%         %filtering and enhancement 
%         grayPET=imgaussfilt(grayPET,0.2); 
%         %remove background
%         se = strel('disk',5);
%         background = imopen(grayPET,se);
%         grayPET=imsubtract (grayPET,background);

        [pixelCount, grayLevels] = imhist(grayPET);
        IPET = im2double(grayPET);
        Imax = max(IPET(:));
        Imean = mean(IPET(:));
        Imin = min(IPET(:));
        Imode = mode(IPET(:));
        BS = Imode/Imax; %background to source ratio
        T = ((0.618*BS+0.32))*(min(IPET(:)) + max(IPET(:))); %%%%%%%%%%%%%%%%%
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
        Tnext = ((0.078/V+0.618*BS+0.32))*(max(IPET(:))); %%%%%%%%%%%%%%%%%%%%%%% %0.24 best if use background subtraction
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
        
        newMask = imdilate(g, true(9)); %dilate the physican mask!!!

        mask_init  = zeros(size(grayPET));
        mask_init(find(newMask)) = 1; %creation of inital contour
        filter=medfilt2(grayPET);
        seg1 = local_AC_MS(filter,mask_init,rad,alpha,num_it,epsilon);

        maxGrayLevel = max(grayPET(:));
        [labeledImage, numberOfBlobs] = bwlabel(g);% find roughly how many lesions in image
        j=1;
        init=zeros(size(grayPET));
        init(find(g))=1;
        seg2=init;
        if g == zeros(size(grayPET)) | seg1 == zeros(size(grayPET))
            [x, y] = find(grayPET == maxGrayLevel);
            [P, J] = regionGrowing(filter,[x,y]);
            J=im2bw(J,0.5);
            J = imdilate(J, true(2)); %dilate the region growing by 2
            seg2 = J; 
            seg1 = seg2;
        else
        for j = 1:numberOfBlobs
            extreme = regionprops(g,grayPET,'MaxIntensity');
            extreme=extreme(j);
            value=extreme.MaxIntensity;
            [x, y] = find(grayPET == value);
            [P, J] = regionGrowing(grayPET,[x,y]);
            seg2 = seg2+J;
            j=j+1;
        end
        end
       
        % inputs: mask1, mask2: binary images
        d1 = bwdist(seg1) - bwdist(~seg1); % dt is the distance transform
        d2 = bwdist(seg2) - bwdist(~seg2);   % ~ is the logical negation
        seg = (d1+d2) > 0;   % output
        seg = ~seg;

        %iteration to get best segmentation 
        D = dice(pbi,seg);
        
        deltaD = 0.05;
        over = false;
        numb=1;
        while ~over
        d1new = bwdist(seg) - bwdist(~seg); % dt is the distance transform
        d2new = bwdist(seg2) - bwdist(~seg2);   % ~ is the logical negation
        segnew = (d1new+d2new) > 0;   % output
        segnew = ~segnew;
        Dnext = dice(pbi,segnew);; %%%%%%%%%%%%%%%%%%%%%%%
        over = numb == 4; %abs(D - Dnext) <= deltaD*D;
        seg = segnew;
        D=Dnext;
        numb=numb+1;
        end
        
        similarity =  dice(pbi,seg); %sum(andImage) / sum(orImage); returns scalar only if images are binary
        Dice = [Dice; similarity];
        JaccardScores = [JaccardScores;jaccard(pbi,seg)]; %similarity./(2-similarity); returns scalar only if images are binary
        Hausdorff = [Hausdorff;HausdorffDist(pbi,seg)];
        Structuralsimilarity=[Structuralsimilarity;ssim(double(pbi),double(seg))];
        
        figure;
        imshow(grayPET,[],'InitialMagnification', 2000);
        hold on
        [B,L]=bwboundaries(seg,'noholes');
%         [B1,L1]=bwboundaries(seg1,'noholes');
%         [B2,L2]=bwboundaries(seg2,'noholes');
        [Z,R] = bwboundaries(pbi,'noholes');
%         for k = 1:length(B1)
%          boundary_res1 = B1{k};
%         plot(boundary_res1(:,2), boundary_res1(:,1), 'y', 'LineWidth', 2)
%         hold on;
%         end
%         for k = 1:length(B2)
%          boundary_res1 = B2{k};
%         plot(boundary_res1(:,2), boundary_res1(:,1), 'm', 'LineWidth', 2)
%         hold on;
%         end
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
a=max(DiceScores);
b=min(DiceScores);
c=mean(DiceScores);
d=max(JaccardScores);
e=min(JaccardScores);
f=mean(JaccardScores);
S= [a,b,c;d,e,f]
