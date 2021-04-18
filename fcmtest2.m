
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('2644mask1.nii');
PET = niftiread('2644PET.nii');

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
        siz=size(grayPET);
        
%         filtering and enhancement 
        grayPET=imgaussian(grayPET,0.2,s); 
%         remove background
        se = strel('disk',4);
        background = imopen(grayPET,se);
        grayPET=imsubtract (grayPET,background);
        
        %%%physcian mask
        grayPImage = imcrop(P,s);
        brightPixels = grayPImage >= 0.5;
        nBP = sum(brightPixels(:));
        volumeP = nBP*4.0728*4.0728*3*0.001;
        Simulationvolume=[Simulationvolume;volumeP];

        pbi=im2bw(grayPImage,0.5); %should be binary Physican mask 
%Matlab fcm function returns an indexed matrix basing on cluster number
        data = grayPET(:);
        [center,member ]= fcm( data, cluster );
        [~,~]= sort( center );
        member = member';
        [~ ,label]= max(member,[],2);
        resize=size(grayPET);
        segmROI = reshape(label,resize);
        px =[];
        py =[];
%         hold on
        in = segmROI (round(9),round(9));
        ROImask = segmROI == in ;
        
        [ n, m1] = size(grayPET);

        A1 = n*m1 ;
        s2 = regionprops(ROImask( ROImask ==1) , 'ConvexArea');
        A2 = s2.ConvexArea ;
        if A2 >0.8* A1 %select complementary image to better visualization
        ROImask = imcomplement( ROImask );
        end
        [L , nobj]= bwlabel( ROImask ,4);
        c = L ;
        se2=strel ('disk' ,1);
        %allow to give the more likely mask also when the image is very degraded
         if nobj >3
            in1 = L ( round(9) , round(9));
            if in1 ==0
            [~ , idx ]= bwdist ( L ~=0);
         d = sub2ind ( size ( L ) , round(9) , round(9));
         ddx = idx ( d );
         [ xc , yc ] = ind2sub ( size ( L ) , ddx );
            in1 = L ( xc , yc );
            end
            c = L == in1 ;
         end
%         hold on
        c = imdilate (c , se2 );
        
        seg=im2bw(c);
        
        %active contours
        
%         newMask = imdilate(g, true(6)); %dilate the physican mask!!!
% 
%         mask_init  = zeros(size(grayPET));
%         mask_init(find(newMask)) = 1; %creation of inital contour
%         filter=medfilt2(grayPET);
%         seg = local_AC_MS(filter,mask_init,rad,alpha,num_it,epsilon);
        

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
        legend('FCM','Physician Mask')
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



