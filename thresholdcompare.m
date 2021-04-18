
clc;    % Clear the command window.
clear all;
close all;
clearvars;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

PhysicianMask = niftiread('k1mask2.nii');
mask40 = niftiread('k1mask2_Rel_thres40.0_Rel_thres40.0.nii');
mask= mask40(:,:,133);
imshow(mask)
%PhysicianMask =imrotate(PhysicianMask ,90);
C = niftiread('k1.nii');
CT = niftiread('k1CT.nii');
%C=imrotate(C,90);
%masked region mask2_Region-2 image2masktotal ->requires 180 rotation

%%%iteration
% slicef =136+1;
% slicei =132+1; %
% N = slicef-slicei;
% 
% Dice = zeros(N,1);
% 
% VolumePhysician = zeros(N,1);
% jaccardscore = zeros(N,1);
% 
% %active contours
% epsilon=0.001;
% num_it = 150;
% rad=9; %2, 29 and 19 are also OK.
% alpha=0.01;
% 
% %% parameters
% cluster=12; % the number of clustering centers
% se=1; % the parameter of structuing element used for morphological reconstruction
% w_size=3; % the size of fitlering window
% 
% slice = slicei;
% 
% i=1;
% 
% for slice = slicei:slicef;
%     
% P = PhysicianMask(:,:,slice);
% J = P; %imrotate(P,180);
% M = mask40(:,:,slice);
% 
% [yCoordinates, xCoordinates] = find(P);
% 
% Y=max(yCoordinates);
% y=min(yCoordinates);
% 
% add=10;
% 
% y1= y-add;
% y2=Y+add;
% 
% X=max(xCoordinates);
% x=min(xCoordinates);
% 
% x1= x-add;
% x2=X+add;
% 
% %s=[40 40 100 60];
% s=[x1 y1 x2-x1 y2-y1];
% 
% % I = C(:,:,slice);
% % grayImage = imcrop(I,s);
% % figure
% % subplot(2, 1, 1);
% % imshow(grayImage,[]);
% % title('Original Image')
% 
% I = C(:,:,slice);
% grayOGImage = imcrop(I,s);
% grayOGImage=im2double(grayOGImage);
% ICT = CT(:,:,slice);
% grayCTImage = imcrop(ICT,s);
% grayCTImage=im2double(grayCTImage);
% 
% % grayImage = wfusimg(grayOGImage,grayCTImage,'sym4',5,'max','max');
% grayImage = 0.3*grayOGImage + 0.7*grayCTImage;
% grayImage=im2double(grayImage);
% % figure
% % imshow(grayImage,[],'InitialMagnification', 1000);
% % title('Original Image')
% 
% grayPImage = imcrop(J,s);
% pbi=im2bw(grayPImage,0.5);
% 
% graymImage = imcrop(M,s);
% mbi=im2bw(graymImage,0.5);
% 
% subplot(2, 2, 2);
% imshowpair(mbi, pbi);
% andImage = mbi & grayPImage;
% orImage  = mbi | grayPImage;
% similarity =  dice(pbi,mbi); %sum(andImage) / sum(orImage); returns scalar only if images are binary
% Dice(i) = similarity;
% jaccardscore(i) = jaccard(pbi,mbi); %similarity./(2-similarity); returns scalar only if images are binary
% title(['Dice Index = ' num2str(similarity)])
% 
% figure;
% imshow(grayImage,[],'InitialMagnification', 2000);
% hold on
% [B,L] = bwboundaries(mbi,'noholes');
% [S,R] = bwboundaries(pbi,'noholes');
% for k = 1:length(B)
%    boundary_res1 = B{k};
%    plot(boundary_res1(:,2), boundary_res1(:,1), 'r', 'LineWidth', 2)
%    hold on;
% end
% for p = 1:length(S)
%    boundary_res3 = S{p};
%    plot(boundary_res3(:,2), boundary_res3(:,1), 'g', 'LineWidth', 2)
%    hold on
% end
% title(['Mask Comparison' ])
% legend('40% SUV','Physician Mask')
% 
% 
% slice = slice + 1;
% 
% i= i+1;
% 
% end 
% 
% DiceScores = Dice 
% JaccardScores = jaccardscore
% Physcianvolume = VolumePhysician


