
clc
close all     
clear all   
%% test a gray image 
I = niftiread('263777599.nii');
PhysicianMask = niftiread('263777599mask.nii');

slicef =97+1;
slicei =85+1; 

N = slicef-slicei;

slice = slicei;

Dice = zeros(N,1);
jaccardscore = zeros(N,1);
i=1;

%% parameters
cluster=10; % the number of clustering centers
se=1; % the parameter of structuing element used for morphological reconstruction
w_size=3; % the size of fitlering window

epsilon=0.001;
num_it = 150;
rad=9; %2, 29 and 19 are also OK.
alpha=0.01;

for slice = slicei:slicef
    
P = PhysicianMask(:,:,slice);

newMask = imdilate(P, true(6)); %dilate the physican mask!!!

Img = I(:,:,slice);
Img = im2double(Img);

P = im2bw(P,0.5); %ground truth
[yCoordinates, xCoordinates] = find(newMask);
%% create the square
Y=max(yCoordinates);
y=min(yCoordinates);

c=10;
d=c;

y1= y-c;
y2=Y+c;

X=max(xCoordinates);
x=min(xCoordinates);

x1= x-d;
x2=X+d;

size1 = [x1 y1 x2-x1 y2-y1];

grayImage = imcrop(Img,size1); %crop image closer to lesion 
P = imcrop(P,size1); 
fn=medfilt2(grayImage);
%% segment an image corrupted by noise

[center1,U1,~,t1]=FRFCM(double(fn),cluster,se,w_size);
f_seg=fcm_image(fn,U1,center1);
f_seg=im2bw(f_seg,0.5);
f_seg = imcomplement(f_seg);

brightPixels = P >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhyscian(i) = volumeP;

newseg = imdilate(f_seg, true(11));

mask_init = zeros(size(fn));
mask_init(find(newseg)) = 1; %creation of inital contour
seg = local_AC_MS(fn,mask_init,rad,alpha,num_it,epsilon);

similarity =  dice(P,seg);
Dice(i) = similarity;
jaccardscore(i) = jaccard(P,seg);

figure;
imshow(grayImage,[],'InitialMagnification', 2000);
hold on
[B,L] = bwboundaries(seg,'noholes');
[S,R] = bwboundaries(P,'noholes');
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
legend('FRFCM+AC','Physician Mask')

slice = slice+1;
i=i+1;
end 

Dice 
jaccardscore
% VolumePhyscian = VolumePhyscian
