
clc
close all     
clear all   
%% test a gray image 
I = niftiread('k8.nii');
PhysicianMask = niftiread('k8mask5.nii');

slicef =275+1;
slicei =275+1; 

N = slicef-slicei;

slice = slicei;

Dice = zeros(N,1);
VolumePhyscian= zeros(N,1);
i=1;

%% parameters
cluster=20; % the number of clustering centers
se=1; % the parameter of structuing element used for morphological reconstruction
w_size=3; % the size of fitlering window

for slice = slicei:slicef
    
P = PhysicianMask(:,:,slice);

newMask = imdilate(P, true(11)); %dilate the physican mask!!!

Img = I(:,:,slice);

P = im2bw(P,0.5);
[yCoordinates, xCoordinates] = find(newMask);
%% create the square
Y=max(yCoordinates);
y=min(yCoordinates);

c=5;
d=c;

y1= y-c;
y2=Y+c;

X=max(xCoordinates);
x=min(xCoordinates);

x1= x-d;
x2=X+d;

size = [x1 y1 x2-x1 y2-y1];

grayImage = imcrop(Img,size); %crop image closer to lesion 
P = imcrop(P,size); 
fn=medfilt2(grayImage);
%% segment an image corrupted by noise
% tic 
[center1,U1,~,t1]=FRFCM(double(fn),cluster,se,w_size);
% Time1=toc;
% disp(strcat('running time is: ',num2str(Time1)))
f_seg=fcm_image(fn,U1,center1);
f_seg=im2bw(f_seg,0.5);
f_seg = imcomplement(f_seg);
similarity =  dice(P,f_seg);
Dice(i) = similarity;

brightPixels = P >= 0.5;
nBP = sum(brightPixels(:));
volumeP = nBP*4.0728*4.0728*3*0.001;
VolumePhyscian(i) = volumeP;

%%%plot
figure;
subplot(2, 1, 1);
imshow(fn,[]),title('Original image');
figure;
subplot(2, 1, 2);
imshow(f_seg);title('segmentated result');
subplot(2, 2, 1);
imshowpair(f_seg, P);
title(['Dice Index = ' num2str(similarity)]);
slice = slice+1;
i=i+1;
end 

Dice = Dice
VolumePhyscian = VolumePhyscian
