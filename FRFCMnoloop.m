info=dicominfo("159Pmask.dcm");
contour = dicomContours(info);
contour.ROIs
C=contour.ROIs.ContourData{1,1}{,end};

imagesc(C)










% % Please cite the paper "Tao Lei, Xiaohong Jia, Yanning Zhang, Lifeng He, Hongying Meng and Asoke K. Nandi, Significantly Fast and Robust
% % Fuzzy C-Means Clustering Algorithm Based on Morphological Reconstruction and Membership Filtering, IEEE Transactions on Fuzzy Systems,
% % DOI: 10.1109/TFUZZ.2018.2796074, 2018.2018"
% 
% % The paper is OpenAccess and you can download the paper freely from "http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8265186."
% % The code was written by Tao Lei in 2017.
% % If you have any problems, please contact me. 
% % Email address: leitao@sust.edu.cn
% 
% clc
% close all     
% clear all   
% %% test a gray image 
% I = niftiread('korea1pet.nii');
% 
% Img=I(:,:,69);
% 
% %% parameters
% cluster=20; % the number of clustering centers
% se=1; % the parameter of structuing element used for morphological reconstruction
% w_size=3; % the size of fitlering window
% %% create the square
% 
% size = [40 40 100 100];
% 
% grayImage = imcrop(Img,size); %crop image closer to lesion 
% fn=medfilt2(grayImage);
% %% segment an image corrupted by noise
% % tic 
% [center1,U1,~,t1]=FRFCM(double(fn),cluster,se,w_size);
% % Time1=toc;
% % disp(strcat('running time is: ',num2str(Time1)))
% f_seg=fcm_image(fn,U1,center1);
% % f_seg=im2bw(f_seg,0.5);
% % f_seg = imcomplement(f_seg);
% 
% %%%plot
% figure;
% subplot(2, 1, 1);
% imshow(fn,[]),title('Original image');
% subplot(2, 1, 2);
% imshow(f_seg);title('segmentated result');
