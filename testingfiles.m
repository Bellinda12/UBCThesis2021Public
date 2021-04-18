
close all;

slicec = 204
slicep = 134;

PhysCT = niftiread('2190CTmask2.nii');
PhysPET = niftiread('2190mask2.nii');

PCT = PhysCT(:,:,slicec);
PPET = PhysPET(:,:,slicep);

[yCoordinatesct, xCoordinatesct] = find(PCT);
[yCoordinatespet, xCoordinatespet] = find(PPET);

Yc=max(yCoordinatesct);
yc=min(yCoordinatesct);

add=10;

y1c= yc-add;
y2c=Yc+add;

Xc=max(xCoordinatesct);
xc=min(xCoordinatesct);

x1c= xc-add;
x2c=Xc+add;
Yp=max(yCoordinatespet);
yp=min(yCoordinatespet);

add=10;

y1p= yp-add;
y2p=Yp+add;

Xp=max(xCoordinatespet);
xp=min(xCoordinatespet);

x1p= xp-add;
x2p=Xp+add;

sc=[x1c y1c x2c-x1c y2c-y1c];
sp=[x1p y1p x2p-x1p y2p-y1p];

CT = niftiread('2190CT.nii');
PET = niftiread('2190PET.nii');
CT = CT(:,:,slicec);
CT = im2double(CT);
CT = imcrop(CT, sc);
CT = imresize(CT, [30 30]);

PET = PET(:,:,slicep);
PET = im2double(PET);
PET = imcrop(PET, sp);
PET = imresize(PET, [30 30]);

Fusion = 0.3*PET+0.7*CT;

figure;
imshow(0.3*PET+0.7*CT,[],'InitialMagnification', 1000);
hold on
[B,L] = bwboundaries(P,'noholes');
for k = 1:length(B)
   boundary_res = B{k};
   plot(boundary_res(:,2), boundary_res(:,1), 'r', 'LineWidth', 2)
   hold on;
end