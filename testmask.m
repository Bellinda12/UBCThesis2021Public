PhysicianMask = dicomread('949141002mask.dcm');

P = PhysicianMask(:,:,82);

imshow(P)