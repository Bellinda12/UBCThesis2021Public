% This function is only suitabe for gray image
function gx=fcm_image(f,U,center)
[m,n]=size(f);
% [~,idx_f]=max(U);
[~,idx_f]=min(U);
% [~,idx_f]=sort(U(:), 'ascend');
input_f=reshape(idx_f,[m n]); 
input_ff=zeros(m,n); %input_ff denotes the classification result based on the cluster center
for k=1:length(center(:,1))
    t=(input_f==k).*center(k);
    input_ff=input_ff+t; 
end
gx=uint8(input_ff);