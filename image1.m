function S3=image1(t,Img,S)
% Img=imread('Lena512.png');
%   Img=im2double(Img);
% %  S=[1,2;3,4];
% %1 S1=[1,2;3,4];
% %I=imread('Lena512.png');
% %Img=im2double(I);
% %[m1,n]=size(Img);
% % hS = [0 0 0; 0 -1 0; 0 1 0];
% % hE = [0 0 0; 0 -1 1; 0 0 0];
% [fx,fy]=gradient(Img);
% % fx = imfilter(Img,hS,'conv');
% % fy = imfilter(Img,hE,'conv');   
% s11=fx.^2;
% s12=fx.*fy;
% s21=fy.*fx;
% s22=fy.^2;
% pt=215;
% pt1=215;
% % S=partial_derivative_to_structure_tensor_form(pt,pt1);
%  [s11x,s11y]=gradient(s11);
%  [s12x,s12y]=gradient(s12);
%  [s21x,s21y]=gradient(s21);
%  [s22x,s22y]=gradient(s22);
%  S1=partial_derivative_to_structure_tensor_form1(pt,pt1);
% % hS = [0 0 0; 0 -1 0; 0 1 0];
% %  hE = [0 0 0; 0 -1 1; 0 0 0];
% % % [fx,fy]=gradient(Img);
% %  fx = imfilter(Img,hS,'conv');
% %  fy = imfilter(Img,hE,'conv');
% %S2=zeros(2,2);
% [m1,n]=size(S1);
C=3.31488;
lambda=1;
m=2;
S=double(S);
kappa=30;
hN = [0 1 0; 0 -1 0; 0 0 0];
hS = [0 0 0; 0 -1 0; 0 1 0];
hE = [0 0 0; 0 -1 1; 0 0 0];
hW = [0 0 0; 1 -1 0; 0 0 0];
nablaN = imfilter(S,hN,'conv');
nablaS = imfilter(S,hS,'conv');   
nablaW = imfilter(S,hW,'conv');
nablaE = imfilter(S,hE,'conv');  
if (norm(nablaN,1)>0)
cN = 1-exp(-C/(norm(nablaN,1)/kappa)^m);
else
    cN=1;
end
if (norm(nablaS,1)>0)
cS = 1-exp(-C/(norm(nablaS,1)/kappa)^m);
else
    cS=1;
end
if (norm(nablaW,1)>0)
cW = 1-exp(-C/(norm(nablaW,1)/kappa)^m);
else
    cW=1;
end
if (norm(nablaE,1)>0)
cE = 1-exp(-C/(norm(nablaE,1)/kappa)^m);
else
    cE=1;
end
% cN = 1./(1 + (norm(nablaN,1)/kappa).^2);
% cS = 1./(1 + (norm(nablaS,1)/kappa).^2);
%             cW = 1./(1 + (norm(nablaW,1)/kappa).^2);
%             cE = 1./(1 + (norm(nablaE,1)/kappa).^2);

% cS = exp(-(norm(nablaS)/kappa)^2);
% cW = exp(-(norm(nablaW)/kappa)^2);
% cE = exp(-(norm(nablaE)/kappa)^2);
C=[0,cN,0;cW,-(cN+cS+cW+cE),cE;0,cS,0];
S3=conv2(Img,C,'same');
S3=S3(:);
end