function S2=structurtensor(t,S,S1)
% I=imread('Lena512.png');
% I=imnoise(I,'gaussian');
% Img=im2double(I);
% %  S=[1,2;3,4];
% %1 S1=[1,2;3,4];
% %I=imread('Lena512.png');
% %Img=im2double(I);
% %[m1,n]=size(Img);
% % hS = [0 0 0; 0 -1 0; 0 1 0];
% % hE = [0 0 0; 0 -1 1; 0 0 0];
% [m1,n]=size(Img);
% [fx,fy]=gradient(Img);
% s11=fx.^2;
% s12=fx.*fy;
% s21=fy.*fx;
% s22=fy.^2;
% S_0=zeros(2*m1,2*n);
% for i=1:m1
%     for j=1:n
%        S_0(i,j)=s11(i,j); 
%         
%     end
% end
% for i=1:m1
%     for j=1:n
%        S_0(i+m1,j)=s21(i,j); 
%         end
% end
% for i=1:m1
%     for j=1:n
%        S_0(i,j+n)=s12(i,j); 
%         
%     end
% end
% 
% for i=1:m1
%     for j=1:n
%        S_0(i+m1,j+n)=s22(i,j);
%         
%     end
% end

% pt=215;
% pt1=215;
% S=partial_derivative_to_structure_tensor_form(pt,pt1);
%  [s11x,s11y]=gradient(s11);
%  [s12x,s12y]=gradient(s12);
%  [s21x,s21y]=gradient(s21);
%  [s22x,s22y]=gradient(s22);
 %S1=partial_derivative_to_structure_tensor_form1(pt,pt1);
% hS = [0 0 0; 0 -1 0; 0 1 0];
%  hE = [0 0 0; 0 -1 1; 0 0 0];
% % [fx,fy]=gradient(Img);
%  fx = imfilter(Img,hS,'conv');
%  fy = imfilter(Img,hE,'conv');
%S2=zeros(2,2);
S1=double(S1);
[m1,n]=size(S1);
C=3.31488;
lambda=1;
lambda1=.1;
m=2;
kappa=1;
hN = [0 1 0; 0 -1 0; 0 0 0];
hS = [0 0 0; 0 -1 0; 0 1 0];
hE = [0 0 0; 0 -1 1; 0 0 0];
hW = [0 0 0; 1 -1 0; 0 0 0];
nablaN = imfilter(S1,hN,'conv');
nablaS = imfilter(S1,hS,'conv');   
nablaW = imfilter(S1,hW,'conv');
nablaE = imfilter(S1,hE,'conv');  
if (norm(nablaN)>0)
cN = 1-exp(-C/(norm(nablaN)/kappa)^m);
else
    cN=1;
end
if (norm(nablaS)>0)
cS = 1-exp(-C/(norm(nablaS)/kappa)^m);
else
    cS=1;
end
if (norm(nablaW)>0)
cW = 1-exp(-C/(norm(nablaW)/kappa)^m);
else
    cW=1;
end
if (norm(nablaE)>0)
cE = 1-exp(-C/(norm(nablaE)/kappa)^m);
else
    cE=1;
end
% cN = exp(-(norm(nablaN)/kappa)^2);
% cS = exp(-(norm(nablaS)/kappa)^2);
% cW = exp(-(norm(nablaW)/kappa)^2);
% cE = exp(-(norm(nablaE)/kappa)^2);
% cN1 = 1./(1 + (norm(nablaN,1)/kappa).^2);
%             cS1 = 1./(1 + (norm(nablaS,1)/kappa).^2);
%             cW1 = 1./(1 + (norm(nablaW,1)/kappa).^2);
%             cE1 = 1./(1 + (norm(nablaE,1)/kappa).^2);
%             cW = norm(nablaW,1).*(1 - (norm(nablaW,1)/kappa).^2).^2;
%             cE = norm(nablaE,1).*(1 - (norm(nablaE,1)/kappa).^2).^2;
%             cS = norm(nablaS,1).*(1 - (norm(nablaS,1)/kappa).^2).^2;
%             cN = norm(nablaN,1).*(1 - (norm(nablaN,1)/kappa).^2).^2;
%S_0=S_0(:);
C=[0,cN,0;cW,-(cN+cS+cW+cE)+1,cE;0,cS,0];
%C=[0,cN1*cN,0;cW1*cW,-(cN1*cN+cS1*cS+cW1*cW+cE1*cE)+1,cE1*cE;0,cS1*cS,0];
%S2=reshape(S1,m1,n);
S2=-S+conv2(S,C,'same');%+lambda1.*(4-cN1-cW1-cE1-cS1).*(S_0-S);
S2=S2(:);
end