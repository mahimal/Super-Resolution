%anisotropic nonlinear structure tensor
function HR = ANLST(I,fx,fy)
%  I=imread('Lena512.png');
% I=imread('Lena512.png');
% I=imnoise(I,'gaussian');
% z=im2double(I);
% H = fspecial('gaussian',[3 3],1);
% K = 4;
% noise_level = 10/255;
% %rng(0)
% %calculate the observed image
% sigma1=1;
% sigma2=.2;
% sigma3=.02;
% sigma4=.2;
% y = imfilter(z,H,'circular');
% y = downsample2(y,K);
% y1=y+ noise_level*randn(size(y));
% %imshow(y);
% I=interp2(y1,'cubic');
%imshow(I)
%Img=[1,2,3;4,5,6;7,8,9];
Img=im2double(I);
[m1,n]=size(Img);
% hS = [0 0 0; 0 -1 0; 0 1 0];
% hE = [0 0 0; 0 -1 1; 0 0 0];
ds=.1;
dt=.1;
C=3.31488;
lambda=1;
m=2;
%  maskSize = max([m1, n]);
%  midpt = ceil(maskSize/2);
[fx,fy]=gradient(Img);
%fx = imfilter(Img,hS,'conv');
%fy = imfilter(Img,hE,'conv');   
s11=fx.^2;
s12=fx.*fy;
s21=fy.*fx;
s22=fy.^2;
S=zeros(2*m1,2*n);
for i=1:m1
    for j=1:n
       S(i,j)=s11(i,j); 
        
    end
end
for i=1:m1
    for j=1:n
       S(i+m1,j)=s21(i,j); 
        end
end
for i=1:m1
    for j=1:n
       S(i,j+n)=s12(i,j); 
        
    end
end

for i=1:m1
    for j=1:n
       S(i+m1,j+n)=s22(i,j);
        
    end
end
r=((s11-s22)./(sqrt((s22-s11)^2+4*s12^2)));
r1=((s12)./(sqrt((s22-s11)^2+4*s12^2)));
mu1=(s11+s22+sqrt((s22-s11)^2+4*s12^2));
mu2=(s11+s22-sqrt((s22-s11)^2+4*s12^2));
lambda2=1;
if (mu1>0)
    lambda1=(1-exp(-C./(mu1./lambda).^m));
else 
    lambda1=1;
end
gtens1=1./(1+lambda1);
gtens2=(1./(1+lambda2)^.5);
a=1/2.*(gtens1+gtens2+(gtens1-gtens2).*(r));
b=((gtens1-gtens2).*(r1));
c=1/2.*(gtens1+gtens2-(gtens1-gtens2).*(r));
kappa=1;
S11=[a,b;b,c];
[s11x,s11y]=gradient(s11);
[s12x,s12y]=gradient(s12);  
[s21x,s21y]=gradient(s21);
[s22x,s22y]=gradient(s22);
Ix1=s11x; 
Ix2=s12x;
Ix3=s22x;
Ix4=s21x; 
Iy1=s11y; 
Iy2=s12y;
Iy3=s22y;
Iy4=s21y;
s11_1=Ix1.*Ix1+2*Ix2.*Ix2+Ix3.*Ix3;
s12_1=Ix1.*Iy1+2*Ix2.*Iy2+Ix3.*Iy3;
s21_1=Ix1.*Iy1+2*Ix2.*Iy2+Ix3.*Iy3;
s22_1=Iy1.*Iy1+2*Iy2.*Iy2+Iy3.*Iy3;
r_1=((s11_1-s22_1)./(sqrt((s22_1-s11_1)^2+4*s12_1^2)));
r1_1=((s12_1)./(sqrt((s22_1-s11_1)^2+4*s12_1^2)));
mu1_1=(s11_1+s22_1+sqrt((s22_1-s11_1)^2+4*s12_1^2));
mu2_1=(s11_1+s22_1-sqrt((s22_1-s11_1)^2+4*s12_1^2));
lambda2_1=1;
if (mu1>0)
    lambda1_1=(1-exp(-C./(mu1_1./lambda).^m));
else 
    lambda1_1=1;
end
gtens1_1=1./(1+lambda1_1);
gtens2_1=(1./(1+lambda2_1)^.5);
a1=1/2.*(gtens1_1+gtens2_1+(gtens1_1-gtens2_1).*(r_1));
b1=((gtens1_1-gtens2_1).*(r1_1));
c1=1/2.*(gtens1_1+gtens2_1-(gtens1_1-gtens2_1).*(r_1));
S1=[a1,b1;b1,c1];
[n1,n11]=size(S1);
t=0;
dtime=.1;
tf=.4;
while(t<tf)
%ds=.1;
%C=[0,cN,0;cW,-(cN+cS+cW+cE),cE;0,cS,0];
tnext = min([tf,t+dtime]);
[t,X]=ode45('structurtensor',[0:tnext/2:tnext],S,S1);
ly = size(X,1);
Z=X(ly,:);
Z=reshape(Z,n11,n11);
S=Z;
%end
%[fx,fy]=gradient(Img);
%fx = imfilter(Img,hS,'conv');
%fy = imfilter(Img,hE,'conv'); 

E=Img;
S=double(S);

[t,F]=ode45('image1',[0:tnext/2:tnext],E,S);%H,K,sigma1,sigma2,sigma3,sigma4);
 ly = size(F,1);
 Img=F(ly,:);
 Img=reshape(Img,m1,n);
fprintf('\rIteration %d\n',t);
t=tnext;
[fx,fy]=gradient(Img);
%fx = imfilter(Img,hS,'conv');
%fy = imfilter(Img,hE,'conv');   
s11=fx.^2;
s12=fx.*fy;
s21=fy.*fx;
s22=fy.^2;
r=((s11-s22)./(sqrt((s22-s11)^2+4*s12^2)));
r1=((s12)./(sqrt((s22-s11)^2+4*s12^2)));
mu1=(s11+s22+sqrt((s22-s11)^2+4*s12^2));
mu2=(s11+s22-sqrt((s22-s11)^2+4*s12^2));
lambda2=1;
if (mu1>0)
    lambda1=(1-exp(-C./(mu1./lambda).^m));
else 
    lambda1=1;
end
gtens1=1./(1+lambda1);
gtens2=(1./(1+lambda2)^.5);
a=1/2.*(gtens1+gtens2+(gtens1-gtens2).*(r));
b=((gtens1-gtens2).*(r1));
c=1/2.*(gtens1+gtens2-(gtens1-gtens2).*(r));
S=[a,b;b,c];
[fx,fy]=gradient(Img);
s11=fx.^2;
s12=fx.*fy;
s21=fy.*fx;
s22=fy.^2;
[s11x,s11y]=gradient(s11);
[s12x,s12y]=gradient(s12);  
[s21x,s21y]=gradient(s21);
[s22x,s22y]=gradient(s22);
Ix1=s11x; 
Ix2=s12x;
Ix3=s22x;
Ix4=s21x; 
Iy1=s11y; 
Iy2=s12y;
Iy3=s22y;
Iy4=s21y;
s11_1=Ix1.*Ix1+2*Ix2.*Ix2+Ix3.*Ix3;
s12_1=Ix1.*Iy1+2*Ix2.*Iy2+Ix3.*Iy3;
s21_1=Ix1.*Iy1+2*Ix2.*Iy2+Ix3.*Iy3;
s22_1=Iy1.*Iy1+2*Iy2.*Iy2+Iy3.*Iy3;
r_1=((s11_1-s22_1)./(sqrt((s22_1-s11_1)^2+4*s12_1^2)));
r1_1=((s12_1)./(sqrt((s22_1-s11_1)^2+4*s12_1^2)));
mu1_1=(s11_1+s22_1+sqrt((s22_1-s11_1)^2+4*s12_1^2));
mu2_1=(s11_1+s22_1-sqrt((s22_1-s11_1)^2+4*s12_1^2));
lambda2_1=1;
if (mu1>0)
    lambda1_1=(1-exp(-C./(mu1_1./lambda).^m));
else 
    lambda1_1=1;
end
gtens1_1=1./(1+lambda1_1);
gtens2_1=(1./(1+lambda2_1)^.5);
a1=1/2.*(gtens1_1+gtens2_1+(gtens1_1-gtens2_1).*(r_1));
b1=((gtens1_1-gtens2_1).*(r1_1));
c1=1/2.*(gtens1_1+gtens2_1-(gtens1_1-gtens2_1).*(r_1));
S1=[a1,b1;b1,c1];
end
%  for i=1:m
%      for j=1:n
%         ST(i,j)=[s11(i,j),s12(i,j)];
%      end
% end
HR=(Img);
%image8Bit = uint8(255 * mat2gray(HR));
       % HR=pwlsig1(HR);
        %l = uint8(255 * mat2gray(HR));
        %figure, subplot 121, imshow(l,[]), subplot 122, imshow(I,[])
end