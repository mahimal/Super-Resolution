clear all;
clc;
I=imread('Lena512.png');
%I=imread('mandrill.png');
%I=rgb2gray(I);
z=im2double(I);
H = fspecial('gaussian',[3 3],1);
D = 4;
noise_level = 10/255;
%rng(0)
%calculate the observed image
sigma1=1;
sigma2=.2;
sigma3=.02;
sigma4=.2;
y = imfilter(z,H,'circular');
y = downsample2(y,D);
y111=y+ noise_level*randn(size(y));
[fx1,fy1]=gradient(y111);
s1=sqrt(fx1.^2+fy1.^2);
s1=upsample2(s1,4);
%imshow(y);
%y111=imresize(y11_1,[512,512]);
HR=imresize(y111,4,'bicubic');
%imshow(HR)
[m,n1]=size(HR);
% [fx,fy]=gradient(Img);
Iterations=200;
sigma11=20;
Mu=20;
t=0;
dtime=1;
tf=2.5;
while(t<tf)
[fx,fy]=gradient(HR);
R=[fx,fy];
Fext=GVFOptimizeImageForces2D(fx,fy, Mu, Iterations, sigma11);
V_gvf=[Fext(:,:,1),Fext(:,:,2)];
y=imfilter(HR,H,'circular');
y=downsample2(y,D);
y11=y+ noise_level*randn(size(y));
y1_1=(y11-y111);
y1_2=upsample2(y1_1,4);
y4=imfilter(y1_2,transpose(H),'circular');
Z=ANLST(HR,fx,fy);
Z=im2double(Z);
s2=sqrt(fx.^2+fy.^2);
y2=s1-s2;
[f1xx,f1xy]=gradient(fx);
[f1xy,f1yy]=gradient(fy);
L1=f1xx.*fy.^2-2.*f1xy.*fx.*fy+f1yy.*fx.^2;
L11=L1./(fx.^2+fy.^2);
L2=f1xx.*fx.^2+2.*f1xy.*fx.*fy+f1yy.*fy.^2;
L22=L2./(fx.^2+fy.^2);
% C3=[0,0,0;1,-2,1;0,0,0];
% C4=[0,1,0;0,-2,0;0,1,0];
% L1=conv2(HR,C4,'same');
% L=conv2(HR,C3,'same');
L=sign(imgaussfilt(L22)).*y2.*sqrt(mtimes(V_gvf,transpose(R)));
tnext = min([tf,t+dtime]);
[t,X]=ode45('lastmain2',[0:tnext/2:tnext],HR,V_gvf,R,y4,L,L11,Z,fx,fy);
ly = size(X,1);
HR = reshape(X(ly,:),m,n1);
fprintf('\rIteration %d\n',t);
t=tnext;
end
 figure, subplot 121, imshow(HR,[]), subplot 122, imshow(y111,[])

