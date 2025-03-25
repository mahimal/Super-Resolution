 function HR2= lastmain2(t,I,V_gvf,R,y4,L,L11,Z,fx,fy)
 %I=imread('Lena512.png');
Img=im2double(I);
H=fspecial('gaussian',[3 3],1);
D=4;
%noise_level = 10/255;
sigma1=1;
sigma2=.2;
sigma3=.02;
sigma4=.2;
m=4;
kappa=1;
C_1=3.31488;
c=1;
[f1x,f1y]=gradient(Z);
[fxx,fxy]=gradient(f1x);
[fxy,fyy]=gradient(f1y);
Z1=fxx.*f1y.^2-2.*fxy.*f1x.*f1y+fyy.*f1x.^2;
Z11=Z1./(f1x.^2+f1y.^2);
Z2=fxx.*f1x.^2+2.*fxy.*f1x.*f1y+fyy.*f1y.^2;
Z22=Z2./(f1x.^2+f1y.^2);
%[fx,fy]=gradient(Img);
J11=cumtrapz(fx.*fx);
J22=cumtrapz(fy.*fy);
J21=cumtrapz(fx.*fy);
newmu1=(J11+J22+sqrt((J22-J11).^2+4.*J21.^2));
newmu2=(J11+J22-sqrt((J22-J11).^2+4.*J21.^2));
if (norm(sqrt(newmu1))>0)
c1 = 1-exp(-C_1./(norm(sqrt(newmu1))./kappa).^m);
else
    c1=1;
end
if (norm(sqrt(newmu2+newmu1))>0)
c2 = 1-exp(-C_1./(norm(sqrt(newmu2+newmu1))./kappa).^m);
else
    c2=1;
end
C=[0,sigma2.*c1,0;sigma3.*c2,-2*(c1.*sigma2+c2.*sigma3),sigma3.*c2;0,sigma2.*c1,0];
C1=[0,1,0;0,-2,0;0,1,0];
%R1=c.*conv2(I,C1,'same');
HR1=real(-sigma1.*y4-sigma4.*L+c.*sigma4.*L11+sigma2.*c1.*Z11+sigma3.*c2.*Z22);%.*sqrt(mtimes(V_gvf,transpose(R)));%+sigma_4.*X;%conv2(Z,C,'same');%-sigma1.*y4;%+%-sigma4.*y_2.*sign(imgaussfilt(Img)).*sqrt(mtimes(V_gvf,transpose(R)))+sigma_4.*X;
HR1=HR1(:);
HR2=HR1;
 end