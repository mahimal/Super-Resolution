function Sdata = partial_derivative_to_structure_tensor_form1(pt,pt1)
% partial_derivative_to_structure_tensor_form - sets up structure tensor form from Ix, Iy, ...%%%%
% 
% Example:
%   pDerData = [Ix, Iy, It];
%   S = partial_derivative_to_structure_tensor_form(pDerData);   %%% sets up the 3D structure tensor
% 
% Author: Shawn Arseneau
% Created: September 20, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pt=1:4
% pt1=1:4
%Sdata=zeros(size(pt),size(pt1));
I=imread('Lena512.png');
Img=im2double(I);
[fx,fy]=gradient(Img);
s11=fx^2;
s12=fx.*fy;
s21=fy.*fx;
s22=fy^2;
[s11x,s11y]=gradient(s11);
[s12x,s12y]=gradient(s12);
[s21x,s21y]=gradient(s21);
[s22x,s22y]=gradient(s22);
IxI1=fx(pt,pt1);            % --------- used in structure tensor code
IyI1=fy(pt,pt1);
Ix1=s11x(pt,pt1); 
Ix2=s12x(pt,pt1);
Ix3=s22x(pt,pt1);
Ix4=s21x(pt,pt1); 
Iy1=s11y(pt,pt1); 
Iy2=s12y(pt,pt1);
Iy3=s22y(pt,pt1);
Iy4=s21y(pt,pt1);
Sdata=[Ix1.*Ix1+2*Ix2.*Ix2+Ix3.*Ix3, Ix1.*Iy1+2*Ix2.*Iy2+Ix3.*Iy3;
      Ix1.*Iy1+2*Ix2.*Iy2+Ix3.*Iy3, Iy1.*Iy1+2*Iy2.*Iy2+Iy3.*Iy3];

