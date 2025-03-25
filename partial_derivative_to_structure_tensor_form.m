function Sdata = partial_derivative_to_structure_tensor_form(pt,pt1)
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
IxI=fx(pt,pt1);            % --------- used in structure tensor code
IyI=fy(pt,pt1);

Sdata=[IxI.*IxI, IxI.*IyI;
      IxI.*IyI, IyI.*IyI];

