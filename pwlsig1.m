function y = pwlsig1(x)
%PWLSIG	piecewise linear sigmoid characteristic.
%	Y = PWLSIG(X) = 0.5 * ABS(X + 1) - 0.5 * ABS(X - 1)

%y = abs(x+1)/2 - abs(x-1)/2;
%y = -(x < -1) + (-1 <= x).*(x <= 1).*x + (1 < x); 
y=1./(1+exp(-x));
%y=255-x;
end

  

