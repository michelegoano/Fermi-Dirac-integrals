function x = FDinv(j, y)
%FDinv inverse complete Fermi-Dirac integral.
%      FDinv(J, Y) returns the value of the inverse of the complete 
%      Fermi-Dirac integral of real order J and real scalar argument Y.
%      The Newton-Raphson algorithm is used to solve by iteration the 
%      equation FD(J, X) = Y.
%
%   M-files ETA, ETAEUL, FD, FDETA, FDEUL, FDPOS, M1KUMM and U1KUMM are 
%   also required.

% References:
%   [1] M. Goano, "Series Expansion of the Fermi-Dirac Integral F_j(x) 
%	Over the Entire Domain of Real j and x", Solid-State 
%	Electronics, vol. 36, n. 2, p. 217-221, 1993.
%   [2] A. Ralston, P. Rabinowitz, "A First Course in Numerical 
%       Analysis", New York: McGraw-Hill Book Company, 1978.

%	Michele Goano, 22/01/1992 - 28/03/1993
%	      revised  08/12/2007 - 11/01/2008

epss = 1.0e-6;
itmax = 100;

% initial guess

x = (y*gamma(j+2))^(1/(j+1));

% Newton-Raphson iteration

for n = 1:itmax
    x_corr = -(FD(j, x) - y) / FD(j-1, x);
    x = x + x_corr;
    if abs(x_corr) <= abs(x) * epss
       break;
    end
end
if n == itmax
   disp('FDinv:  Newton-Raphson algorithm not converging')
end