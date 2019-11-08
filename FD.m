function f = fd(j, x, b)
%FD complete and incomplete Fermi-Dirac integral.
%   FD(J,X) returns the value of the complete Fermi-Dirac integral of real
%   order J and real scalar argument X.  It is a driver routine, that selects
%   FDEUL for X < 0, FDETA for X <= 1.5 and J > -2, FDPOS for X > 1.
%
%   FD(J,X,B) returns the value of the incomplete Fermi-Dirac integral of real
%   order J and real scalar arguments X and B.
%
%   M-files ETA, ETAEUL, FDETA, FDEUL, FDPOS, M1KUMM and U1KUMM are also
%   required.

% References:
%   [1] M. Goano, "Series Expansion of the Fermi-Dirac Integral F_j(x) Over the
%	Entire Domain of Real j and x", Solid-State Electronics, vol. 36,
%	n. 2, p. 217-221, 1993.
%
%   [2] J. S. Blakemore, "Approximation for Fermi-Dirac Integrals, Especially
%	the Function F_1/2(eta) Used to Describe Electron Density in a
%	Semiconductor", Solid-State Electronics, vol. 25, n. 11, p. 1067-1076,
%	1982.

%	Michele Goano, 22/1/1992 - 23/1/1992
%	      revised	3/1/1993 - 28/2/1993
%	      revised  30/1/1996 - 31/1/1996

epss = 3.0e-7;
if nargin < 2  |  nargin > 3
   disp('FD: incorrect number of arguments')
   return
  elseif nargin == 2			% Complete
   if j == NaN | x == NaN
      f = NaN;
     elseif j == Inf | x == Inf
      f = Inf;
     elseif x < 0
      f = FDeul(j, x, epss);
     elseif abs(x) <= 1.5  &  j > -2
      f = FDeta(j, x, epss);
     else
      f = FDpos(j, x, epss);
   end
  elseif nargin == 3			% Incomplete
   if j == NaN | x == NaN | b == NaN
      f = NaN;
     elseif j == Inf | x == Inf
      f = Inf;
     elseif b == Inf
      f = 0;
     elseif b == 0
      f = fd(j, x);
     elseif b < 0
      f = [];
      disp('FD: b < 0')
     else
      itmax = 100;
      sum = 0;
      wksp = zeros(1,itmax);
      bn = b;
      if b < x	     % Series involving Kummer's function M
	 ebmx = -exp(b - x);
	 enbmx = -ebmx;
	 for n = 1:itmax
	    sum_old = sum;
	    add = enbmx * M1kumm(j, bn, epss);
	    if n == 1	   % Euler-van Wijngaarden transformation
	       nterm = 1;
	       wksp(1) = add;
	       sum = 0.5 * add;
	      else
	       tmp = wksp(1);
	       wksp(1) = add;
	       for i = 1:nterm
		  dum = wksp(i + 1);
		  wksp(i + 1) = 0.5 * (wksp(i) + tmp);
		  tmp = dum;
	       end
	       if abs(wksp(nterm + 1)) <= abs(wksp(nterm))
		  sum = sum + 0.5 * wksp(nterm + 1);
		  nterm = nterm + 1;
		 else
		  sum = sum + wksp(nterm + 1);
	       end
	    end
	    if abs(sum - sum_old) <= abs(sum) * epss
	       break;
	    end
	    bn = bn + b;
	    enbmx = enbmx * ebmx;
	 end
	 if n == itmax
	    disp('FD: case b < x: ITMAX too small')
	 end
	 f = fd(j, x) - b^(j + 1) / gamma(j + 2) * (1 - sum);
	else	   % Series involving Kummer's function U
	 exmb = -exp(x - b);
	 enxmb = -exmb;
	 for n = 1:itmax
	    sum_old = sum;
	    add = enxmb * U1kumm(j + 1, bn, epss);
	    if n == 1	   % Euler-van Wijngaarden transformation
	       nterm = 1;
	       wksp(1) = add;
	       sum = 0.5 * add;
	      else
	       tmp = wksp(1);
	       wksp(1) = add;
	       for i = 1:nterm
		  dum = wksp(i + 1);
		  wksp(i + 1) = 0.5 * (wksp(i) + tmp);
		  tmp = dum;
	       end
	       if abs(wksp(nterm + 1)) <= abs(wksp(nterm))
		  sum = sum + 0.5 * wksp(nterm + 1);
		  nterm = nterm + 1;
		 else
		  sum = sum + wksp(nterm + 1);
	       end
	    end
	    if abs(sum - sum_old) <= abs(sum) * epss
	       break;
	    end
	    bn = bn + b;
	    enxmb = enxmb * exmb;
	 end
	 if n == itmax
	    disp('FD: case b >= x: ITMAX too small')
	 end
	 f = b^(j + 1) / gamma(j + 1) * sum;
      end
   end
end