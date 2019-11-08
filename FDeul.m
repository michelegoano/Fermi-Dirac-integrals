function b = fdeul(j, x, epss)
%FDEUL	Fermi-Dirac integral.
%	FDEUL(J,X,EPSS) returns the value of the Fermi-Dirac integral [1] of
%	real order J and negative scalar argument X, approximated with a
%	relative error lower than EPSS.
%
%	Euler transformation [2] in van Wijngaarden implementation [3] is used
%	to sum the alternating series (13) of [1].

% References:
%   [1] J. S. Blakemore, "Approximation for Fermi-Dirac Integrals, Especially
%	the Function F_1/2(eta) Used to Describe Electron Density in a
%	Semiconductor", Solid-State Electronics, vol. 25, n. 11, p. 1067-1076,
%	1982.
%
%   [2] A. Ralston, P. Rabinowitz, "A First Course in Numerical Analysis",
%	McGraw-Hill Book Company, New York, 1978 (for the theory of Euler
%	transformation and further references).
%
%   [3] W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling,
%	"Numerical Recipes. The Art of Scientific Computing", Cambridge
%	University Press, Cambridge, 1986 (for the implementation of Euler
%	transformation).

%	Michele Goano, 30/11/1991 - 13/3/1992
%	      revised	5/ 1/1993 - 28/2/1993

if x == NaN	% Tests on the arguments
   b = NaN;
  elseif x > 0
   b = [];
  else
   itmax = 100;
   sum = 0;
   wksp = zeros(1,itmax);
   ex = -exp(x);
   enx = -ex;
   for n = 1:itmax
      sum_old = sum;
      add = enx / n^(j + 1);
      if n == 1 	 % Euler-van Wijngaarden transformation
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
%     disp([n,nterm])
%     disp([sum,add])
      if abs(sum - sum_old) <= abs(sum) * epss
	 break;
      end
      enx = enx * ex;
   end
   if n == itmax
      disp('FDEUL: ITMAX too small')
   end
   b = sum;
end
