function b = fdeta(j,x,epss)
%FDETA Fermi-Dirac integral.
%	FDETA(J,X,EPSS) returns the value of the Fermi-Dirac integral of real
%	order J and real scalar argument X such that |X| < pi, approximated
%	with a relative error lower than EPSS.
%
%	Taylor series expansion (4) of [1] is used, involving eta function
%	defined in (23.2.19) of [2].

% References:
%   [1] W. J. Cody and H. C. Thacher, Jr., "Rational Chebyshev Approximations
%	for Fermi-Dirac Integrals of Orders -1/2, 1/2 and 3/2", Mathematics of
%	Computation, vol. 21, n. 97, p. 30-40, 1967.
%
%   [2] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions
%	with Formulas, Graphs and Mathematical Tables", National Bureau of
%	Standards, Washington, D.C., 1965.

%	Michele Goano,	2/12/1991 - 13/3/1992
%	      revised	5/ 1/1993 - 28/2/1993
%	      revised  30/ 1/1996

if x == NaN	% Tests on the arguments
   b = NaN;
  else
   itmax = 100;
   b = 0;
   fact = 1;
   ok1 = 0;
   ok2 = 0;
   for n = 1:itmax
      add = eta(j + 2 - n, epss) * fact;
      b = b + add;
%     disp([n])
%     disp([b,add])
      if abs(add) > abs(b) * epss
	 ok1 = 0;
	 ok2 = 0;
	elseif ok1 == 0
	 ok1 = 1;
	elseif ok2 == 1
	 break;
	else
	 ok2 = 1;
      end
      fact = fact * x / n;
   end
   if n == itmax
      disp('FDETA: ITMAX too small')
   end
end
