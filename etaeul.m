function eta = etaeul(s,epss)
%ETAEUL eta function.
%	ETAEUL(S,EPSS) returns the value of the eta function, defined in
%	(23.2.19) of [1] and in (3:3:3) of [2], for real scalar argument S,
%	approximated with a relative error lower than EPSS.
%
%	Euler transformation [3] in van Wijngaarden implementation [4] is used.

% References:
%   [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions
%	with Formulas, Graphs and Mathematical Tables", National Bureau of
%	Standards, Washington, D.C., 1965.
%
%   [2] J. Spanier and K. B. Oldham, "An Atlas of Functions", Hemisphere
%	Publishing Corporation, Washington, 1987.
%
%   [3] A. Ralston, P. Rabinowitz, "A First Course in Numerical Analysis",
%	McGraw-Hill Book Company, New York, 1978.
%
%   [4] W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling,
%	"Numerical Recipes. The Art of Scientific Computing", Cambridge
%	University Press, Cambridge, 1986.

%	Michele Goano, 1/12/1991 - 23/1/1992
%	      revised  5/ 1/1993 - 28/2/1993

if s == NaN	% Tests on the argument
   eta = NaN;
  elseif s == Inf
   eta = 1;
  else
   itmax = 100;
   eta = 0;
   wksp = zeros(1,itmax);
   sign = -1;
   for n = 1:itmax
      eta_old = eta;
      sign = -sign;
      add = sign / n^s;
      if n == 1 	 % Euler-van Wijngaarden transformation
	 nterm = 1;
	 wksp(1) = add;
	 eta = 0.5 * add;
	else
	 tmp = wksp(1);
	 wksp(1) = add;
	 for i = 1:nterm
	    dum = wksp(i + 1);
	    wksp(i + 1) = 0.5 * (wksp(i) + tmp);
	    tmp = dum;
	 end
	 if abs(wksp(nterm + 1)) <= abs(wksp(nterm))
	    eta = eta + 0.5 * wksp(nterm + 1);
	    nterm = nterm + 1;
	   else
	    eta = eta + wksp(nterm + 1);
	 end
      end
      if abs(eta - eta_old) <= abs(eta) * epss
	 break;
      end
   end
   if n == itmax
      disp('ETAEUL: ITMAX too small')
   end
end
