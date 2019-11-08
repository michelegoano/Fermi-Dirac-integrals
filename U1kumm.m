function [b,n] = u1kumm(a,x,epss)
%U1KUMM Kummer's confluent hypergeometric function U(1,1+A,X).
%	U1KUMM(A,X,EPSS) returns the value of Kummer's confluent hypergeometric
%	function U(1,1+A,X), defined in (13.1.3) of [1], for real scalar
%	arguments A and X, approximated with a relative error lower than EPSS.
%
%	The relation with incomplete Gamma function is exploited, by means of
%	(13.6.28) and (13.1.29) of [1].  For A < 0 an expansion in terms of
%	Laguerre polynomials is used [2].  Otherwise the recipe of [3] is
%	followed: series expansion (6.5.29) of [1] if X < A+1, continued
%	fraction (6.5.31) of [1] if X > A+1.
%
%	Notes: (1) EXP(-X)*X^A*U1KUMM(A,X) = GAMMA(A)*(1-GAMMAI(A,X)), where
%		   GAMMA and GAMMAI are the standard Matlab function (limited
%		   to A > 0).
%	       (2) Two bugs in the continued-fraction implementation of GAMMAI
%		   distributed with Matlab have been fixed.

% References:
%   [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions
%	with Formulas, Graphs and Mathematical Tables", National Bureau of
%	Standards, Washington, D.C., 1965.
%
%   [2] P. Henrici, "Computational Analysis with the HP-25 Pocket Calculator",
%	John Wiley & Sons, New York, 1977.
%
%   [3] W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling,
%	"Numerical Recipes. The Art of Scientific Computing", Cambridge
%	University Press, Cambridge, 1986.

%	Michele Goano, 5/1/1992 - 13/3/1992
%	      revised  5/1/1993 - 28/2/1993

if a == NaN | x == NaN	% Tests on the arguments
   b = NaN;
  elseif ~min(size(a)) | ~min(size(x))
   b = [];
  elseif x == 0
   b = gamma(a) * exp(x) / x^a;
  else
   itmax = 100;
   if a <= 0	% Laguerre polynomials
      l_n2 = 0;
      l_n1 = 1;
      gam = 1;
      b = 0;
      for n = 1:itmax
	 l_n= ((n - a - 1) * (l_n1 - l_n2) + (n + x) * l_n1) / n;
	 t = gam / (l_n1 * l_n);
	 b = b + t;
	 if abs(t) < abs(b)*epss
	    break;
	 end
	 gam = gam * (n - a) / (n + 1);
	 l_n2 = l_n1;
	 l_n1 = l_n;
      end
      if n == itmax
	 disp('U1KUMM: Laguerre polynomials method: ITMAX too small')
      end
     elseif x < a + 1	% Series expansion
      ap = a;
      b = 1 / a;
      del = b;
      for n = 1:itmax
	 ap = ap + 1;
	 del = del * x / ap;
	 b = b + del;
	 if abs(del) < abs(b) * epss
	    break;
	 end
      end
      b = gamma(a) * exp(x) / x^a  - b;
      if n == itmax
	 disp('U1KUMM: Series expansion method: ITMAX too small')
      end
     else	% Continued-fraction
      gold = 0.;
      a0 = 1;
      a1 = x;
      b0 = 0;
      b1 = 1;
      fac = 1;
      for n = 1:itmax
	 ana = n - a;
	 a0 = (a1 + a0 * ana) * fac;
	 b0 = (b1 + b0 * ana) * fac;
	 anf = n * fac;
	 a1 = x * a0 + anf * a1;
	 b1 = x * b0 + anf * b1;
	 if a1 ~= 0
	    fac = 1 / a1;
	    b = b1 * fac;
	    if abs(b - gold) / b < epss
	       break;
	    end
	    gold = b;
	 end
      end
      if n == itmax
	 disp('U1KUMM: Continued-fraction method: ITMAX too small')
      end
   end
end
