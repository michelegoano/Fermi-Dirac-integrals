function [M,k] = m1kumm(j,x,epss)
%M1KUMM Kummer's confluent hypergeometric function M(1,2+J,-X).
%	M1KUMM(J,X,EPSS) returns the value of Kummer's confluent hypergeometric
%	function M(1,2+J,-X), defined in (13.1.2) of [1], for real scalar
%	arguments J and X, approximated with a relative error lower than EPSS.
%
%	Continued fraction representation [2] is used.	Renormalization is
%	carried out as proposed in [3].

% References:
%   [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions
%	with Formulas, Graphs and Mathematical Tables", National Bureau of
%	Standards, Washington, D.C., 1965.
%
%   [2] P. Henrici, "Applied and Computational Complex Analysis. Volume 2.
%	Special Functions-Integral Transforms-Asymptotics-Continued
%	Fractions", John Wiley & Sons, New York, 1977.
%
%   [3] W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling,
%	"Numerical Recipes. The Art of Scientific Computing", Cambridge
%	University Press, Cambridge, 1986.

%	Michele Goano, 19/1/1992 - 23/1/1992
%	      revised	3/1/1993 - 28/2/1993

if j == NaN | x == NaN	% Tests on the arguments
   M = NaN;
  elseif ~min(size(j)) | ~min(size(x))
   M = [];
  elseif x == 0
   M = 1;
  else
   itmax = 100;
%
   gold = 0.;		  %
   p1 = 1;		  %
   q1 = 1;		  % Initial conditions of the continued fraction
   p2 = j + 2;		  %
   q2 = x + j + 2;	  %
   b = j + 2;		  %
   fac = 1;		  % Initial value of the normalization factor
   for k = 1:itmax
      a = -k * x;			% Evaluation of a_(2k+1) and b_(2k+1)
      b = b + 1;			%
      p1 = (a * p1 + b * p2) * fac;
      q1 = (a * q1 + b * q2) * fac;
      a = (j + k + 1) * x;		% Evaluation of a_(2k+2) and b_(2k+2)
      b = b + 1;			%
      p2 = b * p1 + a * p2 * fac;
      q2 = b * q1 + a * q2 * fac;
      if q2 ~= 0
	 fac = 1 / q2;			% Normalization factor
	 g = p2 * fac;			% Evaluation of w_(2k+2)
	 if abs(g - gold) <= abs(g) * epss
	    break;
	 end
	 gold = g;
      end
   end
   if k == itmax
      disp('M1KUMM: ITMAX too small')
   end
   M = g;
end
