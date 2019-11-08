function z = eta(s,epss)
%ETA eta function.
%    ETA returns the value of the eta function, defined in (23.2.19) of [1]
%    and in (3:3:3) of [2], for real scalar argument S, approximated with a
%    relative error lower than EPSS.
%
%    For S > -1 Euler transformation [3] in van Wijngaarden implementation [4]
%    is used. Otherwise the reflection formula (23.2.6) of [1] is employed,
%    involving gamma function evaluation, except when S == 0 and in the
%    trivial zeros S == -2N.

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
%	      revised  28/2/1993

if s == NaN
   z = NaN;
  elseif s == Inf
   z = 1;
  elseif s == 0
   z = 0.5;
  elseif s < 0	&  rem(s, 2) == 0
   z = 0;
  elseif s > -1
   z = etaeul(s, epss);
  else
   twotos = 2^s;
   z = (twotos - 2) / (1 - twotos) * pi^(s - 1) * sin(s * pi / 2) *...
					    gamma(1 - s) * etaeul(1 - s, epss);
end
