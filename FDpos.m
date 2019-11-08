function b = fdpos(j, x, epss)
%FDPOS Fermi-Dirac integral.
%      FDPOS(J,X,EPSS) returns the value of the Fermi-Dirac integral of real
%      order J and positive scalar argument X, approximated with a relative
%      error lower than EPSS.

% References:
%   [1] M. Goano, "Series Expansion of the Fermi-Dirac Integral F_j(x) Over the
%	Entire Domain of Real j and x", Solid-State Electronics, vol. 36,
%	n. 2, p. 217-221, 1993.
%
%   [2] J. S. Blakemore, "Approximation for Fermi-Dirac Integrals, Especially
%	the Function F_1/2(eta) Used to Describe Electron Density in a
%	Semiconductor", Solid-State Electronics, vol. 25, n. 11, p. 1067-1076,
%	1982.

%	Michele Goano, 5/1/1992 - 13/3/1992
%	      revised  5/1/1993 - 28/2/1993

if x == NaN	% Tests on the arguments
   b = NaN;
  elseif x < 0
   b = [];
  else
   itmax = 100;
   wksp = zeros(1,itmax);
   sum = 0;
   sign = 1;
   xn = x;
   for n = 1:itmax
      sum_old = sum;
      add = sign * ((j + 1) * U1kumm(j + 1, xn, epss) - M1kumm(j, xn, epss));
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
      if abs(sum - sum_old) <= abs(sum)*epss
	 break;
      end
      sign = -sign;
      xn = xn + x;
   end
   if n == itmax
      disp('FDPOS: ITMAX too small')
   end
   b = x^(j + 1) / gamma(j + 2) * (1 + sum);
end
