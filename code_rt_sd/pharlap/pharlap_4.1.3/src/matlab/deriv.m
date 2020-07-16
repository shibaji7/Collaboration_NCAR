function d = deriv(y,x)
%+
% NAME:
%	DERIV
%
% PURPOSE:
%	Perform numerical differentiation using 3-point, Lagrangian 
%	interpolation.
%
% CATEGORY:
%	Numerical analysis.
%
% CALLING SEQUENCE:
%	Dy = Deriv(Y)	 	%Dy(i)/di, point spacing = 1.
%	Dy = Deriv(Y, X)	%Dy/Dx, unequal point spacing.
%
% INPUTS:
%	Y:  Variable to be differentiated.
%	X:  Variable to differentiate with respect to.  If omitted, unit 
%	    spacing for Y (i.e., X(i) = i) is assumed.
%
% OPTIONAL INPUT PARAMETERS:
%	As above.
%
% OUTPUTS:
%	Returns the derivative.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% PROCEDURE:
%	See Hildebrand, Introduction to Numerical Analysis, Mc Graw
%	Hill, 1956.  Page 82.
%
% MODIFICATION HISTORY:
%	Written, DMS, Aug, 1984.
%       Translated into MATLAB from IDL, M. A. Cervera, HFRD, Jan 1996
%-
%

n = size(y); n=n(2);
if n < 3 
  error('Parameters must have at least 3 points');
  return
end

if exist('x') == 1 
  if n ~= length(x) 
    error('Vectors must have same size')
    return
  end
  d = ([y(2:n),y(1)] - [y(n),y(1:n-1)])./([x(2:n),x(1)] - [x(n),x(1:n-1)]);
  d(1) = (-3.0*y(1) + 4.0*y(2) - y(3))./(x(3)-x(1));
  d(n) = (3.0*y(n) - 4.0*y(n-1) + y(n-2))./(x(n)-x(n-2));
else
  d = ([y(2:n),y(1)] - [y(n),y(1:n-1)])/2.0;
  d(1) = (-3.0*y(1) + 4.0*y(2) - y(3))/2.0;
  d(n) = (3.*y(n) - 4.*y(n-1) + y(n-2))/2.0;
end











