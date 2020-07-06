%
% Name :
%   julday.m
% 
% Purpose :
%   Calculates the Julian Day Number for a given day, month and year.
% 
% Calling Sequence :
%   result = julday(day, month, year)
% 
% Inputs:
%   month -  number of the desired month (1 = january, ..., 12 = december). 
%   day   -  number of day of the month. 
%   year  -  number of the desired year.
% 
%   all inputs are required to be scalar integers
% 
% Outputs:
%   julday - Julian Day Number (which begins at noon) of the specified
%            calendar date.
% 
% Common blocks:
%   None.
% 
% Dependencies:
%   None.
%
% Restrictions:
%   Accuracy using IEEE double precision numbers is approximately
%   1/10000th of a second.
% 
% Modification History:
%   Translated from "Numerical Recipies in C", by William H. Press,
%   Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
%   Cambridge University Press, 1988 (second printing).
% 
% V1.0  Manuel A. Cervera 08/10/1998
% 
% V1.1  L.H. Pederick 18/09/2014
%      Modified to handle vectorized input

function julian = julday(day, month, year)

% % make sure that all inputs are scalar
% if (length(day) > 1 | length(month) > 1 | length(year) > 1)
%   fprintf('\nall inputs are required to be scalar\n')
%   return
% end

% gregorian calender was adopted on oct. 15, 1582
greg = 15 + 31 .* (10 + 12 .* 1582);

% process the input
if (any(year == 0))
  fprintf('\nthere is no year zero\n')
  return
end

year(year < 0) = year(year < 0) + 1;

jy = zeros(size(day));
jm = zeros(size(day));

afterfeb = month > 2;

jy(afterfeb) = year(afterfeb);
jm(afterfeb) = month(afterfeb) + 1;
jy(~afterfeb) = year(~afterfeb) - 1;
jm(~afterfeb) = month(~afterfeb) + 13;

julian = floor(365.25 .* jy) + floor(30.6001 .* jm) + day + 1720995;

% test whether to change to gregorian calendar.
aftergreg = ((day + 31 .* (month + 12 .* year)) >= greg);
ja = fix(0.01 .* jy(aftergreg));
julian(aftergreg) = julian(aftergreg) + 2 - ja + fix(0.25 .* ja);

return
