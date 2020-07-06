%
% Name :
%   chapman.m
% 
% Purpose:
%   Calculates ionospheric plasma frequency profile as a function of height
%   from input ionospheric layer parameters using Chapman layer profiles.
%
% Calling sequence:
%   plasma_freq = chapman(foE, hmE, ymE, foF1, hmF1, ymF1, foF2, ...
%                         hmF2, ymF2, height);
% 
% Inputs:
%   foE  - critical frequency of E layer (scaler)
%   hmE  - height of E layer (scaler)
%   ymE  - semi thickness of E layer (scaler)
%   foF1 - critical frequency of F1 layer (scaler)
%   hmF1 - height of F1 layer (scaler)
%   ymF1 - semi thickness of F1 layer (scaler)
%   foF2 - critical frequency of F2 layer (scaler)
%   hmF2 - height of F2 layer (scaler)
%   ymF2 - semi thickness of F2 layer (scaler)
%   heights - array of heights at which the plasma frequency is to be
%             calculated
%
% Outputs:
%   plasma_freq - array of calculated plasma frequencies
% 
% Modification History:
% 20/10/2005  V1.0  M.A.Cervera Author.
%
% 07/04/2008  V1.1  M.A.Cervera
%   The maximum value of each ionospheric layer is now calculated 
%
% 15/04/2008  V1.2  M.A.Cervera
%   Bug-fix: if maximum value of F1 is -ve then set F1 to zero and
%   recalculate the maximum value of remaining layers (E and F2)
%

function plasma_freq = chapman(foE, hmE, ymE, foF1, hmF1, ymF1, foF2, ...
      hmF2, ymF2, height); 

  % first calculate the ionospheric Chapman layer maximum values 
  a12 = chap_func(0.5, 2.*(hmE - hmF1)./ymF1);
  a13 = chap_func(1.0, sqrt(2).*(hmE - hmF2)./ymF2);
  a21 = chap_func(0.5, 2.*(hmF1 - hmE)./ymE);
  a23 = chap_func(1.0, sqrt(2).*(hmF1 - hmF2)./ymF2);
  a31 = chap_func(0.5, 2.*(hmF2 - hmE)./ymE);
  a32 = chap_func(0.5, 2.*(hmF2 - hmF1)./ymF1);
  A = [[1   a12 a13]; ...
       [a21 1   a23]; ...
       [a31 a32 1  ]];
  
  b = [foE.^2; foF1.^2; foF2.^2];
  
  lmv_sq = inv(A) * b;
  
  % check to see if the square of the F1 layer maximum plasma frequency has a
  % solution which is less than 0. This is not physical, in this case repeat,
  % but construct the A and b matricies with foF1 set to 0.
  if lmv_sq(2) < 0
    A(2,1) = 0;
    A(2,3) = 0;
    b(2) = 0;
    lmv_sq = inv(A) * b;
  end
  lmv_E =  sqrt(lmv_sq(1));
  lmv_F1 = sqrt(lmv_sq(2));
  lmv_F2 = sqrt(lmv_sq(3));

  % plasma frequency contribution from E layer
  plas_E_sq = lmv_E.^2 .* chap_func(0.5, 2 .* (height - hmE) ./ ymE);
  
  % plasma frequency contribution from F1 layer
  plas_F1_sq = lmv_F1.^2 .* chap_func(0.5, 2 .* (height - hmF1) ./ ymF1);
  
  % plasma frequency contribution from F2 layer
  plas_F2_sq = lmv_F2.^2 .* chap_func(1, sqrt(2) .* (height - hmF2) ./ ymF2);
      
  % total
  plasma_freq = sqrt(plas_E_sq + plas_F1_sq + plas_F2_sq);
  
 return
 
 
 %
 % define the Chapman function
 %
 function chap_f = chap_func(C, x)
   
   xt = x;
   xt(x > 30) = 30;
   xt(x < -30) = 30;
   
   chap_f = exp(C .* (1.0 - xt - exp(-xt)));
   
 return
 