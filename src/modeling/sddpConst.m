function out = sddpConst(k)
% SDDPCONST Builds a symbolic constant expression
%
% v = sddpConst(k) builds a variable v with a constant value k
%
% See also SDDPVAR
out = sddpAffine.createConstVariable(k) ;
