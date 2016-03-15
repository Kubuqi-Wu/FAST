function out = sddpVar(varargin)
% SDDPVAR Sddp variables
%
% out = sddpVar() builds one sddp variable that can be used to model
%    problems
% out = sddpVar(size) builds an array of size size filled with independant
%    variables
% out = sddpVar(n1,n2,...) builds an array of size n1 x n2 x ... filled
%    with independant variables
%
% Example
% x = sddpVar()
% y = sddpVar(3, 4, 2)
% z = sddpVar([1 2 3])
%
% See Also SDDP

% varargin =
% (none) : a 1x1 sddpVar is created
% 1, 2, 5, 9 : a 1x2x5x9 sddpVar is created
% [1 2 3] : a 1x2x3 sddpVar is created
% A mixed of the above will create the corresponding array 
if nargin == 0
    out = sddpAffine() ;
else
    dimensions = mergeCells(varargin) ;
    out = sddpAffine.createArray(dimensions) ;
end