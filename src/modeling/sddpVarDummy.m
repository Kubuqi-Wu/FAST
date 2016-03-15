function out = sddpVarDummy(varargin)
% varargin =
% (none) : a 1x1 sddpVar is created
% 1, 2, 5, 9 : a 1x2x5x9 sddpVar is created
% [1 2 3] : a 1x2x3 sddpVar is created
% A mixed of the above will create the corresponding array 
% This array is filled with dumy variables !
if nargin == 0
    out = sddpAffine(1,0,-1) ;
else
    % Flatten out varargin : {[1 2],3,[1 4],6} -> [1 2 3 1 4 6]
    sizes = mergeCells(varargin)' ;
    % Pre-allocate sddpAffine array
    out = repmat(sddpAffine(1,0,-1),sizes);       
end