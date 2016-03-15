function [indepvar1,indepvar2,depvar1,depvar2] = variableDependance(var1,var2)
% From two problem, output
%   - indep1 : variables in problem 1 which do not appear in problem 2.
%   - indep2 : variables in problem 2 which do not appear in problem 1.
%   - depvar : variables which appear in both problems.
%
% Returns index in input space, not in "modeling" space

n1 = length(var1);
n2 = length(var2);
[~,depvar1,depvar2] = intersect(var1,var2);
indepvar1 = setdiff(1:n1,depvar1);
indepvar2 = setdiff(1:n2,depvar2);

% Reshape
indepvar1 = indepvar1(:);
indepvar2 = indepvar2(:);
depvar1 = depvar1(:);
depvar2 = depvar2(:);
