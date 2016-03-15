function [cutCoeffs, cutRHS] = buildCut(cutCoeffs, cutRHS, solution, scenarioPrevious, scenarioCurrent)
% BUILDCUT Build the cuts piece by piece
%
% [cutCoeffs, cutRHS] = buildCut(cutCoeffs, cutRHS, solution, scenarioPrevious, scenarioCurrent)
%    build the cut from cutCoeffs and cutRHS, using the result from
%    solution. solution is the solution from node scenarioCurrent, and the
%    cut is intended for scenarioPrevious.

% Retrieve information
scenarioCurrentIdx = scenarioCurrent.index ; % current
proba = scenarioPrevious.transitionProba(scenarioCurrentIdx) ; % proba previous -> current

T = scenarioCurrent.model.T ;
h = scenarioCurrent.model.h ;

% # const at current stage (to partition dual multipliers)

% # var at previous stage (to know where to send cuts)
nVar = size(scenarioPrevious.model.W,2) ;
forwardIdx = scenarioPrevious.model.forwardIdx ; % We only store/update the cuts from the variables in x(t-1) that appear in x(t)
                                                 % Because variables are
                                                 % always sorted, this
                                                 % works as it
% Partition dual multipliers + boundary case
pi = solution.dualCntr ;    % The part corresponding to the constraints
sigma = solution.dualCuts ; % The part corresponding to the cuts
cutRhsOld = scenarioCurrent.model.cutRHS ; % The old cut RHS
if isempty(sigma)
    sigma = 0 ;
elseif (numel(cutRhsOld) ~= numel(sigma))
    error('Vectors cutRHSOld and cutCoefIdx must have the same length! (This error should not appear : please report this bug)');
end    
if isempty(cutCoeffs)
    cutCoeffs = zeros(1,nVar) ;
end
if isempty(cutRhsOld)
    cutRhsOld = 0 ;
end
% Compute the cut itself, and accumulate in cutCoeffs and cutRHS
cutCoeffsCompact = proba*(T'*pi)';
cutCoeffs(forwardIdx) = cutCoeffs(forwardIdx) + cutCoeffsCompact;
cutRHS = cutRHS + proba*(h'*pi + cutRhsOld'*sigma);

end