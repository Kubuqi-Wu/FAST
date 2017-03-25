function lattice = backwardPass(lattice,solutionForward,params)
% BACKWARDPASS Run all backward passes on the lattice
%
% lattice = BACKWARDPASS(lattice, solutionForward, params)
%   runs all backward passes on the lattice lattice given the
%   solutionForward
%
%   See also SDDP, LATTICE, SDDPSETTINGS, FORWARDPASS

H = lattice.getH() ;
McCount = size(solutionForward, 1) ;
if(size(solutionForward, 2)~= 1)
    error('fast::backwardPass solutionForward should be of size [1 x McCount]') ;
end
for i=1:McCount
    if H ~= size(solutionForward{i}.solutionForwardCells, 1)
        error('fast::backwardPass solutionForward.solutionForwardCells should be of size [H x 1]') ;
    end
end

for time = H:-1:2
    % Retreive nodes at previous stage
    scenarioPreviousCells = lattice.getScenariosCells(time-1) ;
    L = length(scenarioPreviousCells) ;
    cutRHS = zeros(L,McCount) ;
    cutCoeffs = cell(L,McCount) ;
    for Mc = 1:McCount
        scenarioCurrent = [] ;
        while true
            scenarioCurrent = lattice.explore(scenarioCurrent,time);
            if isempty(scenarioCurrent)
                break ;
            end
            displayMessage(sprintf('%d) %d - %d backward pass', Mc, time, scenarioCurrent.getIndex()), params, 2) ;
            solutionBackward = scenarioCurrent.solve(solutionForward{Mc,1}.solutionForwardCells{time-1}, time==H, params) ;
            for idxPrevious = 1:L
                [cutCoeffs{idxPrevious, Mc}, cutRHS(idxPrevious, Mc)] = buildCut(cutCoeffs{idxPrevious, Mc}, cutRHS(idxPrevious, Mc), solutionBackward, scenarioPreviousCells{idxPrevious}, scenarioCurrent);
            end
        end
    end
    % Add cuts to the scenario at previous stage (gaffe parallel)
    for idxPrevious = 1:L
        lattice = lattice.addCuts(scenarioPreviousCells{idxPrevious}, cutCoeffs(idxPrevious,:), cutRHS(idxPrevious,:), params) ;
    end
end

end