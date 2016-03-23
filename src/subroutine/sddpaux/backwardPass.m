function lattice = backwardPass(lattice,solutionForward,params)
% BACKWARDPASS Run all backward passes on the lattice
%
% lattice = FORWARDPASS(lattice, solutionForward, params)
%   runs all backward passes on the lattice lattice given the
%   solutionForward
%
%   See also SDDP, LATTICE, SDDPSETTINGS, FORWARDPASS

H = lattice.getH() ;
McCount = size(solutionForward, 2) ;
if H ~= size(solutionForward, 1)
    error('fast::backwardPass solutionForward should be of size [H x McCount]') ;
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
            solutionBackward = scenarioCurrent.solve(solutionForward{time-1,Mc}, time==H, params) ;
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