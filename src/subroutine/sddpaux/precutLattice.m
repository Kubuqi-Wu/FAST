function lattice = precutLattice(lattice,solutionForward,params)

warning('Most likely to crash')

H = lattice.H;
NSolForward = size(solutionForward,2);

for time = H:-1:2
    scenarioPreviousCells = lattice.getScenariosCells(time-1) ;
    L = length(scenarioPreviousCells) ;
    cutRHS = zeros(L,NSolForward) ;
    cutCoeffs = cell(L,NSolForward) ;
    for j = 1:NSolForward
        scenarioCurrent = [] ;
        while true
            scenarioCurrent = lattice.explore(scenarioCurrent,time);
            if isempty(scenarioCurrent)
                break ;
            end
            displayMessage(sprintf('%d) %d - %d Computing precut',j, time, scenarioCurrent.getIndex()), params, 2) ;   
            [~, constraints] = scenarioCurrent.solve(solutionForward{time-1,j}, time==H, params) ;
            for idxPrevious = 1:L
                [cutCoeffs{idxPrevious, j}, cutRHS(idxPrevious, j)] = buildCut(cutCoeffs{idxPrevious, j}, cutRHS(idxPrevious, j), constraints, scenarioPreviousCells{idxPrevious}, scenarioCurrent);
            end
        end
    end
    % Add cuts to the scenario at previous stage (gaffe parallel)
    for idxPrevious = 1:L
        lattice = lattice.addCuts(scenarioPreviousCells{idxPrevious}, cutCoeffs(idxPrevious,:), cutRHS(idxPrevious,:), params) ;
    end
end