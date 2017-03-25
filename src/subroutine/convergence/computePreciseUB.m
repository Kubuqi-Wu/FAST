function [converged, solutionForwardPrecise] = computePreciseUB(lattice, resultsFolder, iteration, params)

% Precise
UBMC = params.precise.count  ;
displayMessage(sprintf('Computing precise Mean Forward Cost wiht %d Monte-Carlo samples.',UBMC),params,1) ;
H = lattice.H ;
startTimePrecise = tic() ;
solutionForwardPrecise = cell(UBMC,1);            
% Computation   
for McPrecise = 1:UBMC
    [~,~,~,solutionForwardPrecise{McPrecise,1}] = forwardPass(lattice,'random',params) ;    
end                
% Convergence check  
[converged, lowerBound, meanCost, stds] = checkPereiraStd(lattice, solutionForwardPrecise, params) ;
% Save to file
if params.log.saveTempResults
    filename = [resultsFolder '/precise_' params.runId '_' num2str(iteration) '.mat'] ;
    saveResults(solutionForwardPrecise, lowerBound, meanCost, stds, iteration, lattice, params, filename) ;    
    displayMessage(sprintf('Precise results saved to %s', filename), params, 1) ;                
    displayMessage(sprintf('Precise computation done in %f s.',toc(startTimePrecise)), params, 1) ;     
end

end