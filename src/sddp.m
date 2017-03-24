function output = sddp(lattice,params)  
%% SDDP Stochastic Dynamic Dual Programming
%
% output = SDDP(lattice)
% output = SDDP(lattice,params) solves the problem depicted by lattice
% params is optionnal, and should be the result of sddpSettings().
%
% Returns output, a struct containing:
% * output.lattice, the lattice updated with the cuts
% * output.lowerBounds, a vector of size [iteration x 1] containing all the lowerBounds
% * output.meanCost, a vector of size [iteration x 1] containing all the mean cost
% * output.stds, a vector of size [iteration x 1] containg the std of the
%   mean cost at each iteration
% * output.runningTime, the running time of the algorithm in seconds
% * output.params, the parameter given to sddp. If no parameter have been
%   provided, output.params = sddpSettings().
% * output.solution is a cell array of size [H x McCount]. Each entry is a
%   structure containing the result of each Nlds :
%     - solution.primal, the primal solution
%     - solution.trials, the trials sent to the next node
%     - solution.dualCntr, the dual of the constraints
%     - solution.dualCuts, the dual of the cuts
%     - solution.costWithoutTheta, the cost
%     - solution.costWithTheta, the cost + theta
%
% See also SDDPSETTINGS, LATTICE, WAITANDSEE

%% Pre-Sddp Phase
if nargin == 1
    params = sddpSettings() ;
end

% Logging
if params.log.useDiary
    diary(params.log.logFile) ;
end
displayMessage(sprintf('FAST (version %s-%s) is now running.',params.version.id,params.version.type),params, 1) ;
displayMessage(params.version.warning,params,1) ;
displayMessage(params.version.copyright,params,1) ;

% Displaying params
displayParams(params) ;

% Creating results folder if it doesn't exist
if params.log.saveTempResults    
    resultsFolder = GetFullPath(params.log.resultsFolder) ;
    resultsFolder = strrep(resultsFolder, '\', '/');
    if ~ exist(resultsFolder,'dir')
        mkdir(resultsFolder) ;
        displayMessage(sprintf('Results folder created at %s.', resultsFolder), params, 2) ;
    end
else
    resultsFolder = [] ;
end

% Handling NLDS : put McCount to nScenarios
if params.algo.deterministic
    displayMessage('Deterministics SDDP (= NLDS) has been choosen. McCount is discarded.', params, 1) ;
    params.algo.McCount = lattice.getNScenarios() ;
    displayMessage(sprintf('Deterministic SDDP will do %d "simulation" at each pass.',params.algo.McCount),params, 1) ;
    lattice = lattice.generateScenarioTable() ;
end

%% Running SDDP
converged = false ;
iterationMax = params.stop.iterationMax ;
H = lattice.getH() ;
McCount = params.algo.McCount ;
startTimeAlgo = tic() ;

% Display generated models
if params.debug
    lattice.displayModels() ;
end

% Diagnostic
iteration      = 1 ;
lowerBounds    = zeros(iterationMax, 1) ;
meanCost       = zeros(iterationMax, 1) ;
stds           = zeros(iterationMax, 1) ;

while ~ converged
    
    % Iteration info
    displayMessage(sprintf('-------------\nIteration %d\n-------------', iteration), params, 1) ;
    startTime = tic() ;
    solutionForward = cell(McCount,1);
    
    % Forward pass
    for Mc = 1:McCount
        displayMessage(sprintf('%d) forward pass',Mc), params, 2) ;        
        if params.algo.deterministic 
            path = lattice.deterministicPath(Mc) ;
        else
            path = lattice.randomPath() ;
        end              
        [~,~,~,solutionForward{Mc,1}] = forwardPass(lattice, path, params) ;
    end
    
    % Backward pass
    lattice = backwardPass(lattice, solutionForward, params) ;    
    
    % Convergence and stuff
    [convergedRegular, lowerBounds(iteration), meanCost(iteration), stds(iteration)] = checkPereiraStd(lattice, solutionForward, params) ;
    if params.precise.compute % If we want to compute precise UB
        convergedPrecise = computePreciseUB(lattice, resultsFolder, iteration, params) ;
    else
        convergedPrecise = false ;
    end
    
    % Check stopping criterion
    if iteration < params.stop.iterationMin || toc(startTimeAlgo) < params.stop.timeMin
        converged = false ;
    elseif iteration < params.stop.iterationMax && toc(startTimeAlgo) < params.stop.timeMax
        if params.stop.regular && convergedRegular
            converged = true ;
            displayMessage('Convergence in the regular setting.', params, 1) ;
        elseif params.stop.precise && convergedPrecise
            converged = true ;
            displayMessage('Convergence in the precise setting.', params, 1) ;
        else
            converged = false ;
        end
    elseif params.stop.iterationMax <= iteration
        displayMessage(sprintf('Maximun iteration reached (%d/%d). Aborting.',iteration,params.stop.iterationMax), params, 1) ;
        converged = true ;
    elseif params.stop.timeMax <= toc(startTimeAlgo)
        displayMessage(sprintf('Maximun time allowed reached (%d/%d). Aborting.',toc(startTimeAlgo),params.stop.timeMax), params, 1) ;
        converged = true ;
    end
    
    % Purge cuts
    if params.algo.purgeCuts
        lattice = lattice.purgeCuts(params);
    end
    
    % Save to file
    if params.log.saveTempResults
        filename = [resultsFolder '/temp_' params.runId '_' num2str(iteration) '.mat'] ;
        saveResults(solutionForward, lowerBounds(1:iteration), meanCost(1:iteration), stds(1:iteration), iteration, lattice, params, filename) ;
        displayMessage(sprintf('Temp results saved to %s.', filename), params, 1) ;
    end
    
    % Display infos
    iterationDuration = toc(startTime) ;
    displayMessage(sprintf('This iteration took %f s.', iterationDuration), params, 1) ;
    
    % Iteration ++
    iteration = iteration + 1 ;
end

% Compute Precise UB one last time, as well as final solution
if params.precise.computeEnd
    displayMessage('Computing precise solution one last time.', params, 1) ;
    [~, solutionFinal] = computePreciseUB(lattice, resultsFolder, iteration, params) ;
    output.solution    = solutionFinal ;
else
    output.solution    = solutionForward ;
end

% Output stuff
output.lattice = lattice ;
output.lowerBounds = lowerBounds(1:iteration-1) ;
output.meanCost    = meanCost(1:iteration-1) ;
output.stds        = stds(1:iteration-1) ;
runningTime = toc(startTimeAlgo) ;
output.runningTime = runningTime ;
output.params = params ;

displayMessage(sprintf('The whole algorithm took %f s.', runningTime), params, 1) ;
displayMessage('FAST stopped.', params, 1) ;

if params.log.useDiary
    diary off ;
end

end