function [converged, lowerBound, meanCost, stdVal] = checkPereiraStd(lattice, solutionForward, params)

stopWhen = params.stop.stopWhen ; % 'pereira', 'std', 'pereira and std', 'never'
pereiraCoef = params.stop.pereiraCoef ;
stdMcCoef = params.stop.stdMcCoef ; 

McCount = size(solutionForward, 1) ;
H       = size(solutionForward{1}.solutionForwardCells, 1) ;

% LowerBound
costFirstStage = zeros(McCount, 1) ;
for m = 1:McCount
    costFirstStage(m) = solutionForward{m,1}.solutionForwardCells{1,1}.costWithTheta ; % If we have only one node at first stage, they should all be the same
end
if params.algo.deterministic
    lowerBound = dot(lattice.scenarioTableProba,costFirstStage) ;
else
    lowerBound = mean(costFirstStage) ; % They should all be the same
end

% UpperBound & Std
cost = zeros(H, McCount) ;
for t = 1:H
    for m = 1:McCount
        cost(t,m) = solutionForward{m,1}.solutionForwardCells{t,1}.costWithoutTheta ;
    end
end
cost = sum(cost,1) ;
stdVal = std(cost) ; 
if params.algo.deterministic
    meanCost = dot(lattice.scenarioTableProba,cost) ;
else
    meanCost = mean(cost) ;
end

% Check criterions
peireiraLb = meanCost - pereiraCoef * stdVal / sqrt(McCount) ;
peireiraUb = meanCost + pereiraCoef * stdVal / sqrt(McCount) ;
pereiraCriterionMet = (peireiraLb <= lowerBound) && (lowerBound <= peireiraUb) ;
stdMcCriterionMet = stdVal / sqrt(McCount) <= stdMcCoef * abs(lowerBound) ;

% Return wether, according to stopWhen, we should stop or not
checkPereira = strcmp(stopWhen,'pereira') || strcmp(stopWhen,'pereira and std') ;
checkStdMc   = strcmp(stopWhen,'std')     || strcmp(stopWhen,'pereira and std') ;
if strcmp(stopWhen,'pereira')
    converged = pereiraCriterionMet ;
elseif strcmp(stopWhen,'std')
    converged = stdMcCriterionMet ;
elseif strcmp(stopWhen,'pereira and std')
    converged = pereiraCriterionMet && stdMcCriterionMet ;
elseif strcmp(stopWhen,'never')
    converged = false ;
else
    error('Wrong params.stop.stopWhen parameter') ;
end

% Display stuff
if params.verbose >= 1
    fprintf('LowerBound                                 : %e\n', lowerBound) ;
    fprintf('Mean(ForwardCosts)   (K = %3.1d)             : %e\n', McCount, meanCost) ; 
    fprintf('Std(ForwardCosts)    (K = %3.1d)             : %e\n', McCount, stdVal) ;
    fprintf('95 pc confidence interval around mean cost : [%e   %e]\n', meanCost-2*stdVal/sqrt(McCount), meanCost+2*stdVal/sqrt(McCount)) ;
    fprintf('95 pc confidence interval for solution     : [%e   %e]\n', lowerBound, meanCost+2*stdVal/sqrt(McCount)) ;
    fprintf('Confidence interval desired (coef %2.1e) : [%e   %e]\n', pereiraCoef, peireiraLb, peireiraUb) ;
    fprintf('Pereira''s  criterion (%17s)   : %s\n',checkOrNot(checkPereira), metOrNot(pereiraCriterionMet)) ;
    fprintf('StdMc      criterion (%17s)   : %s\n',checkOrNot(checkStdMc),    metOrNot(stdMcCriterionMet)) ;   
end

end

function str = metOrNot(boolean)
if boolean
    str = 'met' ;
else
    str = 'not met' ;
end
end

function str = checkOrNot(boolean)
if boolean
    str = 'to be checked' ;
else
    str = 'not to be checked' ;
end
end