function [objectiveVal,costByStage,solutionForward] = forwardPassFixedDecision(lattice,path,decision,params)
% We suppose that the lattice is already compiled
% decision (cell{double[]}) are values of variables. We will fix only the fowardDep
% and backwardDep variables (i.e. decisions variables).

warning('To recode')

H = lattice.getH();
solutionForward = cell(H,1);
costByStage = zeros(1,H);
for time = 1:H      
    scenarioCurrent = lattice.graph{time}{path(time)} ;              
    if(time > 1)
        [solutionForward{time}, constraints] = scenarioCurrent.solveUnderDecision(decision{time-1}, time==H, decision{time}, params) ;
    else
        [solutionForward{1},    constraints] = scenarioCurrent.solveUnderDecision([], time==H, decision{time}, params) ;
    end            
    solutionForward{time}.dualsMatrix = constraints ;
    scenarioCurrent.assignValueToVariables(solutionForward{time}.x) ;
    costByStage(time) = solutionForward{time}.costWithoutTheta;
end
objectiveVal = sum(costByStage);
fixYalmipVariables(lattice,solutionForward)