function [objectiveVal,costByStage,lbByStage,solutionForward] = forwardPass(lattice,path,params)
% FORWARDPASS Run a forward pass on the lattice
%
% [objectiveVal,costByStage,lbByStage,solutionForward] = FORWARDPASS(lattice,path,params)
%   runs one forward pass on the lattice lattice. 
%   - lattice needs to be compiled (i.e., to contain the models) and, hopefully, some cuts (though this is
%     not formally required).
%   - path can either be 'random' or a sequence of integers such that path(i)
%     is the desired node at stage i.
%   - params is the output of sddpSettings.
%   
%   The function returns
%   - objectiveVal, the objective value of the given forward pass
%   - costByStage, the details for each stage. objectiveVal =
%     sum(costByStage)
%   - lbByStage is the lower-bound at each stage
%   - solutionForward is a cell-array of size Hx1 where solutionForward{t}
%     is a struct() with the solution of stage t (with the same values as
%     in sddp)
%
%   See also SDDP, LATTICE, SDDPSETTINGS, BACKWARDPASS
 
H = lattice.getH();
costByStage = zeros(H,1);
lbByStage = zeros(H,1);
solutionForward = cell(H,1);
scenarioCurrent = [] ;
for time = 1:H   
    if isnumeric(path)
        scenarioCurrent = lattice.graph{time}{path(time)} ;              
    elseif strcmp(path,'random')
        scenarioCurrent = lattice.nextRandomScenario(scenarioCurrent) ;
    else    
        error('fast::forwardPass::path should either be a vector of index or ''random''.') ;
    end
    if(time > 1)
        solutionForward{time} = scenarioCurrent.solve(solutionForward{time-1}, time==H, params) ;
    else
        solutionForward{1}    = scenarioCurrent.solve([], time==H, params) ;
    end
    costByStage(time) = solutionForward{time}.costWithoutTheta ;
    lbByStage(time) = solutionForward{time}.costWithTheta ;
end
objectiveVal = sum(costByStage);