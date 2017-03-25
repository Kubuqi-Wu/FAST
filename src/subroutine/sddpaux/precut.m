function lattice = precut(lattice, variables, values, params)
%PRECUT Precut the lattice
%
% lattice = precut(lattice, variables, values) use the solution stored in
%   the (variables, values) pair to "precut" the lattice. More precisely,
%   it uses such solution (note that it should be feasible at each stage)
%   to run a backward pass on the lattice, and create cuts.
%   variables is an array of sddpVar, and values are the corresponding
%   values
%   All the variables of the model should be present in variables, in any
%   order. Additional variables can be present. In this case, their
%   corresponding value is ignored.
%   For variables from stage t that do not appear in stage t+1, the value
%   can be arbitrary.
%
% See also SDDP BACKWARDPASS

if any(size(variables) ~= size(values))
    error('fast::precut variables and values should have the same size')
end

% We need to build a "solutionForward" like structure. Basically, we need
% the trials
varAllStages = [variables.ids]' ;
H = lattice.getH() ;
solutionForward = cell(1, 1) ;
solutionForward{1}.solutionForwardCells = cell(H, 1);
for t = 1:H-1
    % Get the x at stage i
    model = lattice.graph{t}{1}.model ; % Structure is uniform accross nodes at a given stage
    varStageI = model.modelVarIdx ;     % Indexes of the variables in modeling space
    valStageI = zeros(size(varStageI)) ;
    for i = 1:numel(valStageI)
        idxi = find(varAllStages == varStageI(i)) ;
        if isempty(idxi)
            error('fast::precut a variable is missing in the variables array') ;
        end
        valStageI(i) = values(idxi) ;
    end   
    solutionForward{1}.solutionForwardCells{t}.trials = model.extractXTrial(valStageI) ;
end

% Run the backwardPass !
lattice = backwardPass(lattice, solutionForward, params) ;

