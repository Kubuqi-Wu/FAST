function [objOpt,diagnostics,solutionForward] = waitAndSee(lattice,path,params)
%  function WAITANDSEE
%   
%  wsValue = WAITANDSEE(lattice,path)   Computes the wait-and-see value 
%   (wsValue) of the precised path. The wait-and see value is given by 
%   optimising the problem with a perfect knowledge ofthe future.
%
%   WAITANDSEE(lattice,path,param)      Run the function WAITANDSEE with
%   some specified parameters. See SDDPSETTINGS for more information.
%
%   [wsValue,diag] = WAITANDSEE(...)    Returns also the diaggnostic of the
%   solver. See SOLVELP for more information.
%   
%   [...,solutionForward] = waitAndSee  Return also the cell array
%   solutionForward of size (H,1). This cell array contains the optimal
%   value of the variables of the problem.
%
% See also SDDP

H = lattice.H;
if(length(path) ~= H)
    error('The length of the path is not equal to field lattice.H.')
end

% Mem alloc
nnzMax = 0; % nnz: simply non-zero entries of W and T. (sum over H ofc)
nvar = 0;
nconstr = 0;
for i=1:H
    scenarioCurrent = lattice.graph{i}{path(i)} ;
    [nnzW,nnzT,nvarW,~,nconstrW,~] = scenarioCurrent.model.getSizeBlocks();
    nnzMax = nnzW + nnzT + nnzMax;
    nvar = nvar + nvarW;
    nconstr = nconstr + nconstrW;
end

% A = spalloc(nconstr,nvar,nnzMax);
% instead of A(i,j) = val, we store
%   Ai      = i
%   Aj      = j
%   Aval    = val
Ai = zeros(nnzMax,1);
Aj = zeros(nnzMax,1);
Aval = zeros(nnzMax,1);
b = zeros(nconstr,1);
c = zeros(nvar,1);
k = 0 ;

% Create WS problem :
% A = [
%   W1 0  0  0  ...
%   T2 W2 0  0  ...
%   0  T3 W3 0  ...
%   0  0  T4 W4 ...
%   ...              ];
%
% Where A is the constr matrix

scenarioCurrent = lattice.graph{1}{path(1)} ;

% First scenario : we don't have T
[c_ws,k_ws,W,h_ws,~] = scenarioCurrent.model.getWSblocks(0,[]);
[previousNContr,previousNVar] = size(W);

trialIdx = scenarioCurrent.model.forwardIdx;

% A(1:previousNContr,1:previousNVar) = W;
[Aitemp,Ajtemp,Avaltemp] = find(W);
Ai(1:length(Aitemp)) = Aitemp;
Aj(1:length(Aitemp)) = Ajtemp;
Aval(1:length(Aitemp)) = Avaltemp;

% Store b, c and k
b(1:previousNContr) = h_ws;
c(1:previousNVar) = c_ws;
k = k+k_ws ;


idxSwitchConstr = 0;
idxSwitchVar = 0;
counterSparseIdx = length(Aitemp);


for i=2:H
    scenarioCurrent = lattice.graph{i}{path(i)} ;   
    [c_ws,k_ws,W,h_ws,T] = scenarioCurrent.model.getWSblocks(previousNVar, trialIdx);    
    [NContr,NVar] = size(W);
    idxSwitchConstr = idxSwitchConstr + previousNContr;
    
    
    if(~isempty(T))
        % Add block T
        % A(idxSwitchConstr + (1:NContr),idxSwitchVar+(1:previousNVar)) = T;
        [Aitemp,Ajtemp,Avaltemp] = find(T);
        Ai(counterSparseIdx + (1:length(Aitemp))) = idxSwitchConstr + Aitemp;
        Aj(counterSparseIdx + (1:length(Aitemp))) = idxSwitchVar + Ajtemp;
        Aval(counterSparseIdx + (1:length(Aitemp))) = Avaltemp;
        counterSparseIdx = counterSparseIdx + length(Aitemp);
        idxSwitchVar = idxSwitchVar + previousNVar;
    end
    
    
    % Add block W
    % A(idxSwitchConstr + (1:NContr),idxSwitchVar+(1:NVar)) = W;
    [Aitemp,Ajtemp,Avaltemp] = find(W);
    Ai(counterSparseIdx + (1:length(Aitemp))) = idxSwitchConstr + Aitemp;
    Aj(counterSparseIdx + (1:length(Aitemp))) = idxSwitchVar + Ajtemp;
    Aval(counterSparseIdx + (1:length(Aitemp))) = Avaltemp;
    counterSparseIdx = counterSparseIdx + length(Aitemp);
    
    trialIdx = scenarioCurrent.model.forwardIdx;

    b(idxSwitchConstr + (1:NContr)) = h_ws;
    c(idxSwitchVar + (1:NVar)) = c_ws;
    k = k+k_ws ;
 
    [previousNContr,previousNVar] = size(W);
end

% Solve it
A = sparse(Ai,Aj,Aval,idxSwitchConstr+NContr,idxSwitchVar+NVar,nnzMax);
[x, ~, objOpt, diagnostics] = solveLP(A, b, c, params);
objOpt = objOpt + k ;

% Store solutionForward
switchIdxSol = 0;
solutionForward = cell(H,1);
for i=1:H
    solutionForward{i} = struct();
    scenarioCurrent = lattice.graph{i}{path(i)} ;
    nVar = size(scenarioCurrent.model.W,2);
    solutionForward{i}.x = x(switchIdxSol+(1:nVar));
    solutionForward{i}.xTrial = lattice.graph{i}{path(i)}.model.extractXTrial(solutionForward{i}.x);
    solutionForward{i}.costWithoutTheta = c(switchIdxSol+(1:nVar))'*solutionForward{i}.x;
    switchIdxSol = switchIdxSol + nVar;
end