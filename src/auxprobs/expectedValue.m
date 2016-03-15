function [EV,diagnostics,solutionExpectedValue] = expectedValue(lattice,params)

H = lattice.H;

% Prepare memory (can be big)
nnzMax = 0;
nvar = 0;
nconstr = 0;
for i=1:H
    scenarioCurrent = lattice.expectedGraph{i} ;
    [nnzW,nnzT,nvarW,~,nconstrW,~] = scenarioCurrent.model.getSizeBlocks();
    nnzMax = nnzW + nnzT + nnzMax;
    nvar = nvar + nvarW;
    nconstr = nconstr + nconstrW;
end
% A = spalloc(nconstr,nvar,nnzMax);
Ai = zeros(nnzMax,1);
Aj = zeros(nnzMax,1);
Aval = zeros(nnzMax,1);
b = zeros(nconstr,1);
c = zeros(nvar,1);
k = 0 ;

scenarioCurrent = lattice.expectedGraph{1}  ;

% First scenario : we don't have T
[c_ws,k_ws,W,h_ws,~] = scenarioCurrent.model.getWSblocks(0,[]);
[previousNContr,previousNVar] = size(W);

k = k+k_ws ;

trialIdx = scenarioCurrent.model.forwardIdx;

% A(1:previousNContr,1:previousNVar) = W;
[Aitemp,Ajtemp,Avaltemp] = find(W);
Ai(1:length(Aitemp)) = Aitemp;
Aj(1:length(Aitemp)) = Ajtemp;
Aval(1:length(Aitemp)) = Avaltemp;
b(1:previousNContr) = h_ws;
c(1:previousNVar) = c_ws;

idxSwitchConstr = 0;
idxSwitchVar = 0;
counterSparseIdx = length(Aitemp);

% Build the big big matrix
for i=2:H
    scenarioCurrent = lattice.expectedGraph{i} ;   
    [c_ws,k_ws,W,h_ws,T] = scenarioCurrent.model.getWSblocks(previousNVar, trialIdx);
    k = k+k_ws ;
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
 
    [previousNContr,previousNVar] = size(W);
end

% Solve
A = sparse(Ai,Aj,Aval,idxSwitchConstr+NContr,idxSwitchVar+NVar,nnzMax);
[x, ~, EV, diagnostics] = solveLP(A, b, c, params);
EV = EV+k ;

% Retreive solution
switchIdxSol = 0;
solutionExpectedValue = cell(H,1);
for i=1:H
    solutionExpectedValue{i} = struct();
    scenarioCurrent = lattice.expectedGraph{i} ;
    nVar = size(scenarioCurrent.model.W,2);
    solutionExpectedValue{i}.primal = x(switchIdxSol+(1:nVar));
    solutionExpectedValue{i}.trials = lattice.expectedGraph{i}.model.extractXTrial(solutionExpectedValue{i}.primal);
    solutionExpectedValue{i}.costWithoutTheta = c(switchIdxSol+(1:nVar))'*solutionExpectedValue{i}.primal;
    switchIdxSol = switchIdxSol + nVar;
end