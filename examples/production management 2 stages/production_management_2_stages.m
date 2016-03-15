% How to create a simple lattice
clc ; clear all ; close all ;

% Creating a simple 2 stages lattice with 4 nodes at second stage
% Demand is a random variables, with scenarios 2, 3, 4 or 5
nNodes = 4;
H = 2;
demandVector = 1+(1:nNodes);
lattice = Lattice.latticeEasy(H, nNodes, @(t,i) demand(t,i,demandVector,nNodes)) ;
lattice = lattice.initExpectedLattice(@(t,i) demand(t,i,demandVector,nNodes)) ;


% Run SDDP
params = sddpSettings('stop.iterationMax',20,...
                      'algo.McCount',200,...
                      'stop.pereiraCoef',1e-6,...                     
                      'algo.deterministic',true,...
                      'verbose',1,...
                      'algo.minTheta',-20,...
                      'solver','gurobi') ;
x = sddpVar(2,1) ; % 1 product to produce
s = sddpVar(2,1) ; % and to sell
lattice = compileLattice(lattice,@(scenario)nlds(scenario,x,s),params);
lattice = compileExpectedLattice(lattice,@(scenario)nlds(scenario,x,s),params);

output = sddp(lattice,params) ;
lattice = output.lattice ;

plotOutput(output) ;
ylim([-10 10]) ;

% Getting solution
nForward = 5 ;
objVec = zeros(nForward,1);
xVec = zeros(nForward,2);
sVec = zeros(nForward,2);
for  i = 1:nForward
    [objVec(i),~,~,solution] = forwardPass(lattice,lattice.randomPath(),params) ;    
    xVec(i,:) = lattice.getPrimalSolution(x, solution) ;
    sVec(i,:) = lattice.getPrimalSolution(s, solution) ;
end

% Computing EV
[EV,~,solutionEV] = expectedValue(lattice,params);
EV

% Plotting value function
lattice.graph{1}{1}.plotCuts(1:2,[0 0],[5 5],false) ;
plotVx() ;




