% In this example we noticed that linprog wasn't robust enough
% It is advised to use another solver (e.g. gurobi) ; 
% you may still change the solver below (around line 45) and use linprog if
% desired
clc ; clear all ; close all ;

% The number of stages
H = 4 ;

% The number of scenarios
N = 3 ;

% Cost of production
C = rand(N, 1) ;

% Sale price
S = 2 + rand(N, 1) ;

% Creating the matrix of transition probabilities.
% At each stage, there is uniform(50,80)% chance of staying in the same stage
% or jump uniformly to another stage
randVector = 0.5 + 0.3*rand(N,1);
transitionProba = ones(N,N)/(N-1);
for i=1:N
    transitionProba(i,:) = transitionProba(i,:)*(1-randVector(i));
    transitionProba(i,i) = randVector(i);
end

% We start at the second node
firstNode = 2;

% Creating a simple H stages lattice with 2 nodes at second stage
% Demand is a random variable with value 2 or 10
lattice = Lattice.latticeEasyMarkov(H, transitionProba, firstNode) ;

% Visualisation
figure ;
lattice.plotLattice() ;

% Run SDDP
params = sddpSettings('algo.McCount',25,...
                      'stop.pereiraCoef',0.1,...
                      'stop.iterationMax',10,...
                      'verbose',1,...
                      'algo.minTheta',-1e3,...
                      'solver','gurobi') ; % Linprog is not robust enough for this example
var.x = sddpVar(N,H-1) ;
var.s = sddpVar(N,H-1) ;
lattice = lattice.compileLattice(@(scenario)production_management_markov_nlds(scenario,var,H,C,S),params) ; 


output = sddp(lattice, params) ;
plotOutput(output);

% Getting forward solution
lattice = output.lattice ;
nForward = 10 ;
for i = 1:nForward
    [~,~,~,solution] = forwardPass(lattice,lattice.randomPath(),params) ;  
    xVal = lattice.getPrimalSolution(var.x, solution) ;
end


