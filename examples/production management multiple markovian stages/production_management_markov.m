% How to create a simple lattice
clc ; clear all ; close all ;

% The number of stages
H = 4 ;

% The number of product
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

% We star at the second node
firstNode = 2;

% Creating a simple H stages lattice with 2 nodes at second stage
% Demand is a random variable with value 2 or 10
lattice = Lattice.latticeEasyMarkov(H, transitionProba, firstNode) ;

% Visualisation
figure ;
lattice.plotLattice(@(data) ['d = ' num2str(data)]) ;

% Run SDDP
params = sddpSettings('algo.McCount',4,...
                      'stop.pereiraCoef',0.0000,...
                      'verbose',1,...
                      'algo.minTheta',-1e3,...
                      'solver','linprog') ;
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


