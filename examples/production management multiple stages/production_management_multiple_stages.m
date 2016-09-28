% How to create a simple lattice
clc ; clear all ; close all ;
rng(100); % To avoid some issues with linprog

% The number of stages
H = 3 ;

% The number of product
N = 2 ;

% Cost of production
C = rand(N, 1) ;

% Sale price
S = 2 + rand(N, 1) ;

% Creating a simple H stages lattice with 2 nodes at second and more stage
% Demand is a random variable with value 2 or 10
lattice = Lattice.latticeEasy(H, 2) ;

% Visualisation
figure ;
lattice.plotLattice(@(data) ['d = ' num2str(data)]) ;

% Run SDDP
params = sddpSettings('algo.McCount',25,...
                      'stop.pereiraCoef',0.1,...
                      'stop.iterationMax',10,...
                      'verbose',1,...
                      'algo.minTheta',-1e3,...
                      'solver','linprog') ;
var.x = sddpVar(N,H-1) ;
var.s = sddpVar(N,H-1) ;
lattice = lattice.compileLattice(@(scenario)production_management_multiple_stages_nlds(scenario,var,H,C,S),params) ; 


output = sddp(lattice, params) ;
plotOutput(output);

% Getting forward solution
lattice = output.lattice ;
nForward = 10 ;
for i = 1:nForward
    [~,~,~,solution] = forwardPass(lattice,lattice.randomPath(),params) ;  
    xVal = lattice.getPrimalSolution(var.x, solution) ;
end


