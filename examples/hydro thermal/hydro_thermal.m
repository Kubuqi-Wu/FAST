% HYDRO THERMAL example where uncertainty is the amount of rainfall at each
% stage. Rainfall can either be high (10) or small (2).
% We can use some fuel at a price C (5) to meet demand d (6).
% From stage to stage, water can been stored in the reservoir, but there is
% a tank limit V (8).

% With these values, in the case of a two stage problem, the value function
% V(x1) at stage 1 is
%   V(x1) = 30 - 0.5 ( 5 min(x1+2,6) + 5 min(x1+10,6) )
%         = 15 - 5/2 min(x1+2,6)
% The problem at stage 1 is then
%   min 5 p1 + V(x1) s.t. x1 <= 8, x1 <= 6 - y1, p1 + y1 >= 6
% The solution at stage 1 can then trivially be found to be
%   x1 = 0, y1 = 6, p1 = 0
% for a cost of 10.
% This is (approximately) the value you will find by running this example
% If you want the exact solution, use the option 
%   'algo.deterministic',true
% This will iterate over all possible samples in order to build an exact
% meanCost, hence providing the exact solution.

% How to create a simple lattice
clc ; close all ; clear all ;
rng(100); % To avoid some issues with linprog

H = 5 ;

% Creating a simple 5 stages lattice with 2 nodes at second stage
lattice = Lattice.latticeEasy(H, 2, @rainfall) ;

% Visualisation
figure ;
lattice.plotLattice(@(data) num2str(data)) ;

% Run SDDP
params = sddpSettings('algo.McCount',25, ...
                      'stop.iterationMax',10,...                      
                      'stop.pereiraCoef',2,...                   
                      'solver','linprog') ;
var.x = sddpVar(H) ; % The reservoir level at time t
var.y = sddpVar(H) ; % For how much we use the water at time t
var.p = sddpVar(H) ; % For how much we use the fuel generator at time t                  
lattice = compileLattice(lattice,@(scenario)hydro_thermal_nlds(scenario,var),params) ;                                    
output = sddp(lattice,params) ;

% Visualise output
plotOutput(output) ;

% Forward passes
lattice = output.lattice ;
nForward = 5 ;
objVec = zeros(nForward,1);
x = zeros(nForward,H);
y = zeros(nForward,H);
p = zeros(nForward,H);
dataForward = cell(nForward,1);
for  i = 1:nForward
    [objVec(i),~,~,solution] = forwardPass(lattice,'random',params) ; 
    x(i,:) = lattice.getPrimalSolution(var.x, solution) ;
    y(i,:) = lattice.getPrimalSolution(var.y, solution) ;
    p(i,:) = lattice.getPrimalSolution(var.p, solution) ;
    dataForward{i} = lattice.getDataSolution(solution);
end
