% HYDRO THERMAL example where uncertainty is the amount of rainfall at each
% stage. Rainfall can either be high (10) or small (2).
% We can use some fuel at a price C (5) to meet demand d (6).
% From stage to stage, water can been stored in the reservoir, but there is
% a tank limit V (8).
%
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

H = 2 ;

% Creating a simple 6 stages lattice with 2 nodes at second stage
lattice = Lattice.latticeEasy(H, 2, @rainfall) ;

% Visualisation
figure ;
lattice.plotLattice(@(data) num2str(data)) ;

% Run SDDP
params = sddpSettings('algo.McCount',2, ...
                      'stop.iterationMax',10,...                      
                      'stop.pereiraCoef',0.01,...
                      'algo.deterministic',true,...
                      'solver','gurobi') ;
var.x = sddpVar(H) ; % The reservoir level at time t
var.y = sddpVar(H) ; % For how much we use the water at time t
var.p = sddpVar(H) ; % For how much we use the fuel generator at time t                  
lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;                                    
output = sddp(lattice,params) ;

% Visualise output
plotOutput(output) ;

% Visualise cuts
plotCuts(output.lattice.graph{1}{1},1,0,10,false) ;