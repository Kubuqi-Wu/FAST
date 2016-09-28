% HYDRO THERMAL example where uncertainty is the amount of rainfall at each
% stage. 

% How to create a simple lattice
clc ; close all ; clear all ;

H = 6 ;

% Creating a simple 6 stages lattice with 2 nodes at second stage
lattice = Lattice.latticeEasy(H, 2, @rainfallUnif) ;

% Visualisation
figure ;
lattice.plotLattice(@(data) num2str(data)) ;

% Run SDDP
params = sddpSettings('algo.McCount',5, ...
                      'stop.iterationMax',10,...                      
                      'stop.pereiraCoef',0.1,...
                      'solver','linprog') ;
var.x = sddpVar(H) ;
var.p = sddpVar(H) ;
var.r = sddpVar(H) ;
var.q = sddpVar(H) ;
var.l = sddpVar(H) ;
lattice = compileLattice(lattice,@(scenario)hydro_thermal_ar_nlds(scenario,var),params) ;                                    
output = sddp(lattice,params) ;

% Visualise output
plotOutput(output) ;