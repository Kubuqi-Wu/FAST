% How to create a simple lattice
clc ; close all ; clear all ;

% Creating a simple 2 stages lattice with 2 nodes at second stage
lattice = Lattice.latticeEasy(2, 2, @optimize_stock_demand) ; % 2 stages with 2 nodes at each stage
lattice = lattice.initExpectedLattice(@optimize_stock_demand);

% Visualisation
figure ;
lattice.plotLattice(@(data) num2str(data)) ;

% Run SDDP
params = sddpSettings('algo.McCount',2, ...
                      'stop.iterationMax',10,...                      
                      'precise.computeEnd',true,...
                      'precise.count',20,...
                      'solver','linprog') ;
var.x = sddpVar(1,1) ;
var.s = sddpVar(1,1) ;                 
lattice = compileLattice(lattice,@(scenario)optimize_stock_nlds(scenario,var),params) ;     
lattice = compileExpectedLattice(lattice,@(scenario)optimize_stock_nlds(scenario,var),params) ;                  
output = sddp(lattice,params) ;

% Expected and obtained V(x) and cuts
output.lattice.graph{1}{1}.plotCuts(1,0,10,false) ;
x = linspace(0, 5, 100) ;
z = - min(x, 2) - min(x, 3) ;
plot(x, z, '-g') ;

% Solution should be - 2.5
lattice = output.lattice;
wsVal = 0;
nWS = 100 ;
for i = 1:nWS
    ws = waitAndSee(lattice,lattice.randomPath(),params);
    wsVal = wsVal + ws/nWS;
end
display(wsVal)

% ExpectedValue
[EV_val,~,EV_sol] = expectedValue(lattice,params);


