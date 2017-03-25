% How to create a simple lattice
clc ; clear all ; close all ;

% Creating a simple 2 stages lattice with 4 nodes at second stage
% Demand is a random variables, with scenarios 2, 3, 4 or 5
nNodes = 4;
H = 2;
demandVector = 1+(1:nNodes);
lattice = Lattice.latticeEasy(H, nNodes, @(t,i) production_management_2_stages_demand(t,i,demandVector,nNodes)) ;
lattice = lattice.initExpectedLattice(@(t,i) production_management_2_stages_demand(t,i,demandVector,nNodes)) ;


% Run SDDP
params = sddpSettings('stop.iterationMax',20,...
                      'algo.McCount',10,...
                      'stop.pereiraCoef',1e-1,...                     
                      'algo.deterministic',false,...
                      'verbose',1,...
                      'algo.minTheta',-20,...
                      'solver','linprog') ;
x = sddpVar(2,1) ; % 1 product to produce
s = sddpVar(2,1) ; % and to sell
lattice = compileLattice(lattice,@(scenario)production_management_2_stages_nlds(scenario,x,s),params);
lattice = compileExpectedLattice(lattice,@(scenario)production_management_2_stages_nlds(scenario,x,s),params);

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
fprintf('EV = %f\n', EV) ;

% Getting cuts
[cutCoeffs, cutRHS] = lattice.getCuts(1, 1, x) ;
% Plotting cuts + objective
% since coeffs * x + theta >= rhs
% we have theta >= rhs - coeffs^t x
% and we can thus plot that
[x1, x2] = meshgrid(0:0.1:5, 0:0.1:5) ;
vx = zeros(size(x1)) ;
obj = zeros(size(x1)) ;
for i = 1:size(x1, 1)
    for j = 1:size(x1, 2) 
        vxall = cutRHS - cutCoeffs * [x1(i,j) ; x2(i,j)] ;        
        vx(i,j) = max(vxall) ; % Since we take the max
        obj(i,j) = x1(i,j) + 2*x2(i,j) ;
    end
end
figure ; 
surf(x1, x2, vx + obj) ;

% Plotting the real analytical value
plotVx() ;




