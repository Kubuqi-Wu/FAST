% DEMO OPTIMIZE STOCK
%
% Assume we want to solve the following problem :
% We work on a two-step timeline.
% At time 1, we can decide to produce some amount of a product.
% At time 2, we can sell this product, up to some demand.
% But the problem is that demand at time 2 is unknown and random.
%
% Let's define by C the cost of production, P the sale price.
% Let's denote by x the amount of product we produce at time 1, and by s
% the amount of product we sell at time 2.
% Finally, let d be the (random) demand at time 2.
%
% Formally, the problem at time 1 can be stated as
%           min C x + E(V(x))
%           s.t.  x >= 0
% where E(V(x)) is the expectation of the cost of future stages.
% At stage 2, we thus have 
%           V(x) = min - P s
%                  s.t.  s >= 0
%                        s <= d(xi)
%                        s <= x
% where d(xi) is a random variable.
% Now let's see how to solve this problem using FAST.

clc ; close all ; clear all ;

%% 1. Lattice creation
% First, let's build the lattice.
% This lattice is defined over H = 2 stages
H = 2 ;
% using nNodes = 2 nodes at each stage, except a the first one where we
% often have only one node.
nNodes = 2 ;
lattice = Lattice.latticeEasy(H, nNodes, @demand) ;
% demand(t,i) is the function that will be used to store the random demand
% at each node. t is the time (1 or 2) and i the index (1 or 2)
% The random demand will be either 2 or 3, with equal probabilities.
% For that, de fine the demand function. See demand.m

% Let's see what the lattice looks like.
% The following 2 lines will display the lattice, as well as transition
% probabilities (1/2 by default) and the data stored, using the @(data) num2str(data)
% to convert the output of the demand function into string. This is
% optional, and the following two lines will work too
% figure ;
% lattice.plotLattice() ;
figure ;
lattice.plotLattice(@(data) num2str(data)) ;

%% 2. Variables declaration
% Now, we need to defined the variables we will use in our problem.
% Here, we only need x (of size 1 x 1) and y (1 x 1).
% Se we define then like that. Note that the 1,1 are optional, but it is to
% emphasize that x is a scalar, and so is y.
x = sddpVar(1,1) ;
s = sddpVar(1,1) ;

%% 3. NDLS declaration
% Finally, we need to create the function that, at each node, return the
% NDLS problem
% This is defined in the nlds function. See nlds.m.
lattice = lattice.compileLattice(@(scenario)nlds(scenario,x,s)) ;    

%% 4. Run SDDP !
% The followings are the settings.
% You can leave the default one, but it's a good idea to take a look at the
% documentation to see what all of them have to do
params = sddpSettings('algo.McCount',25, ...
                      'stop.iterationMax',10,...
                      'stop.pereiraCoef',2,...
                      'solver','gurobi') ; % You need to adapt this for your solver.                  
output = sddp(lattice,params) ;
lattice = output.lattice ;

% Visualise the evolution of the mean cost and the lower bound during the
% algorithm.
% They should meet at the end
plotOutput(output) ;
% We also see that the cost at the end is -2.

%% 5. Finally, retreive the solution.
% To do so, we do some 'forwardPass', and for each forward pass we extract
% the solution x and s.
nForward = 10 ;
for i = 1:nForward
    [~,~,~,solution] = forwardPass(lattice,'random',params) ;  
    xVal(i) = lattice.getPrimalSolution(x, solution) ;
    sVal(i) = lattice.getPrimalSolution(s, solution) ;
end
disp('x') ;
disp(xVal) ;
disp('s') ;
disp(sVal) ;