%% Example of modeling
% Shows what already works
clc ; close all ; clear all ;

A = [1 2 3 4 ; 7 8 9 10] ;
b = [1 ; 4] ;
x = sddpVar(4,1) ;
% min |Ax - b|_1 <=> min sum x_i, x >= 0, Ax == b

cntr1 = x >= 0 ;
cntr2 = A*x == b ;
cntr = [cntr1 ; cntr2] ;
obj = sum(x) ;


milpModel = MILPModel(cntr, obj) ;

%% Let's now call Gurobi
gurobiModel.obj = milpModel.c ;
gurobiModel.A = milpModel.A ;
gurobiModel.rhs = milpModel.b ;
gurobiModel.lb = - inf * ones(size(milpModel.c)) ;
gurobiModel.sense = '<' ;
vtypes(1:size(A,2)) = 'C' ;
vtypes(milpModel.binaryIdx) = 'C' ;% 'B' ; % Otherwise no dual variables
vtypes(milpModel.integerIdx) = 'C' ;% 'I' ;
gurobiModel.vtype = char(vtypes) ;
gurobiParams.outputflag = 0 ;
result = gurobi(gurobiModel, gurobiParams) ;
if strcmpi(result.status,'OPTIMAL')
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('Error while solving problem using gurobi') ;
    disp(result) ;
end

% We want the results :-)
xVal = milpModel.primalValue(x, result.x) ;

% We also want dual variables :-)
cntr1Dual = milpModel.dualValue(cntr1, result.pi) ;
cntr2Dual = milpModel.dualValue(cntr2, result.pi) ;

