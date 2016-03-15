%% Example of modeling
% Shows what already works
clc ; close all ; clear all ;

%% Let's now call Gurobi
gurobiModel.obj = -5 ;
gurobiModel.A = sparse([2 ; -2]) ;
gurobiModel.rhs = [3 ; -3] ;
gurobiModel.lb = - inf ;
gurobiModel.sense = '<' ;
vtypes(1) = 'C' ;
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

result.x

result.pi