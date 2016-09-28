function  [x, duals, objOpt, diagnostics] = solveGurobi(A, b, obj, gurobiParams)

% gurobi solves
% min c' x      s.t.
%     Ax >= b

gurobiModel.obj = obj ;
gurobiModel.A = sparse(A) ;
gurobiModel.rhs = full(b) ;
gurobiModel.lb = - inf * ones(size(obj)) ;
gurobiModel.sense = '>' ;
if(isempty(opts))
    gurobiParams.outputflag = 0 ;
end

result = gurobi(gurobiModel, gurobiParams) ;

diagnostics = [] ;

if strcmpi(result.status,'OPTIMAL')
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('fast:solveGurobi:notSolved','Error while solving problem using gurobi.\nThe solver returned:') ;
    disp(result) ;
end

x = result.x ;
duals = result.pi ;

diagnostics.explanation = result ;
objOpt = result.objval ;
diagnostics.objective = objOpt ;