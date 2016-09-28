function  [x, duals, objOpt, diagnostics] = solveLinprog(A, b, obj, opts)

% linprog solves
% min c' x     s.t.
%   Ax <= b
if(isempty(opts))
    opts = optimset('display','off','MaxIter',1e6) ;
end
[x,objOpt,exitflag,output,lambda] = linprog(obj,-A,-full(b),[],[],[],[],[],opts) ;
duals = lambda.ineqlin ;
if exitflag == 1
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('fast:solveLinprog:notSolved','Error while solving problem using linprog. Exitflag %d. Try using another solver (like cplex or gurobi) or changing the solver options.',exitflag) ;
end
diagnostics.explanation = exitflag ;
diagnostics.objective = objOpt ;

