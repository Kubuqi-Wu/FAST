function  [x, duals, objOpt, diagnostics] = solveLinprog(A, b, obj, opts)

% linprog solves
% min c' x     s.t.
%   Ax <= b
if(isempty(opts))
    opts = optimoptions(@linprog,'display','off','MaxIter',1e6) ;
end
opts = optimoptions(@linprog,opts,'algorithm','Dual-Simplex') ;
[x,objOpt,exitflag,output,lambda] = linprog(obj,-A,-full(b),[],[],[],[],[],opts) ;

if(exitflag<=0)
    warning('Linprog with dual-simplex failed. Trying interior-point.')
    exitflag
    opts = optimoptions(@linprog,opts,'algorithm','interior-point') ;
    [x,objOpt,exitflag,output,lambda] = linprog(obj,-A,-full(b),[],[],[],[],[],opts) ;
end

if(exitflag<=0)
    warning('Linprog with dual-simplex failed. Trying legacy interior-points.')
    exitflag
    opts = optimoptions(@linprog,opts,'algorithm','interior-point-legacy') ;
    [x,objOpt,exitflag,output,lambda] = linprog(obj,-A,-full(b),[],[],[],[],[],opts) ;
end


duals = lambda.ineqlin ;
if exitflag == 1
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('fast:solveLinprog:notSolved','Error while solving problem using linprog. Exitflag %d. Try using another solver (like cplex or gurobi) or changing the solver options.',exitflag) ;
end
diagnostics.explanation = exitflag ;
diagnostics.objective = objOpt ;

