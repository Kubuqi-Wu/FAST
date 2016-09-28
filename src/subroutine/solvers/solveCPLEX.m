function  [x, duals, objOpt, diagnostics] = solveCPLEX(A, b, obj, opts)

% CPLEX solves
% min c' x     s.t.
%   Ax <= b
if(isempty(opts))
    opts = cplexoptimset('display','off') ;
end
[x,objOpt,exitflag,output,lambda] = cplexlp(obj,-A,-full(b),[],[],[],[],[],opts) ;
duals = lambda.ineqlin ;
if exitflag == 1
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('fast:solveCplex:notSolved','Error while solving problem using CPLEX. Exitflag %d',exitflag) ;
end
diagnostics.explanation = exitflag ;
diagnostics.objective = objOpt ;

