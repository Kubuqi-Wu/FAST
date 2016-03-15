function  [x, duals, objOpt, diagnostics] = solveLinprog(A, b, obj)

% linprog solves
% min c' x     s.t.
%   Ax <= b
options = optimset('display','off','MaxIter',1e6) ;
[x,objOpt,exitflag,output,lambda] = linprog(obj,-A,-full(b),[],[],[],[],[],options) ;
duals = lambda.ineqlin ;
if exitflag == 1
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('fast:solveLinprog:notSolved','Error while solving problem using linprog. Exitflag %d',exitflag) ;
end
diagnostics.explanation = exitflag ;
diagnostics.objective = objOpt ;

