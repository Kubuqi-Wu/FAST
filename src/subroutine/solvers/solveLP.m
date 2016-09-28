function [x, duals, objOpt, diagnostics] = solveLP(A, b, c, params)
% SolveLP Solves a Linear Program
%
% [x, duals, objOpt, diagnostics] = SOLVELP(A, b, c, params) solves the
%   linear program 
%       min c'x s.t.
%           Ax >= b
%   using the solver params.solver ('gurobi','linprog','mosek','glpk','cplex').
%   x is the optimal primal solution
%   duals is the dual solution
%   objOpt is the objective value at the optimal solution
%   diagnostics is the diagnostic output by the solver.
%   
%   The function is most likely to crash in case of infeasibility or
%   unboundedness.

% Solve :
%   min c'x s.t.
%   Ax >= b

if size(b, 2) ~= 1 || size(c, 2) ~= 1
    error('fast:solveLP:notColumn','b and c should be column vectors.') ;
end
if size(A, 1) ~= size(b, 1)
    error('fast:solveLP:dimensionsMismatch', 'A and b should have the same number of rows') ;
end
if size(A, 2) ~= size(c, 1)
    error('fast:solveLP:dimensionsMismatch', 'A and c should have the same number of columns') ;
end

opts = params.solverOpts;

if strcmpi(params.solver,'linprog')
    [x, duals, objOpt, diagnostics] = solveLinprog(A, b, c, opts) ;
elseif strcmpi(params.solver,'gurobi')    
    [x, duals, objOpt, diagnostics] = solveGurobi(A, b, c, opts) ;
elseif strcmpi(params.solver,'mosek')    
    [x, duals, objOpt, diagnostics] = solveMosek(A, b, c, opts) ;
elseif strcmpi(params.solver,'glpk')    
    [x, duals, objOpt, diagnostics] = solveGlpk(A, b, c, opts) ;
elseif strcmpi(params.solver,'cplex')    
    [x, duals, objOpt, diagnostics] = solveCPLEX(A, b, c, opts) ;
else
    error('fast:solveLp:nonSupportedSolver','Non supported solver !') ;
end             
% Sanity check
if abs(objOpt - c'*x) > 1e-6 * abs(c'*x)   
    warning('fast:solveLp::objectiveNumericalError','Objective seems different from obj''x. Check why.') ;
end    
% Warning if not solved
if ~ diagnostics.solved
    warning('fast:solveLp:notSolved','LP not solved ! Explanation (as from the solver):') ;
    disp(diagnostics.explanation) ;
end
% Returning solution
if params.debug > 0
    fprintf('Solution :\n') ;
    disp(x) ;
end
