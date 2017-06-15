function  [x, duals, objOpt, diagnostics] = solveMosek(A, b, obj, opts)

if(~ isempty(opts))
    error('fast:solveMosek:optsNonEmpty','non-empty opts not supported for Mosek');
end

options = mskoptimset; % default options
[cmd,verb,paramMosek] = msksetup(1,options);
[~,b,cmosek,l,u] = mskcheck('linprog',verb,length(obj),size(A),b,[0,0],[],[],[]);
[numineq,~] = size(A);
prob.a = -A;
prob.c = obj;
prob.blc    = [-inf*ones(size(b));cmosek];
prob.buc    = [-b;cmosek];
prob.blx    = l;
prob.bux    = u;

diagnostics = struct() ;

% Mosek can stalls when interior points is used, not with
% simplex
paramMosek.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_FREE_SIMPLEX';
[rcode,res] = mosekopt(cmd,prob,paramMosek);
res.sol.itr = res.sol.bas;  % Hack to avoid a bug in mskeflag :
% It displays a warning
% using res.sol.itr but this
% field does not exist.
exitflag = mskeflag(rcode,res);

if(exitflag ~= 1) % In general it works when simplex fails
    warning('fast:solveMosek:changingSolver','Unable to solve the problem. Trying interior-points.')
    paramMosek.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
    [rcode,res] = mosekopt(cmd,prob,paramMosek);
    exitflag = mskeflag(rcode,res);
    % do not type res.sol.itr = res.sol.bas;
    % because we are using interior-points.
end

if(exitflag ~= 1)
    % If we are here, I think we can give up ...
    warning('fast:solveMosek:changingSolver','Unable to solve the problem. Trying primal-dual simplex.')
    paramMosek.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX';
    [rcode,res] = mosekopt(cmd,prob,paramMosek);
    res.sol.itr = res.sol.bas;
    exitflag = mskeflag(rcode,res);
end

if(exitflag ~= 1)
    warning('fast:solveMosek:changingSolver','Unable to solve the problem. Trying primal simplex.')
    paramMosek.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_PRIMAL_SIMPLEX';
    [rcode,res] = mosekopt(cmd,prob,paramMosek);
    res.sol.itr = res.sol.bas;
    exitflag = mskeflag(rcode,res);
end

if(exitflag ~= 1)
    warning('fast:solveMosek:changingSolver','Unable to solve the problem. Trying dual simplex.')
    paramMosek.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_DUAL_SIMPLEX';
    [rcode,res] = mosekopt(cmd,prob,paramMosek);
    res.sol.itr = res.sol.bas;
    exitflag = mskeflag(rcode,res);
end

if(exitflag ~= 1)
    diagnostics.solved = false ;    
    warning('fast:solveMosek:notSolved','Error while solving problem using mosek. Exitflag %d',exitflag) ;
else
    diagnostics.solved = true ;
end
diagnostics.explanation = res ;
diagnostics.objective = res.sol.itr.pobjval ;


if ( isfield(res,'sol') )    
    if ( isfield(res.sol,'itr') )
        x = res.sol.itr.xx;
        objOpt = res.sol.itr.pobjval;
        duals = -res.sol.itr.y(1:numineq);
    else
        x = res.sol.bas.xx;
        objOpt = res.sol.bas.pobjval;
        duals = -res.sol.bas.y(1:numineq);
    end
else
    error('fast:solveMosek:notSolved','Error : field ''sol'' does not exist.');
end
