function  [x, duals, objOpt, diagnostics] = solveGlpk(A, b, obj, param)

% Glpk solves
% min c' x      s.t.
%    Ax <= b

lb = [];
ub = [];
ctype = [];
vartype = [];
sense = 1;
% param.dual = 0;
% param.scale = 0;
% param.presol = 0;
% param.lpsolver = 2;
% param.dual = 2;
if(isempty(param))
    param = struct();
end
[x, objOpt, status, extra] = glpk(obj,-A,full(-b),lb,ub,ctype,vartype,sense,param);

if (status == 5)
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;    
    warning('fast:solveGlpk:notSolved','Error while solving problem using glpk.\nStatus and extras:') ;
    disp(status) ;
    disp(extra)
end
diagnostics.explanation = status ;
diagnostics.objective = objOpt ;

duals =  abs(extra.lambda); % Y U NOT POSITIVE? <- ?