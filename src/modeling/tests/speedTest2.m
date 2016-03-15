%% Our
clc ; close all ; clear all ;

tic ;
x = sddpVar(20, 20) ;
y = sddpVar(20, 20) ;
z = sddpVar(20, 20) ;
p = sddpVar() ;
E = rand(20, 20) ;

% Model stuff
cntr1 = x + y + z <= E  ;
cntr2 = p == sum(y) ;
cntr = [cntr1 ;
        cntr2 ;        
        x >= 0 ;
        y >= 0 ;
        z >= 0 ;
        p >= 0 ] ;
    
[AA, bb, binary, integer, backMapVar, backMapCntr] = export(cntr) ;

% Objective
obj = z(:)' * [1:400]' ;
size(obj)
[cp, k, ids] = export(obj) ;
c = zeros(1, size(AA, 2)) ;
c(forwardMapping(backMapVar, obj.ids)) = cp ;

% Call gurobi
gurobiModel.obj = c ;
gurobiModel.A = AA ;
gurobiModel.rhs = bb ;
gurobiModel.lb = - inf * ones(size(c)) ;
gurobiModel.sense = '<' ;
gurobiParams.outputflag = 0 ;
result = gurobi(gurobiModel, gurobiParams) ;
if strcmpi(result.status,'OPTIMAL')
    diagnostics.solved = true ;
else
    diagnostics.solved = false ;
    warning('Error while solving problem using gurobi') ;
    disp(result) ;
end
toc ;

%% Yalmip
tic ;
clear all ;
x = sdpvar(20, 20) ;
y = sdpvar(20, 20) ;
z = sdpvar(20, 20) ;
p = sdpvar() ;
E = rand(20, 20) ;

% Model stuff
cntr1 = x + y + z <= E  ;
cntr2 = p == sum(y) ;
cntr = [cntr1 ;
        cntr2 ;        
        x >= 0 ;
        y >= 0 ;
        z >= 0 ;
        p >= 0 ] ;
    
obj = z(:)' * [1:400]' ;    
opt = sdpsettings('solver','mosek','verbose',0) ;    

optimize(cntr, obj, opt) ;
toc ;