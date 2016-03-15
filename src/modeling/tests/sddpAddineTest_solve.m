% Test solve
clc ; close all ; clear all ;

a = sddpVar(5, 5) ; % Useless but useful to check indexing
z = sddpVar(1, 3) ;
b = sddpVar() ; % Same
y = sddpVar(1, 3) ;
x = sddpVar(1, 3) ;
p = sddpVar() ;
b = sddpVar() ;
E = rand(1, 3) ;
P = rand(1, 1) ;

% Model stuff
cntr1 = x + 2.*y + 3*z <= E  ;
cntr2 = p == sum(y) ;
cntr3 = p <= b*P ;
cntr = [cntr1 ;
        cntr2 ;        
        cntr3 ;
        isBinary(b) ;        
        x >= 0 ;
        y >= 0 ;
        z >= 0 ;
        p >= 0 ] ;
    
[AA, bb, binary, integer, backMapVar, backMapCntr] = export(cntr) ;

% Objective
obj = z * [1 2 3]' ;

toString(obj)
[cp, k, ids] = export(obj) ;
c = zeros(1, size(AA, 2)) ;
c(forwardMapping(backMapVar, obj.ids)) = cp ;

% Call gurobi
gurobiModel.obj = c ;
gurobiModel.A = AA ;
gurobiModel.rhs = bb ;
gurobiModel.lb = - inf * ones(size(c)) ;
gurobiModel.sense = '<' ;
vtypes(1:size(AA,2)) = 'C' ;
vtypes(binary) = 'C' ;% 'B' ;
vtypes(integer) = 'C' ;% 'I' ;
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

% We want the results :-)
xVal = double(x, backMapVar, result.x) ;
yVal = double(y, backMapVar, result.x) ;
zVal = double(z, backMapVar, result.x) ;
pVal = double(p, backMapVar, result.x) ;
bVal = double(b, backMapVar, result.x) ;

% We also want dual variables :-)
cntr1Dual = dual(cntr1, backMapCntr, result.pi) ;
cntr2Dual = dual(cntr2, backMapCntr, result.pi) ;
cntr3Dual = dual(cntr3, backMapCntr, result.pi) ;
