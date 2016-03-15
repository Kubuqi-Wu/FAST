%% Example of modeling
% Shows what already works
clc ; close all ; clear all ;

dummyVar = sddpVar(4,3,5) ; % Just to show that dummy variable dont impact anything
z = sddpVar(1, 3) ;
y = sddpVar(1, 3) ;
x = sddpVar(1, 3) ;
dummyVar2 = sddpVar(4,3,9) ; % Same reason
p = sddpVar() ;
q = sddpVar() ;
E = [1 2 4] ;
P = 4 ;

% Modeling stuff
cntr1 = x + 2.*y + 3*z >= E  ;
cntr2 = p == sum(y) ;
cntr3 = p <= q*P ;
dummyCntr = x+2.*y == E ; % Same reason
cntr = [cntr1 ;
        cntr2 ;        
        cntr3 ;
        setBinary(q) ;        
        x >= 0 ;
        y >= 0 ;
        z >= 0 ;
        p >= 0 ] ;
obj = z * [1 2 3]' ;

% Modeling -> matrix part
% This gives A, b such that Ax <= b
[A, b, binaryIdx, integerIdx, backMapVar, backMapCntr] = export(cntr) ;
% This gives the objective + constant term
[cp, k, ids] = export(obj) ;
% Since the above export does not know the rest of the variable, another
% mapping is necessary
c = zeros(1, size(A, 2)) ;
c(forwardMapping(backMapVar, obj.ids)) = cp ;

% We now have all we want
full(A)
b
c
binaryIdx

%% Let's now call Gurobi
gurobiModel.obj = c ;
gurobiModel.A = A ;
gurobiModel.rhs = b ;
gurobiModel.lb = - inf * ones(size(c)) ;
gurobiModel.sense = '<' ;
vtypes(1:size(A,2)) = 'C' ;
vtypes(binaryIdx) = 'C' ;% 'B' ; % Otherwise no dual variables
vtypes(integerIdx) = 'C' ;% 'I' ;
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
bVal = double(q, backMapVar, result.x) ;

% We also want dual variables :-)
cntr1Dual = dual(cntr1, backMapCntr, result.pi) ;
cntr2Dual = dual(cntr2, backMapCntr, result.pi) ;
cntr3Dual = dual(cntr3, backMapCntr, result.pi) ;

