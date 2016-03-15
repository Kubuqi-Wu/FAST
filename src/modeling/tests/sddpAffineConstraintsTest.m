%% Small example
clear all ;
x = sddpVar() ;
y = sddpVar() ;
z = sddpVar() ;
p = sddpVar() ;

cntr = [x + y <= 0 ; ...
        y - 3.*z >= 1 ;
        x - p == 3] ;
    
% The corresponding A * x <= b matrix is
Acheck = [1  1  0 0  ;
          0  -1 3 0  ;
          -1 0  0 1  ;
          1  0  0 -1 ] ;
bcheck = [0 -1 -3 3]' ;

[A, b, mapping] = export(cntr) ;

Acheck
A

bcheck
b