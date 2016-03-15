%% Small example
clear all ;
x = sddpVar() ;
y = sddpVar() ;
z = sddpVar() ;
p = sddpVar() ;

a = 0.33 * ([x y] * [1 ; 2] + 1);

a.toString


A = [1 2 ; 3 4] ;
b = A * a ;
s = b.toString ;
s{:}
    
clc ;

x = sddpVar(2, 1) ;
y = A*x ;
s = y.toString ;
s{:}

A = rand(500, 500) ;
x = sddpVar(500, 1) ;
tic ;
y = A*x ;
toc ;


x = sddpVar(1, 1000) ;
y = rand(1000, 1) ;

% Method 1
tic ;
z1 = dot(x, y) ;
toc ;

% Method 2
tic ;
z2 = x * y ;
toc ;
