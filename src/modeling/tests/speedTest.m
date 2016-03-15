function test1
clc ; close all ; clear all ;

% Our method

tic ;
x = sddpVar(10, 10, 10, 10) ;
y = sddpVar(10, 10, 10, 10) ;
z = sddpVar(10, 10, 10, 10) ;
toc ;
E = rand(10, 10, 10, 10) ;
p = sddpVar(1) ;
F = rand(1) ;

tic ;
expr1 = x + y + z  ; 
expr2 = x + y + z ;
expr3 = F + p ;
expr4 = sum(expr2) ;
cntr1 = expr1 <= 0 ;
cntr2 = expr1 >= E ;
cntr3 = expr4 <= expr3 ;
cntr = [cntr1 ; cntr2 ; cntr3] ;
toc ;

tic ;
[A, b, binary, integer, mapVar, mapCntr] = export(cntr) ;
toc ;

% % Yalmip
% clear all ;
% close all ;
% tic ;
% x = sdpvar(10, 10, 10, 10) ;
% y = sdpvar(10, 10, 10, 10) ;
% z = sdpvar(10, 10, 10, 10) ;
% toc ;
% E = rand(10, 10, 10, 10) ;
% p = sdpvar(1) ;
% F = rand(1) ;
% 
% tic ;
% expr1 = x + y + z  ; 
% expr2 = x + y + z ;
% expr3 = F + p ;
% expr4 = sum(expr2) ;
% cntr1 = expr1 <= 0 ;
% cntr2 = expr1 >= E ;
% cntr3 = expr4 <= expr3 ;
% cntr = [cntr1 ; cntr2 ; cntr3] ;
% toc ;

