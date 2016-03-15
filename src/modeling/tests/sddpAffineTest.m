% Tests
clc ; close all ; clear all ;

% Modeling
A = sddpVar(3, 2) ;
B = sddpVar(3, 2) ;
C = A + B ;
D = A - B ;
E = ones(size(A))+eye(size(A)) ;
F = E .* A ;
G = E * B ;
H = A ./ E ;
I = H + E ;

A
B

% Constraint
cntr1 = A <= B ;
cntr1

cntr2 = A == B ;
cntr2

cntr3 = [cntr1 ; cntr2] ;
cntr3

cntr4 = [cntr1 cntr2 ; cntr1 cntr2] ;
cntr4

cntr5 = A(1:2,:) <= B(1:2,:) ;
cntr5

cntr6 = [cntr1 ; cntr5] ;
cntr6

% Export constraints
[AA, bb, map] = export(A <= B) ;
AA
bb
map

[AA, bb, map] = export(A == B) ;
AA
bb
map

[AA, bb, map] = export(A <= B + E) ;
AA
bb
map
