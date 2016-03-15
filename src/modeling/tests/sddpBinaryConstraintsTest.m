%% Small example
clear all ;
x = sddpVar() ;
y = sddpVar() ;
z = sddpVar() ;
p = sddpVar() ;

cntr1 = isInteger(x) ;
cntr2 = isBinary(y) ;
cntr3 = x + y <= 0 ;

cntr1
cntr2
cntr3

cntr = [cntr1 ; cntr2 ; cntr3] ;
[A1, b1, binary1, integer1] = cntr.export() ;
A1
b1
integer1
binary1

cntr = [x+y >= 1 ; isInteger(z)] ;
[A2, b2, binary2, integer2] = cntr.export() ;
A2
b2
integer2
binary2

cntr = [x+y <= 1 ; 
        x+z == 2 ;
        x+3.*y >= 3 ;
        isBinary(x) ;
        isInteger(p) ]

[A3, b3, binary3, integer3] = cntr.export() ;
A3
b3
integer3
binary3
    