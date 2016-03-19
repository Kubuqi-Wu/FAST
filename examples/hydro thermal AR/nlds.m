function [cntr, obj] = nlds(scenario, var)

t = scenario.getTime() ;

x = var.x ;
p = var.p ;
q = var.q ;
r = var.r ;
l = var.l ;

if t > 1
    xt1 = x(t-1) ;
    rt1 = r(t-1) ;
else    
    xt1 = 0 ;
    rt1 = 0 ;
end

obj = 1000 * l(t) + 25 * p(t) ;
cntr = [ l(t) + p(t) + q(t) >= 1000 ;
         p(t) <= 500 ;
         q(t) <= xt1 + r(t) ;
         x(t) == xt1 + r(t) - q(t) ;
         x(t) <= 1000 ;
         r(t) <= rt1 + scenario.data ;
         l(t) >= 0 ;
         p(t) >= 0 ] ;
     
end

