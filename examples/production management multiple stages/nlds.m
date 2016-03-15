function [cntr, obj] = nlds(scenario, var, H, productionCost, salePrice)

x = var.x ; 
s = var.s ; 
t = scenario.getTime() ;
i = scenario.getIndex() ;

demand = (i == 1) * 2 + (i == 2) * 10 ;   

if(t == 1)
	obj = productionCost' * x(:, 1) ;
	cntr1 = x(:, 1) >= 0 ;
    cntr = cntr1 ;
elseif(t > 1 && t < H)    
    obj = productionCost' * x(:, t) - salePrice' * s(:, t-1) ; 
    cntr1 = x(:, t)   >= 0 ;
    cntr2 = s(:, t-1) >= 0 ;
    cntr3 = s(:, t-1) <= x(:, t-1) ;
    cntr4 = sum(s(:, t-1)) <= demand ;
    cntr = [cntr1 cntr2 cntr3 cntr4] ;
elseif(t == H)
    obj = - salePrice' * s(:, H-1) ;
    cntr1 = s(:, H-1) >= 0 ;
    cntr2 = s(:, H-1) <= x(:, H-1) ;
    cntr3 = sum(s(:, H-1)) <= demand ;
    cntr = [cntr1 cntr2 cntr3] ;
else
    error('t not consistent') ;
end

end
