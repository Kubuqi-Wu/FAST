function [cntr, obj] = optimize_stock_nlds(scenario, var)

costP 	= 1 ;
costS 	= 2 ; 

x = var.x ;
s = var.s ;

if(scenario.getTime() == 1)
	cntr = (x>=0);
	obj = costP*x;
end

if(scenario.getTime() == 2)
	obj = -costS*s;
	notExceedStock = ( s <= x );
	notExceedDemand = ( s <= scenario.data );
	pos = ( s >= 0 );
	cntr = [notExceedDemand ; notExceedStock ; pos];
end

end