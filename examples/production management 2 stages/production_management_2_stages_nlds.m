function [cntr, obj] = production_management_2_stages_nlds(scenario, x, s)

if(scenario.getTime() == 1)
	obj = x(1) + 2*x(2);
	cntr = x >= 0 ;
end

if(scenario.getTime() == 2)
	obj = - 2 * s(1) - 3*s(2);
	notExceedStock  = s <= x ;
	notExceedDemand = s <= scenario.data ;
	pos = s >= 0 ;
	cntr = [notExceedDemand ; notExceedStock ; pos];
end

end