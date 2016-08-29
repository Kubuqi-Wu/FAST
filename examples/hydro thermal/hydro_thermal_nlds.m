function [cntr, obj] = hydro_thermal_nlds(scenario, var)

t = scenario.getTime() ;

C = 5 ;
V = 8 ;
d = 6 ;

x = var.x ;
y = var.y ;
p = var.p ;

% Fuel cost
fuel_cost = C * p(t) ;
% Maximum reservoir level
reservoir_max_level = x(t) <= V ;
% Meet demand
meet_demand = p(t) + y(t) >= d ;
% Positivity
positivity = [x(t) >= 0, y(t) >= 0, p(t) >= 0] ;
% Take rain into account
if t == 1
    reservoir_level = x(1) + y(1) <= scenario.data ;
else
    reservoir_level = x(t) - x(t-1) + y(t) <= scenario.data ;
end

obj = fuel_cost ;
cntr = [reservoir_max_level, ...
        meet_demand, ...
        positivity, ...
        reservoir_level] ;
end

