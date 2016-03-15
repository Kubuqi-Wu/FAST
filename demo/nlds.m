function [cntr, obj] = nlds(scenario, x, s)
C 	= 1 ;
S 	= 2 ;

% It at stage 1, the problem is
%           min C x + E(V(x))
%           s.t.  x >= 0
% So we write the following constraints :
%           x >= 0
% and the objective function, without the E(V(x)) : the toolbox
% automatically takes care of that
%           min C * x
if(scenario.time == 1)
    % Constraint
    cntr = x >= 0 ;
    % Objective
    obj  = C * x ;
end

% If at stage 2, the problem is
%           min - P s
%           s.t.  s >= 0
%                 s <= d(xi)
%                 s <= x
% (by default we minimize, so the '-' is needed)
if(scenario.time == 2)
    % Constraints
    pos             = s >= 0 ;
    notExceedDemand = s <= scenario.data ; % data is the random demand store into each node (scenario)
    notExceedStock  = s <= x ;
    % Concatenate all constraints together
    cntr = [notExceedDemand ; notExceedStock ; pos];
    % The objective
    obj = -S * s;
end
end