function out = demand(t,i)
% We don't have any randomnes at stage 1, no we don't store anything
if t == 1 
    out = [] ;
% At time 2, the demand is random, either 2 or 3, depending on the scenario:    
else 
    if i == 1
        out = 2 ;
    elseif i == 2
        out = 3 ;
    end
end
end