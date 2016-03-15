function out = optimize_stock_demand(t,i)
if t == 1
    out = [] ;
else
    if(i == 1)
        out = 2;
    elseif (i == 2)
        out = 3;
    elseif (i == -1) % for expected lattice
        out = 2.5;
    end
end
end
