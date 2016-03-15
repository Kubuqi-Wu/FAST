function out = demand(t,i,demandVector,nNodes)
if t == 1
    out = [] ;
else
    if(i >= 1)
        out = demandVector(i);
    elseif(i == -1)
        out = mean(demandVector(1:nNodes));
    else
        error('Invalid idx')
    end
end
end