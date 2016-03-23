function out = rainfall(t,i)
high = 10 ;
low = 2 ;
if t == 1
    out = (high + low)/2 ;
else
    if i == 1
        out = low ;
    else
        out = high ;
    end
end
end