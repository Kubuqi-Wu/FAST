function out = rainfall(t,i)
high = 10 ;
low = 2 ;
if t == 1
    out = (high + low)/2 ;
else
    out = (i == 1) * low + (i == 2) * high ;
end
end