function out = rainfallUnif(t,i)
min = 1 ;
max = 20 ;
rainUnif = min:max ;
if t == 1
    out = mean(rainUnif) ;
else
    out = rainUnif(i) ;
end
end