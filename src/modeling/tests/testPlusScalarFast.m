% test [ids, coefs] = plusScalarFast(id1, id2, coefs1, coefs2)

clc ; close all ; clear all ;

i = 1 ;

while true ;
    id1 = unique(randi(1000000, [10 1])) ;
    id2 = unique(randi(1000000, [8 1])) ;
    coefs1 = rand(size(id1)) ;
    coefs2 = rand(size(id2)) ;
        
    tic ;
    [ids, coefs] = plusScalarFast(id1, id2, coefs1, coefs2) ;
    time1 = toc ;

    tic ;
    [common, i1, i2] = intersect(id1, id2, 'stable') ;
    only1 = setdiff(1:numel(id1), i1) ;
    only2 = setdiff(1:numel(id2), i2) ;
    idsCheck = [common ; id1(only1) ; id2(only2)] ;
    coefsCheck = [coefs1(i1)+coefs2(i2) ; coefs1(only1) ; coefs2(only2)] ;
    time2 = toc ;
    
    if ~ (all(ids == idsCheck) && all(coefs == coefsCheck))
        error('Error !') ;
    elseif mod(i, 100) == 0
        fprintf('%d - %f vs %f\n', i, time1, time2) ;
    end
    i = i+1 ;
end
