% Test intersect fast
clc ; close all ; clear all ;
k = 0 ;
while true    
    id1 = unique(randi(30, 30, 1)) ;
    id2 = unique(randi(20, 25, 1)) ;  
    tic ;
    [~, idCommon1Check, idCommon2Check] = intersect(id1, id2) ;
    idOnly1Check = setdiff(1:length(id1), idCommon1Check)' ;
    idOnly2Check = setdiff(1:length(id2), idCommon2Check)' ;
    time1 = toc ;
    tic ;
    [idCommon1, idCommon2, idOnly1, idOnly2] = intersectFast(id1, id2) ;  
    time2 = toc ;
    
    check1 = all(idCommon1Check == idCommon1) ;
    check2 = all(idCommon2Check == idCommon2) ;
    check3 = all(idOnly1Check == idOnly1) ;
    check4 = all(idOnly2Check == idOnly2) ;
    
    if ~ (check1 && check2 && check3 && check4)
        error('Error') ;
    else
        k = k+1 ;
        fprintf('%d - %f vs %f\n', k, time1, time2) ;
    end
 end



