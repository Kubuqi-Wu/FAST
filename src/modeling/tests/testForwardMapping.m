clc ; close all ; clear all 

k = 0 ;
while true

    backward = (randi(100, [150, 1])) ;
    ids = (randi(100, [30, 1])) ;
    tic ;
    forward = forwardMapping(backward, ids) ;
    time = toc ;
    
    forwardCheck = zeros(size(ids)) ;
    for i = 1:length(forwardCheck)
        res = find(backward == ids(i)) ;
        if ~isempty(res)
            forwardCheck(i) = res(1) ;
        end
    end
    k = k+1 ;
    
    if any(forward ~= forwardCheck)
        error('Error !') 
    elseif mod(k, 100) == 0
        fprintf('%d - %f s. - Ok\n', k, time) ;
    end
    
end
