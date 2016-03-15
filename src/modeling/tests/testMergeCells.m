% test merge cells
clc ; close all ; clear all ;


while true
    
    m = 10 ;
    n = 10 ;
    C = cell(m, n) ;
    for i = 1:m
        for j = 1:n
            C{i, j} = randi(100, [randi(10) randi(10)]) ;
        end
    end
    
    % Our function
    tic ;
    merged = mergeCells(C) ;
    time1 = toc ;
    
    % Stupid implementation
    tic ;
    mergedCheck = [] ;
    for i = 1:numel(C)
        data = C{i} ;
        mergedCheck = [mergedCheck ; data(:)] ;
    end
    time2 = toc ;
    
    if any(merged ~= mergedCheck)
        error('Error !') ;
    else
        fprintf('%f vs %f\n', time1, time2) ;
    end
    
end


