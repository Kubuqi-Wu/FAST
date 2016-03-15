function merged = mergeCells(C)
% Non-mex alternative for mergeCells
merged = [] ;
for i = 1:numel(C)
    data = C{i} ;
    merged = [merged ; data(:)] ;
end
end