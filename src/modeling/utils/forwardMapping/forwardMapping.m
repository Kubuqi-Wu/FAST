function forward = forwardMapping(backward, ids)
% Non-mex alternative for forwardMapping

forward = zeros(size(ids)) ;
for i = 1:length(forward)
    res = find(backward == ids(i)) ;
    if ~isempty(res)
        forward(i) = res(1) ;
    end
end