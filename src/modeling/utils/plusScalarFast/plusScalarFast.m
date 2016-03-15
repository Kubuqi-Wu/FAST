function [ids, coefs] = plusScalarFast(id1, id2, coefs1, coefs2)
% Non-mex alternative for plusScalarFast

% [common, i1, i2] = intersect(id1, id2, 'stable') ;

% Alternative to intersect(..., ..., 'stable')
[common, i1, i2] = intersect(id1, id2);
[i1,idx_sort] = sort(i1);
common = common(idx_sort);
i2 = i2(idx_sort);

only1 = setdiff(1:numel(id1), i1) ;
only2 = setdiff(1:numel(id2), i2) ;
ids = [common ; id1(only1) ; id2(only2)] ;
coefs = [coefs1(i1)+coefs2(i2) ; coefs1(only1) ; coefs2(only2)] ;

end
