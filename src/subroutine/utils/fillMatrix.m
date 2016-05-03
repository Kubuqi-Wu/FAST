function matrix = fillMatrix(vector) 
if numel(vector) == 1
    matrix = (1:vector(end))' ;
else
    size = prod(vector) ; % rows
    l = vector(end) ; % columns        
    matrix = zeros(size, length(vector)) ;
    step = size / l ;
    for i = 1:l
        matrix((i-1)*step+1:i*step,end) = i ;
        mat = fillMatrix(vector(1:end-1)) ;
        matrix((i-1)*step+1:i*step,1:end-1) = mat ;
    end             
end
end
