function displayParams(params)
if params.verbose >= 1
    fprintf('                                   Param -                Value\n') ;
    displayParamsStruct(params,'') ;
end
end

function displayParamsStruct(params,prefix)
fields = fieldnames(params) ;
for i = 1:numel(fields)
    field = params.(fields{i}) ;
    if isstruct(field)
        displayParamsStruct(field,[prefix fields{i} '.']) ;
        continue ;
    elseif isa(field, 'function_handle')
        value = func2str(field) ;
    elseif ismatrix(field) || ischar(field)
        value = num2str(field) ;
    else
        error('Param type not supported') ;
    end
    fprintf('%40s - %20s\n',[prefix fields{i}],value) ;
end
end