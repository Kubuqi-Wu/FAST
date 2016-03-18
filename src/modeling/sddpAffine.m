classdef sddpAffine % Linear + Affine expression
    
    % This class is used to describe an affine expression
    % expression = [a1 a2 ... am] * [x1 x2 ... xm]' + b
    
    properties
        % constraint = coefs * var' + indep
        coefs ; % [n x 1] real
        indep ; % scalar
        ids   ; % [n x 1] integer, >= 1
    end
    
    methods
        
        %% Constructor
        function out = sddpAffine(coefs, indep, ids)
            % Coefs
            if nargin == 0
                coefs = 1 ;
                out.coefs = coefs ;
                nVar = 1 ;
            else
                nVar = numel(coefs) ;
                out.coefs = reshape(coefs, [nVar 1]) ; % Coefs
            end
            % Indep term
            if nargin <= 1
                indep = 0 ;
            end
            out.indep = indep ;
            % Id
            if nargin <= 2 % Fresh id
                ids = sddpAffine.getNewId(nVar) ;                
            end
            out.ids = ids ;
        end
        
        %% Export to c' x + k format, and return the indices
        
        function [c, k, idx] = export(self)
            if numel(self) ~= 1
                error('sddpAffine:export:notdefined','Export not defined for arrays of sddpAffine.') ;
            end
            c = reshape(self.coefs, [1 numel(self.coefs)]) ;
            k = self.indep ;
            idx = self.ids ;
        end
        
        %% Evaluate expression with a vector of value
        
        function out = eval(self, x)
            out = dot(self.coefs, x) + self.indep ;
        end
        
        %% Double
        
        % All the variables should be pressent in backwardMapping
        function out = double(self, backwardMapping, x)
            out = zeros(size(self)) ;
            for i = 1:numel(self)
                idx = forwardMapping(backwardMapping, self(i).ids) ;
                if ~ all(idx)
                    error('one of the variable is not present in backwardMapping') ;
                end
                out(i) = eval(self(i), x(idx)) ;
            end
        end
        
        function out = doubleBatch(self, backwardMapping, x)
            out = nan(size(self)) ;
            for i = 1:numel(self)
                idx = forwardMapping(backwardMapping, self(i).ids) ;
                if ~ all(idx)
                    % skip and leave nan
                else
                    out(i) = eval(self(i), x(idx)) ;
                end
            end            
        end
        
        %% Multiplication
        
        function out = times(in1, in2) % A .* B
            out = applyScalarNonScalarFunction(in1, in2, @timesScalar) ;
        end
        
        function out = rdivide(in1, in2)
            % in2 should be a real scalar or matrix
            if isreal(in2) && all(in2(:) ~= 0) && all(size(in1) == size(in2))
                out = times(in1, 1./in2) ; % Nice :-)
            else
                error('sddpAffine:rdivide:wrongInputArgument', 'Second argument should be a real non-zero matrix, with the same size as the first argument.') ;
            end
        end
        
        function out = timesScalar(in1, in2)
            if isreal(in1) && isa(in2, 'sddpAffine')
                out = in2 ;
                out.coefs = out.coefs * in1 ;
                out.indep = out.indep * in1 ;
            elseif isa(in1, 'sddpAffine') && isreal(in2)
                out = in1 ;
                out.coefs = out.coefs * in2 ;
                out.indep = out.indep * in2 ;
            else
                error('sddpAffine:timesScalar:wrongTypes', ...
                    'One argument should be a real number and the other a sddpAffine object.') ;
            end
        end
        
        function out = dot(in1, in2) % in1 and in2 should be vectors ; the orientation does not matter
            if ~ (isvector(in1) && isvector(in2))
                error('sddpAffine:dot:wrongSize1','The two arguments should be vectors, as specified by isvector') ;
            end
            if numel(in1) ~= numel(in2)
                error('sddpAffine:dot:wrongSize2','The two arguments should be have the same length') ;
            end
            if ~ (isreal(in1) || isreal(in2))
                error('sddpAffine:dot:wrongType','One of the two argument should be real') ;
            end
            
            if isreal(in1)
                inReal = in1 ;
                inAff = in2 ;
            else
                inReal = in2 ;
                inAff = in1 ;
            end
            allIds = inAff.collectIds ;
            coefs = zeros(size(allIds)) ;
            indep = 0 ;
            for i = 1:numel(in1)
                idsPos = forwardMapping(allIds, inAff(i).ids) ;
                coefs(idsPos) = coefs(idsPos) + inReal(i) * inAff(i).coefs ;
                indep = indep + inReal(i) * inAff(i).indep ;
            end
            out = sddpAffine(coefs, indep, allIds) ;
        end
        
        function out = mtimes(in1, in2) % A * B
            if isscalar(in1) || isscalar(in2)
                out = times(in1, in2) ;
            elseif length(size(in1)) == length(size(in2)) && length(size(in1)) == 2
                out = sddpVarDummy(size(in1, 1), size(in2, 2)) ;
                for i = 1:size(in1, 1)
                    for j = 1:size(in2, 2)
                        out(i,j) = dot(in1(i,:), in2(:,j)) ;
                    end
                end
            else
                error('sddpAffine:mtimes:wrongDimension', 'The two arguments should have sizes (m x p) and (p x n).') ;
            end
        end
        
        %% Addition
        
        function out = plus(in1, in2)
            out = applyScalarNonScalarFunction(in1, in2, @plusScalar) ;
        end
        
        function out = minus(in1, in2)
            out = plus(in1, times(in2, -1)) ; % Nice :-)
        end
        
        function out = uminus(in1)
            out = times(in1, -1) ;
        end
        
        function out = sum(in) 
            out = dot(in(:), ones(numel(in), 1)) ; % Faster than accumulated sums :-)
        end
        
        % out = in1 + in2 (e.g. (5x + 3y + 2) + (3z + 4))
        function out = plusScalar(in1, in2)
            if ~(isscalar(in1) && isscalar(in2))
                error('sddpAffine:plusScalar:wrongSizes', 'Both argument should be 1x1.') ;
            end
            if isa(in1,'sddpAffine') && isa(in2,'sddpAffine')
                %id1 = in1.ids ;
                %id2 = in2.ids ;
                %[idCommon1, idCommon2, idOnly1, idOnly2] = intersectFast(id1, id2) ; % Pre-compiled version
                %coefs =   [in1.coefs(idCommon1) + in2.coefs(idCommon2) ; in1.coefs(idOnly1)   ; in2.coefs(idOnly2)    ] ;
                %ids =     [in1.ids(idCommon1)                          ; in1.ids(idOnly1)     ; in2.ids(idOnly2)      ] ;
                [idsNew, coefsNew] = plusScalarFast(in1.ids, in2.ids, in1.coefs, in2.coefs) ;
                indepNew = in1.indep + in2.indep  ;
                out = sddpAffine(coefsNew, indepNew, idsNew) ;
            elseif isa(in1,'sddpAffine') && isreal(in2)
                out = in1 ;
                out.indep = out.indep + in2 ;
            elseif isa(in2, 'sddpAffine') && isreal(in1)
                out = in2 ;
                out.indep = out.indep + in1 ;
            else
                error('sddpAffine:plus:wrongTypes', 'At least one of the argument should be a sddpAffine variable. If only one, the other should be a real scalar.') ;
            end
        end
        
        %% Affine constraints
        
        function out = eq(in1, in2)
            out = sddpConstraint(in1, in2, '==') ;
        end
        
        function out = le(in1, in2)
            out = sddpConstraint(in1, in2, '<=') ;
        end
        
        function out = ge(in1, in2)
            out = sddpConstraint(in1, in2, '>=') ;
        end
        
        %% Integrality constraints
        
        function out = setInteger(in)
            if ~ isVariableOnly(in)
                error('sddpAffine:setInteger','isInteger support variable-argument only ; no expressions allowed.') ;
            end
            out = sddpConstraint(in, [], 'integer') ;
        end
        
        function out = setBinary(in)
            if ~ isVariableOnly(in)
                error('sddpAffine:setBinary','isBinary  support variable-argument only ; no expressions allowed.') ;
            end
            out = sddpConstraint(in, [], 'binary') ;
        end
        
        function yesOrNot = isVariableOnly(in)
            yesOrNot = true ;
            for i = 1:numel(in)
                if numel(in(i).coefs) > 1 || in(i).indep ~= 0
                    yesOrNot = false ;
                    break ;
                end
            end
        end
        
        %% Auxiliary methods
        
        function out = applyScalarNonScalarFunction(in1, in2, func)
            s1 = size(in1) ;
            s2 = size(in2) ;
            if numel(s1) == numel(s2) && all(s1 == s2) % 1 and 2 non-scalar
                out = sddpVarDummy(s1) ;
                for i = 1:prod(s1)
                    out(i) = func(in1(i), in2(i)) ;
                end
            elseif isscalar(in1) && ~isscalar(in2) % 1 scalar, 2 non-scalar
                out = sddpVarDummy(s2) ;
                for i = 1:prod(s2)
                    out(i) = func(in1, in2(i)) ;
                end
            elseif ~isscalar(in1) && isscalar(in2) % 1 non-scalar, 2 scalar
                out = sddpVarDummy(s1) ;
                for i = 1:prod(s1)
                    out(i) = func(in1(i), in2) ;
                end
            else
                error('sddpAffine:func:wrongSizes', 'Input argument should either have the same size, or one at least should be a scalar (1x1).') ;
            end
        end
        
        function out = isVariable(self)
            if numel(self.ids) ~= 1
                out = false ; 
                return ;
            elseif self.indep ~= 0
                out = false ;
                return ;
            elseif self.coefs ~= 1
                out = false ;
                return ;
            else
                out = true ;
                return ;
            end
        end
        
        function disp(self)
            s = size(self) ;
            str = num2str(s(1)) ;
            for i = 2:length(s)
                str = [str 'x' num2str(s(i))] ;
            end
            fprintf('%s array of affine expressions\n',str) ;
        end
        
        function idsUnique = collectIds(self)
            ids = {self.ids} ; % Collect ids in a n1*n2*...*nn cell, where each cell contains an array of indices
            idsArray = mergeCells(ids) ;
            % We now have an array with all unique ids
            idsUnique = unique(idsArray) ;
        end
        
        function s = toString(self)
            if numel(self) == 1
                s = toStringScalar(self) ;
            else
                s = cell(size(self)) ;
                for i = 1:numel(self)
                    s{i} = sprintf('%s\n', toStringScalar(self(i))) ;
                end
            end
        end
        
        function s = toStringScalar(self)
            if numel(self) ~= 1
                error('sddpAffine:plot:notdefined', 'toStringScalar is not defined for arrays of sddpAffine') ;
            end
            s = sprintf('%2.2f x%d', self.coefs(1), self.ids(1)) ;
            for i = 2:length(self.coefs)
                s = sprintf('%s + %2.2f x%d', s, self.coefs(i), self.ids(i)) ;
            end
            s = sprintf('%s + %2.2f\n', s, self.indep) ;
        end
        
    end
    
    methods (Static = true)
        
        %% Auxiliary constructor
        function out = createArray(dimensions)  
            % Because repmat(...,3) will create a 3x3 array ...
            if isscalar(dimensions)
                dimensions = [dimensions(1);1] ;
            end
            out = repmat(sddpAffine(1,0,-1),dimensions');
            newIds = sddpAffine.getNewId(numel(out)) ;
            % Allocate id
            for i = 1:numel(out)
                out(i).ids = newIds(i) ;
            end                        
        end
        
        function out = createConstVariable(k)
            out = sddpAffine(double.empty(0,1),k,double.empty(0,1)) ;
        end
        
    end
    
    methods (Access = private, Static = true)
        function out = getNewId(nbr) % Returns an integer from 1 to oo
            % -1 is the 'trash' integer, for internal
            % use only
            if nargin == 0
                nbr = 1 ;
            end
            persistent idMax;
            if isempty(idMax)
                idMax = 1 ;
            end
            out = zeros(nbr, 1) ;
            out(1:nbr) = idMax + ((1:nbr)-1) ;
            idMax = out(end)+1 ;            
        end
    end
    
end