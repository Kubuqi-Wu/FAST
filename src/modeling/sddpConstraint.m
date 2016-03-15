classdef sddpConstraint
    
    % The idea is that everything is stored as entered, and is converted
    % just before the 'export' called.
    % This allow better checks from the user, if needed
    properties
        lhs ; % An expression (or an array of -)
        rhs ; % Another expression (or an array of -) or []
        type ; % The meaning of                % may be :
               %    '=='
               %    '<='
               %    '>='the inequality
               %    'integer'
               %    'binary'
        id ; % A unique ID

    end
    
    methods
        
        %% Constructor
        
        function out = sddpConstraint(lhs, rhs, type)
            out.lhs = lhs ;
            out.rhs = rhs ;
            out.type = type ;
            out.id = sddpConstraint.getNewId() ;
        end   
                       
        %% At the end, for the 3 first types, we want <= inequalities
                
        % Return A, b and forwardMapping such that A*x >= b
        % where backwardMappingVariables map new and old variables indexes       
        % and where backwardMappingCntrs map new and old constraints
        % indexes
        function [A, b, binary, integer, ...
                  backMapVars, backMapCntrs] = export(self)
                    
            
            % Linear constraints
            idsCells = cell(numel(self), 1) ;
            idCntrCells = zeros(numel(self), 1) ;
            exprs = cell(numel(self), 1) ; 
            nExprs = 0 ; 
            nCntr = 0 ;
            % Binary constraints
            idBinarys = cell(numel(self),1) ; nBinarys = 0 ;
            % Integers constraints
            idIntegers = cell(numel(self),1) ; nIntegers = 0 ;
            % Collect everything and get indices
            for i = 1:numel(self)
                switch self(i).type
                    case '<='
                        nExprs = nExprs + 1 ;
                        exprs{nExprs} = self(i).rhs - self(i).lhs ;
                        idsCells{nExprs} = exprs{nExprs}.collectIds() ;
                        idCntrCells(nExprs) = self(i).id ;
                        nCntr = nCntr + numel(exprs{nExprs}) ;
                    case '>='
                        nExprs = nExprs + 1 ;                        
                        exprs{nExprs} = self(i).lhs - self(i).rhs ;
                        idsCells{nExprs} = exprs{nExprs}.collectIds() ;
                        idCntrCells(nExprs) = self(i).id ;
                        nCntr = nCntr + numel(exprs{nExprs}) ;                        
                    case '=='
                        nExprs = nExprs + 1 ;                        
                        exprs{nExprs} = [self(i).lhs - self(i).rhs ;
                                         self(i).rhs - self(i).lhs ] ;
                        idsCells{nExprs} = exprs{nExprs}.collectIds() ;
                        idCntrCells(nExprs) = self(i).id ;
                        nCntr = nCntr + numel(exprs{nExprs}) ;                                                
                    case 'binary'                        
                        idBinarys{nBinarys+1} = self(i).lhs.collectIds() ;
                        nBinarys = nBinarys + 1 ;
                    case 'integer'                        
                        idIntegers{nBinarys+1} = self(i).lhs.collectIds() ;
                        nIntegers = nIntegers + 1 ;
                end
                
            end          

            binaryOld = sort(unique(mergeCells(idBinarys))) ;
            integerOld = sort(unique(mergeCells(idIntegers))) ; 
            % We merge linear constraints to have the unique variables
            affineOld = mergeCells(idsCells) ;
            % We add binary and integer
            idsArrayOld = [affineOld ; binaryOld ; integerOld] ;
            backMapVars = sort(unique(idsArrayOld)) ;
            % Find binary and integer variables in the forwardMapping
            binary = zeros(size(binaryOld)) ;
            for i = 1:numel(binary)
                binary(i) = find(binaryOld(i) == backMapVars) ;
            end
            integer = zeros(size(integerOld)) ;
            for i = 1:numel(integer)
                integer(i) = find(integerOld(i) == backMapVars) ;
            end
            % Prepare backward mapping of affine constraints
            backMapCntrs = zeros(nCntr, 1) ;
            
            % We count the number of variables and constraints
            nVar = numel(backMapVars) ;
            A = sparse(nCntr, nVar) ;          
            
            
            % Secondly, we iterate over all the constraints and build the
            % matrix
            % This could be done more efficiently
            idCntr = 1 ;
            for i = 1:numel(self)
                exprCurrent = exprs{i} ;  
                for j = 1:numel(exprCurrent)
                    exprCntr = exprCurrent(j) ;
                    idsOld = exprCntr.ids ;
                    idsVar = forwardMapping(backMapVars, idsOld) ;
                    A(idCntr, idsVar) = exprCntr.coefs ;
                    b(idCntr, 1) = - exprCntr.indep ; % Because the b is on the right here
                    backMapCntrs(idCntr) = idCntrCells(i) ;
                    idCntr = idCntr + 1 ;
                end                                
            end                                   
        end
        
        function out = dual(self, backwardMapping, lambda)
            if numel(self) == 1
                out = dualScalar(self, backwardMapping, lambda) ;
            else
                error('sddpConstraint:dual:notSupported','dual only handle scalar input (i.e. not an array of sddpConstraint)') ;
            end            
        end
        
        function out = dualScalar(self, backwardMapping, lambda)
            if numel(self) ~= 1
                error('sddpConstraint:dualScalar:wrongSize','self should be scalar') ;
            end
            
            out = zeros(size(self.lhs)) ;  
            idsCheck = backwardMapping == self.id ;
            if strcmp(self.type, '<=') 
                out(:) = - lambda(idsCheck(:)) ;    
            elseif strcmp(self.type, '>=')
                out(:) =  lambda(idsCheck(:)) ;            
            elseif strcmp(self.type, '==')
                L = lambda(idsCheck) ;
                L1 = L(1:end/2) ; % corresponds to the modeling of >=
                L2 = L(end/2+1:end) ; % <=
                out(:) = L1 - L2 ; % We take 1 out of 2
            else
                error('sddpConstraint:dualScalar:notSupported', 'dual not supported for integrality constraints') ;
            end           
        end                                
        
        %% Auxiliary functions
        
        function disp(self)
            sizeSelf = size(self) ;
            size1 = num2str(sizeSelf(1)) ;
            for i = 2:length(sizeSelf)
                size1 = [size1 'x' num2str(sizeSelf(i))] ;
            end
            sizeFirst = size(self(1,1).lhs) ;
            size2 = num2str(sizeFirst(1)) ;
            for i = 2:length(sizeFirst)
                size2 = [size2 'x' num2str(sizeFirst(i))] ;
            end
            strType = self(1,1).type ;
            strSize = ['size ' size2] ;
            for i = 2:numel(self)
                if ~strcmp(self(i).type,strType)
                    strType = 'mixed-type' ;                    
                end
                if any(size(self(i).lhs) ~= sizeFirst)
                    strSize = 'mixed-size' ;
                end
            end            
            fprintf('%s array of %s constraint of %s\n',size1,strType,strSize) ;
        end
        
    end
    
    methods (Access = private, Static = true)
        function out = getNewId() % Returns an integer from 1 to oo
            % -1 is the 'trash' integer, for internal
            % use only
            persistent idCntr;
            if isempty(idCntr)
                idCntr = 1 ;
            else
                idCntr = idCntr + 1 ;
            end
            out = idCntr ;
        end
    end
    
end