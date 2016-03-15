classdef MILPModel % MILP Model
        
    % Models min c'x+k, Ax >= b, x(binaryIdx) binary, x(interIdx) integer
    properties
        A ;
        b ;
        c ;
        k ;
        binaryIdx ;   % ids of binary variable in A's domain
        integerIdx ;  % same for integer
        backMapVar ;  % original ids of the variables
        backMapCntr ; % original ids of the constraints        
    end
    
    methods
        
        %% Constructor
        function out = MILPModel(cntr, obj)   
            if nargin == 1
                obj = sddpConst(0) ;
            end
            if numel(obj) ~= 1
                error('MILPModel::MILPModel - obj need to be 1x1.') ;
            end
            % This gives A, b such that Ax <= b
            [out.A, out.b, out.binaryIdx, out.integerIdx, ...
                out.backMapVar, out.backMapCntr] = export(cntr) ;
            [cObj, out.k, backMapVarObj] = export(obj) ;
            if ~ all(ismember(backMapVarObj,out.backMapVar))
                error('MILPModel::MILPModel - Some of the variables of the objective are not present in the constraints.') ;
            end
            out.c = zeros(size(out.A, 2),1) ;
            out.c(forwardMapping(out.backMapVar, backMapVarObj)) = cObj ;                        
        end
        
        function primalValue = primalValue(self, variable, x)
            if any(size(x) ~= size(self.backMapVar))
                error('MILPModel::primalValue - Size mismatch') ;
            end
            primalValue = double(variable, self.backMapVar, x) ;
        end
        
        function dualValue = dualValue(self, cntr, lambda)
            if any(size(lambda) ~= size(self.backMapCntr))
                error('MILPModel::dualValue - Size mismatch') ;
            end
            dualValue = dual(cntr, self.backMapCntr, lambda) ;
        end  
                        
    end
   
end