classdef NldsModel
    % General problem:
    %       min c'x(t) + t + k
    % s.t.  Wx(t) >= h - Tx(t-1)
    %       x(t-1) fixed
    %      [t >= cuts]
    %       k is a constant
    
    properties
        c ; % [double] N x 1
        k ; % [double] 1 x 1                        
        W ; % [[double]] M x N
        h ; % [double] M x 1
        T ; % [[double]] M x P
        forwardIdx ; % [double]         -> the index in NldsModel space of the values to be sent to the next stage in the forward pass                
        modelVarIdx ; % [double] N x 1  -> the index of the variables in modeling space
        modelCntrIdx ; % [double] M x 1 -> the index of the constraints in modeling space
        cutCoeffs ; % [double] k x N    -> initially empty
        cutRHS ; % [double] k x 1       -> initially empty
    end
    
    methods
        
        function model = NldsModel(W,h,T,c,k,forwardIdx,modelVarIdx,modelCntrIdx)
            model.W = W ;
            model.c = c ;
            model.k = k ;
            model.h = h ;
            model.T = T ;
            model.forwardIdx = forwardIdx ;
            model.modelVarIdx = modelVarIdx ;
            model.modelCntrIdx = modelCntrIdx ;
            model.cutCoeffs = [];
            model.cutRHS = [];                        
        end
        
        function [nnzW,nnzT,nvarW,nvarT,nconstrW,nconstrT] = getSizeBlocks(model)
            nnzW = nnz(model.W);
            nnzT = nnz(model.T);
            [nconstrW,nvarW] = size(model.W);
            [nconstrT,nvarT] = size(model.T);
        end
        
        % TODO : to recode
        function [c,k,W,h,T] = getWSblocks(model,nvarTBlock, idxTBlock)
            c = model.c;
            k = model.k;
            h = model.h;
            W = model.W;
            nconstrTBlock = size(model.T,1);
            
            if(max(idxTBlock) > nvarTBlock)
                error('Unconsistent dimensions and index for block T.')
            end
            
            if(~isempty(model.T) )
                [i,j,val] = find(model.T);
                T = sparse(i,idxTBlock(j),val,nconstrTBlock, nvarTBlock,nnz(model.T));
            else
                T = [];
            end
        end
        
        % cutCoeffs  % [double] k x N
        % cutRHS  % [double] k x 1
        function model = addCut(model, cutCoeffs, cutRHS, params)
            % Check size
            if(~ iscolumn(cutRHS) || numel(cutRHS) ~= size(cutCoeffs, 1))
                error('Wrong size for the cuts to be added') ;
            end
            % Remove redondant cuts
            if params.algo.checkRedondance
                cutCoeffsKept = zeros(size(cutCoeffs)) ;
                cutRHSKept = zeros(size(cutRHS)) ;
                cutKept = 0 ;
                for i = 1:size(cutCoeffs, 1)
                    keepIt = true ;
                    for k = 1:size(model.cutCoeffs, 1)
                        if norm([cutCoeffs(i,:) cutRHS(i)] - [model.cutCoeffs(k,:) model.cutRHS(k)]) ...
                                < params.algo.redondanceTol * norm([cutCoeffs(i,:) cutRHS(i)]) % If cuts are too close
                            keepIt = false ;
                            break ;
                        end
                    end
                    if keepIt
                        cutKept = cutKept + 1 ;
                        cutCoeffsKept(cutKept, :) = cutCoeffs(i,:) ;
                        cutRHSKept(cutKept)       = cutRHS(i)      ;

                    end
                end  
                cutCoeffs  = cutCoeffsKept(1:cutKept,:) ;
                cutRHS     = cutRHSKept(1:cutKept) ;
            end
            if params.debug > 0
                fprintf('Adding cuts at node with coefs and rhs\n') ;
                disp(cutCoeffs) ;
                disp(cutRHS) ;
            end
            model.cutCoeffs = [model.cutCoeffs ; cutCoeffs];
            model.cutRHS = [model.cutRHS ; cutRHS];
        end
        
        function [A,b,obj,k] = getOptiProb(model,xTrial,withoutTheta,params)
            % Return the folowing optimisation problem :
            %
            %   min c'x+k   s.t.
            %   Ax >= b
            %
            % Basically, x = [x(t) t]
            % and        b = [h - Txtrial ; 
            %                 cuts ;
            %                 theta lower bound]
            %
            % Warning! xTrial is supposed to be the output of function
            % model.extract_xTrial(x), used with the model at the previous
            % stage.
            % xTrial can be empty
            
            % Warning : the constraint on theta must be always at the end,
            % because we remove it from the duals vector.
            
            if(~islogical(withoutTheta))
                error('Argument withoutTheta must be logical.')
            end
            
            if(withoutTheta)
                obj = [model.c ; 0]; % Because V(x) will be moved by mintheta.
            else
                obj = [model.c ; 1];
            end
                                  
            A = [model.W                    sparse(size(model.W,1),1)           ; ...
                 model.cutCoeffs            ones(size(model.cutCoeffs,1),1)     ; ...
                 zeros(1,size(model.W,2))   1                                   ] ;

            if ~ isempty(xTrial)
                b = [model.h - model.T*xTrial  ; ...
                     model.cutRHS              ; ...
                     params.algo.minTheta    ] ;
            else
                if ~ isempty(model.T)
                    error('model.T should be empty since xTrial is empty!') ;
                end
                b = [model.h              ; ...
                     model.cutRHS         ; ...
                     params.algo.minTheta ] ;
            end
            k = model.k ;
            
        end
        
        % The numerical values to be sent to the next node
        function xTrial = extractXTrial(model,x)
            xTrial = x(model.forwardIdx);
        end
        
        function [solution, diagnostics] = solve(model, solutionPreviousTime, withoutTheta, params)            
            if ~ isempty(solutionPreviousTime)
                [A,b,obj,k] = model.getOptiProb(solutionPreviousTime.trials,withoutTheta,params) ;
            else
                [A,b,obj,k] = model.getOptiProb([],withoutTheta,params) ;
            end       
            
            [x, duals, ~, diagnostics] = solveLP(A, b, obj, params);
            
            nVar    = size(model.W,2) ;
            nConstr = size(model.W,1) ;
            nDuals  = size(model.cutRHS,1) ;                                    
            
            solution.primal = x(1:nVar) ;
            solution.trials = model.extractXTrial(x) ;
            solution.costWithTheta = obj'*x + k ;
            solution.costWithoutTheta = solution.costWithTheta - obj(end)*x(end);
            solution.dualCntr = duals(1:nConstr) ;
            solution.dualCuts = duals(nConstr+1:nConstr+nDuals) ;
        end
        
        % TODO: recode
        function [solution, constraints, diagnostics] = solveUnderDecision(model, solutionPreviousTime, withoutTheta,solutionThisTime, params)
            warning('To recode') ;
            if ~ isempty(solutionPreviousTime)
                [A,b,obj] = model.getOptiProb(solutionPreviousTime.xTrial,withoutTheta,params) ;
            else
                [A,b,obj] = model.getOptiProb([],withoutTheta,params) ;
            end       
            
            n = length(model.forwardIdx);
            fixationMatrixRow = (1:(2*n))' ;
            fixationMatrixCol = [model.forwardIdx ; model.forwardIdx];
            fixationMatrixVal = [ones(n,1) ; -ones(n,1) ];
            fixationMatrix = sparse(fixationMatrixRow,fixationMatrixCol,fixationMatrixVal,2*n,size(A,2),2*n);
            A = [A ; fixationMatrix];
            b = [b ; solutionThisTime.xTrial ; -solutionThisTime.xTrial];
            
            [x, duals, objOpt, diagnostics] = solveLP(A, b, obj, params);
            
            solution.x = x ;
            solution.xTrial = model.extractXTrial(x) ;
            solution.costWithTheta = obj'*x;
            solution.costWithoutTheta = solution.costWithTheta - obj(end)*x(end);
            % Dual multipiers
            constraints = duals(1:end-1) ;                                    
        end   
        
        function primalSolution = updatePrimalSolution(model, primalSolution, variable, solution)
            if any(size(primalSolution) ~= size(variable))
                error('primalSolution shoudl have the same size as x variable') ;
            end
            if ~ isstruct(solution)
                error('solution should be a struct() with field primal') ;
            end
            newPrimalSolution = doubleBatch(variable, model.modelVarIdx, solution.primal) ;
            idxToUpdate = ~ isnan(newPrimalSolution) ;
            primalSolution(idxToUpdate) = newPrimalSolution(idxToUpdate) ;            
        end
        
        function [cutCoeffs, cutRHS] = getCutsCoeffs(model, variables)
            modelVarIdx = model.modelVarIdx ;
            varsVarIdx = [variables.ids]' ;
            cutRHS = model.cutRHS ;
            cutCoeffs = zeros(size(model.cutCoeffs)) ;
            if any(size(modelVarIdx) ~= size(varsVarIdx)) || any(sort(modelVarIdx) ~= sort(varsVarIdx))
                error('Variables should be the same than the one in the model, up to a reordering') ;
            end
            for i = 1:size(cutCoeffs, 2)
                j = find(modelVarIdx == varsVarIdx(i)) ;
                cutCoeffs(:, i) = model.cutCoeffs(:, j) ;
            end            
        end
                                                     
        function model = clearCuts(model)
            model.cutCoeffs = [];
            model.cutRHS = [];
        end       
        
    end
    
end
