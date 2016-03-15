classdef ModelingTests < matlab.unittest.TestCase
    
    methods (Test)
        
        % Full problem
        function testResults1(testCase)
            
            % Solve a simple LP
            % Check model, result and dual
            
            A = [1 2 3 4 ; 7 8 9 10] ;
            b = [1 ; 4] ;
            x = sddpVar(4,1) ;
            % Min norm sol ==> min sum x_i, x >= 0, Ax == b
            cntr1 = x >= 0 ;
            cntr2 = A*x == b ;
            cntr = [cntr1 ; cntr2] ;
            obj = sum(x)+2 ;
            milpModel = MILPModel(cntr, obj) ;
            
            % Check matrices
            Acheck = sparse([eye(4) ; A ; -A]) ;
            bcheck = [zeros(4,1) ; b ; -b] ;
            ccheck = [1 1 1 1]' ;
            kcheck = 2 ;
            testCase.verifyEqual(milpModel.A,Acheck) ;
            testCase.verifyEqual(milpModel.b,bcheck) ;
            testCase.verifyEqual(milpModel.c,ccheck) ;
            testCase.verifyEqual(milpModel.k,kcheck) ;
            
            % Check indexes
            backMapVarCheck = [1 2 3 4]' ;
            backMapCntrCheck = [1 1 1 1 2 2 2 2]' ;
            testCase.verifyEqual(milpModel.backMapVar,backMapVarCheck) ;
            testCase.verifyEqual(milpModel.backMapCntr,backMapCntrCheck) ;
            
            % Check integer & co
            binaryCheck = double.empty(0,1) ;
            integerCheck = double.empty(0,1) ;
            testCase.verifyEqual(milpModel.binaryIdx,binaryCheck) ;
            testCase.verifyEqual(milpModel.integerIdx,integerCheck) ;
            
            % Run Gurobi
            gurobiModel.obj = milpModel.c ;
            gurobiModel.A = milpModel.A ;
            gurobiModel.rhs = milpModel.b ;
            gurobiModel.lb = - inf * ones(size(milpModel.c)) ;
            gurobiModel.sense = '>' ;
            vtypes(1:size(milpModel.A,2)) = 'C' ;
            vtypes(milpModel.binaryIdx) = 'C' ;% 'B' ; % Otherwise no dual variables
            vtypes(milpModel.integerIdx) = 'C' ;% 'I' ;
            gurobiModel.vtype = char(vtypes) ;
            gurobiParams.outputflag = 0 ;
            result = gurobi(gurobiModel, gurobiParams) ;
            if strcmpi(result.status,'OPTIMAL')
                diagnostics.solved = true ;
            else
                diagnostics.solved = false ;
                warning('Error while solving problem using gurobi') ;
                disp(result) ;
            end
            
            % We want the results :-)
            xVal = milpModel.primalValue(x, result.x) ;
            
            % We also want dual variables :-)
            cntr1Dual = milpModel.dualValue(cntr1, result.pi) ;
            cntr2Dual = milpModel.dualValue(cntr2, result.pi) ;
            
            testCase.verifyEqual(xVal, [1/3; 0; 0; 1/6], 'AbsTol', 1e-14) ;
            testCase.verifyEqual(cntr1Dual, [0;0;0;0], 'AbsTol', 1e-14) ;
            testCase.verifyEqual(cntr2Dual, [-1/6; 1/6], 'AbsTol', 1e-14) ;
        end
        
        % Check of duals (signs) + model
        function testResult2(testCase)
            % Model
            x = sddpVar(2,1) ;
            cntr1 = x >= 0 ;
            cntr2 = 2*x(1)+x(2) <= 2 ;
            cntr3 = x(1)+2*x(2) <= 2 ;
            obj = -sum(x) ;
            milpModel = MILPModel([cntr1, cntr2, cntr3], obj) ;
            
            % Check matrices
            Acheck = sparse([eye(2) ; -2 -1 ; -1 -2]) ;
            bcheck = [0 ; 0 ; -2 ; -2] ;
            ccheck = -[1 1]' ;
            kcheck = 0 ;
            testCase.verifyEqual(milpModel.A,Acheck) ;
            testCase.verifyEqual(milpModel.b,bcheck) ;
            testCase.verifyEqual(milpModel.c,ccheck) ;
            testCase.verifyEqual(milpModel.k,kcheck) ;
            
            % Run Gurobi
            gurobiModel.obj = milpModel.c ;
            gurobiModel.A = milpModel.A ;
            gurobiModel.rhs = milpModel.b ;
            gurobiModel.lb = - inf * ones(size(milpModel.c)) ;
            gurobiModel.sense = '>' ;
            vtypes(1:size(milpModel.A,2)) = 'C' ;
            vtypes(milpModel.binaryIdx) = 'C' ;% 'B' ; % Otherwise no dual variables
            vtypes(milpModel.integerIdx) = 'C' ;% 'I' ;
            gurobiModel.vtype = char(vtypes) ;
            gurobiParams.outputflag = 0 ;
            result = gurobi(gurobiModel, gurobiParams) ;
            if strcmpi(result.status,'OPTIMAL')
                diagnostics.solved = true ;
            else
                diagnostics.solved = false ;
                warning('Error while solving problem using gurobi') ;
                disp(result) ;
            end                        
            testCase.verifyEqual(result.x, [2/3;2/3], 'AbsTol',1e-12) ;            
            xVal = milpModel.primalValue(x, result.x) ;
            cntr1Dual = milpModel.dualValue(cntr1, result.pi) ;
            cntr2Dual = milpModel.dualValue(cntr2, result.pi) ;
            cntr3Dual = milpModel.dualValue(cntr3, result.pi) ;            
            testCase.verifyEqual(xVal, [2/3;2/3], 'AbsTol', 1e-12) ;
            testCase.verifyEqual(cntr1Dual, [0;0]) ;            
            testCase.verifyEqual(cntr2Dual, -1/3, 'AbsTol', 1e-12) ;
            testCase.verifyEqual(cntr3Dual, -1/3, 'AbsTol', 1e-12) ;
        end
        
        % No objective check, mix of ineq/eq, check of dummy ineqs
        function testResult3(testCase)
            x = sddpVar([2 1]) ;            
            dummy = sddpVar(4,5,3,2) ;
            y = sddpVar() ;
            cntr1 = [x == 0] ;
            cntr2 = [y <= x] ;
            dummy = [x+y == 3] ;
            cntr3 = [y >= x - 3] ;     
            cntr = [cntr1, cntr2, cntr3] ;
            milpModel = MILPModel(cntr) ;    
            
            Acheck = sparse([1 0  0 ; 0 1  0 ; -1 0 0 ; 0 -1 0 ; ...
                             1 0 -1 ; 0 1 -1 ; -1 0 1 ; 0 -1 1]) ;
            bcheck = [0 0 0 0 0 0 -3 -3]' ;
            ccheck = zeros(3,1) ;
            kcheck = 0 ;
            testCase.verifyEqual(milpModel.A, Acheck) ;
            testCase.verifyEqual(milpModel.b, bcheck) ;
            testCase.verifyEqual(milpModel.c, ccheck) ;
            testCase.verifyEqual(milpModel.k, kcheck) ;            
        end
        
    end
end