classdef PrecutTests < matlab.unittest.TestCase
    
    methods (Test)
        
        function testResults1(testCase)
                        
            lattice = Lattice.latticeEasy(2, 2, @demand) ; 
            params = sddpSettings('solver','gurobi',...
                'verbose',0) ;
            var.trash = sddpVar(3, 2) ;
            var.x = sddpVar() ;
            var.s = sddpVar() ;
            lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;  
            
            lattice = precut(lattice, [var.x, var.s], [1.5 nan], params) ;
            [cutsCoeffs, cutsRHS] = lattice.getCuts(1, 1, [var.x]) ;
            
            testCase.verifyEqual(cutsCoeffs, 1) ;
            testCase.verifyEqual(cutsRHS, -1) ;                                    
            testCase.verifyError(@()precut(lattice,[var.x], [1.5 nan], params),?MException) ;
            testCase.verifyError(@()precut(lattice,[sddpVar()], [1.5], params),?MException) ;
            
            function out = demand(t,i)
                if t == 1
                    out = [] ;
                else
                    out = (i == 1) * 1 + (i == 2) * 2 ;
                end
            end            
            
            function [cntr, obj] = nlds(scenario, var)
                costP 	= 1 ;
                costS 	= 2 ;
                x = var.x ;
                s = var.s ;
                if(scenario.getTime() == 1)
                    cntr = (x>=0);
                    obj = costP*x;
                end
                if(scenario.getTime() == 2)
                    obj = -costS*s;
                    notExceedStock = ( s <= x );
                    notExceedDemand = ( s <= scenario.data );
                    pos = ( s >= 0 );
                    cntr = [notExceedDemand ; notExceedStock ; pos];
                end
            end  
            
        end 
        
        function testResults2(testCase)
            
            params = sddpSettings('solver','gurobi',...
                'verbose',0) ;
            var.x = sddpVar() ;
            var.s = sddpVar() ;
            
            lattice = Lattice.latticeEasy(2, 5, @demand) ;             
            lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;              
            lattice = precut(lattice, [var.x, var.s], [1.5 nan], params) ;
            [cutsCoeffs, cutsRHS] = lattice.getCuts(1, 1, [var.x]) ;            
            testCase.verifyEqual(cutsCoeffs, 8/5) ;
            testCase.verifyEqual(cutsRHS, -2/5) ;                                    
            testCase.verifyError(@()precut(lattice,[var.x], [1.5 nan], params),?MException) ;
            testCase.verifyError(@()precut(lattice,[sddpVar()], [1.5], params),?MException) ;
            
            lattice2 = Lattice.latticeEasy(2, 5, @demand) ;             
            lattice2 = compileLattice(lattice2,@(scenario)nlds(scenario,var),params) ;              
            lattice2 = precut(lattice2, [var.x, var.s], [2.234 nan], params) ;
            [cutsCoeffs, cutsRHS] = lattice2.getCuts(1, 1, [var.x]) ;            
            testCase.verifyEqual(cutsCoeffs, 6/5, 'Reltol',1e-14) ;
            testCase.verifyEqual(cutsRHS, -6/5, 'Reltol',1e-14) ;                                    
            
            function out = demand(t,i)
                if t == 1
                    out = [] ;
                else
                    out = i ;
                end
            end            
            
            function [cntr, obj] = nlds(scenario, var)
                costP 	= 1 ;
                costS 	= 2 ;
                x = var.x ;
                s = var.s ;
                if(scenario.getTime() == 1)
                    cntr = (x>=0);
                    obj = costP*x;
                end
                if(scenario.getTime() == 2)
                    obj = -costS*s;
                    notExceedStock = ( s <= x );
                    notExceedDemand = ( s <= scenario.data );
                    pos = ( s >= 0 );
                    cntr = [notExceedDemand ; notExceedStock ; pos];
                end
            end  
            
        end 
        
        function testResults3(testCase)
            
            params = sddpSettings('solver','gurobi',...
                'verbose',0) ;
            var.x = sddpVar(2) ;
            var.s = sddpVar(2) ;
            
            lattice = Lattice.latticeEasy(3, 2, @demand) ;             
            lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;              
            lattice = precut(lattice, [var.x ; var.s], [1.2 ; 2.3 ; nan ; nan], params) ;
            
            [cutsCoeffs, cutsRHS] = lattice.getCuts(2, 1, [var.x(2) var.s(1)]) ;            
            testCase.verifyEqual(cutsCoeffs, [3/2 0]) ;
            testCase.verifyEqual(cutsRHS, -3) ;   
            [cutsCoeffs, cutsRHS] = lattice.getCuts(2, 2, [var.x(2) var.s(1)]) ;           
            testCase.verifyEqual(cutsCoeffs, [3/2 0]) ;
            testCase.verifyEqual(cutsRHS, -3) ;   
            
            [cutsCoeffs, cutsRHS] = lattice.getCuts(1, 1, var.x(1)) ;          
            testCase.verifyEqual(cutsCoeffs, [1]) ;
            testCase.verifyEqual(cutsRHS, -3) ;   
            
            
            function out = demand(t,i)
                if t == 1
                    out = [] ;
                elseif t == 2
                    out = i ;
                elseif t == 3
                    out = 2*i ;
                end
            end            
            
            function [cntr, obj] = nlds(scenario, var)                                
                x = var.x ;
                s = var.s ;
                if(scenario.getTime() == 1)
                    cntr = [x(1) >= 0] ;
                    obj  = x(1) ;
                elseif(scenario.getTime() == 2)
                    cntr = [x(2) >= 0 ; s(1) >= 0 ; s(1) <= x(1) ; s(1) <= scenario.data] ;
                    obj  = x(2) - 2*s(1) ;
                elseif(scenario.getTime() == 3)
                    cntr = [s(2) >= 0 ; s(2) <= x(2) ; s(2) <= scenario.data] ;
                    obj  = - 3*s(2) ;
                end
            end  
            
        end
        
        
    end
end