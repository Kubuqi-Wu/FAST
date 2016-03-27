classdef DualsTests < matlab.unittest.TestCase
    
    methods (Test)
        
        function testResults1(testCase)
                        
            lattice = Lattice.latticeEasy(2, 2, @demand) ; 
            params = sddpSettings('solver','gurobi',...
                'verbose',0) ;
            var.x = sddpVar(2) ;
            var.s = sddpVar(2) ;
            cntrTarget1 = var.s(1) <= var.x(1) ;
            cntrTarget2 = var.s(2) <= var.x(1) + var.x(2)  ; % node (t, n) = (2, 2)
            cntrTarget3 = var.s >= 0 ;
            cntrTarget4 = var.s <= [1 ; 2] ;
            lattice = compileLattice(lattice,@(scenario)nlds(scenario,var,cntrTarget1,cntrTarget2,cntrTarget3,cntrTarget4),params) ;  
            [~,~,~,sol] = forwardPass(lattice,[1 ; 2],params) ;
            
            warning('off','all') ;
            duals1 = lattice.getDualSolution(cntrTarget1, sol) ;            
            duals2 = lattice.getDualSolution(cntrTarget2, sol) ;            
            duals3 = lattice.getDualSolution(cntrTarget3, sol) ;            
            duals4 = lattice.getDualSolution(cntrTarget4, sol) ;            
            warning('on','all') ;
            
            testCase.verifyEqual(duals1,-1) ;
            testCase.verifyEqual(duals2,-2) ;
            testCase.verifyEqual(duals3,[0;0]) ;
            testCase.verifyEqual(duals4,[0;0]) ;
                          
            function out = demand(t,i)
                if t == 1
                    out = [] ;
                else
                    if i == 1
                        out = [1 ; 1] ;
                    elseif i == 2
                        out = [1 ; 2] ;
                    else
                        error('wrong index') ;
                    end
                end
            end            
            
            function [cntr, obj] = nlds(scenario, var, cntrTarget1, cntrTarget2, cntrTarget3, cntrTarget4)
                x = var.x ;
                s = var.s ;
                if(scenario.getTime() == 1)
                    cntr = x >= 0 ;
                    obj = x(1) + x(2) ;
                end
                if(scenario.getTime() == 2)
                    obj = - s(1) - 2*s(2) ;
                    if scenario.getIndex() == 1                        
                        cntr = [s(1) <= x(1) ; ...
                                s(2) <= x(1) + x(2) ; ...
                                s >= 0 ;
                                s <= scenario.data] ;
                    else                         
                        cntr = [cntrTarget1 ;
                                cntrTarget2 ;
                                cntrTarget3 ;
                                cntrTarget4] ;
                    end                    
                end
            end  
            
        end 
        
        
        
        
    end
end