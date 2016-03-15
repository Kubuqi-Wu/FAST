classdef EVTest < matlab.unittest.TestCase
    
    methods (Test)
        function test1(testCase)  
            % min x - 2s
            % s.t.   x >= 0 
            %        s >= 0
            %        s <= x
            %        s <= 2.5
            % Sol is s = 2.5, x = 2.5, obj = -2.5
            lattice = Lattice.latticeEasy(2, 2, @optimize_stock_demand) ; 
            lattice = lattice.initExpectedLattice(@optimize_stock_demand);
            % Run SDDP
            params = sddpSettings('solver','gurobi') ;
            x = sddpVar(1,1) ;
            s = sddpVar(1,1) ;
            lattice = compileExpectedLattice(lattice,@optimize_stock_nlds,params) ;
            % ExpectedValue
            EV_val = expectedValue(lattice,params);
            testCase.verifyEqual(EV_val,-2.5) ;            
            function out = optimize_stock_demand(t,i)
                if t == 1
                    out = [] ;
                else
                    if(i == 1)
                        out = 2;
                    elseif (i == 2)
                        out = 3;
                    elseif (i == -1) % for expected lattice
                        out = 2.5;
                    end
                end
            end            
            function [cntr, obj] = optimize_stock_nlds(scenario)
                if(scenario.getTime() == 1)
                    cntr = x>=0;
                    obj = 1*x;
                end                
                if(scenario.getTime() == 2)
                    obj = -2*s;
                    cntr = [s <= scenario.data ; s <= x ; s >= 0];
                end                
            end            
        end
        
        % Check more complex case with constant part
        function test2(testCase)  
            % min x1 (+1) + 2x2 - 3s1 (+1.5) + 1s2 (-0.75) yep, it costs to sell s2
            % s.t.   x1,2 >= 0 
            %        s1,2 >= 0
            %        s1 <= x1
            %        s2 <= x2
            %        s1 <= 2.25
            %        s2 <= 1.5
            % Sol is s1 = 2.25, s2 = 0, x1 = 2.25, x2 = 0, obj = -4.5
            % (+1.75) = -2.75
            lattice = Lattice.latticeEasy(3, 2, @optimize_stock_demand) ; 
            lattice = lattice.initExpectedLattice(@optimize_stock_demand);
            % Run SDDP
            params = sddpSettings('solver','gurobi') ;
            x = sddpVar(2) ;
            s = sddpVar(2) ;
            lattice = compileExpectedLattice(lattice,@optimize_stock_nlds,params) ;
            % ExpectedValue
            EV_val = expectedValue(lattice,params);
            testCase.verifyEqual(EV_val,-2.75) ;            
            function out = optimize_stock_demand(t,i)
                if t == 1
                    out = [] ;
                elseif t == 2
                    if(i == 1)
                        out = 2;
                    elseif (i == 2)
                        out = 2.5;
                    elseif (i == -1) % for expected lattice
                        out = 2.25;
                    end
                elseif t == 3
                    if(i == 1)
                        out = 1;
                    elseif(i == 2)
                        out = 2 ;
                    elseif(i == -1) 
                        out = 1.5 ;
                    end
                end
            end            
            function [cntr, obj] = optimize_stock_nlds(scenario)
                if(scenario.getTime() == 1)
                    cntr = x(1)>=0;
                    obj = 1*x(1) + 1;
                end                
                if(scenario.getTime() == 2)
                    obj = -3*s(1) + 2*x(2) + 1.5;
                    cntr = [s(1) <= scenario.data ; s(1) <= x(1) ; s(1) >= 0 ; x(2) >= 0];
                end  
                if(scenario.getTime() == 3)
                    obj = s(2) - 0.75;
                    cntr = [s(2) <= scenario.data ; s(2) <= x(2) ; s(2) >= 0];
                end  
            end                
        end
        
    end
end