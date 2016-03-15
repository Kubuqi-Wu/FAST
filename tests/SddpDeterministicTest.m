classdef SddpDeterministicTest < matlab.unittest.TestCase
    
    methods (Test)
        function testResults(testCase)
            lattice = Lattice.latticeEasy(2, 2, @demand) ; 
            params = sddpSettings('stop.iterationMax',10,...
                'algo.deterministic',true,...
                'stop.pereiraCoef',0,...
                'solver','gurobi',...
                'verbose',0) ;
            var.x = sddpVar(1,1) ;
            var.s = sddpVar(1,1) ;
            warning('off','all') ;            
            lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;  
            output = sddp(lattice,params) ;
            warning('on','all') ;
            testCase.verifyEqual(output.meanCost(end),-2) ;
            testCase.verifyEqual(output.lowerBounds(end),-2) ;            
            function out = demand(t,i)
                if t == 1
                    out = [] ;
                else
                    out = (i == 1) * 2 + (i == 2) * 3 ;
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
            lattice = Lattice.latticeEasyProbaConst(2, 2, [1/4 3/4]', @demand) ;
            params = sddpSettings('stop.iterationMax',100,...
                'algo.deterministic',true,...
                'stop.pereiraCoef',0,...
                'solver','gurobi',...
                'verbose',0) ;
            var.x = sddpVar(1,1) ;
            var.s = sddpVar(1,1) ;
            warning('off','all') ;            
            lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;  
            output = sddp(lattice,params) ;
            warning('on','all') ;           
            testCase.verifyEqual(output.meanCost(end),-2.5) ;
            testCase.verifyEqual(output.lowerBounds(end),-2.5) ;            
            function out = demand(t,i)
                if t == 1
                    out = [] ;
                else
                    out = (i == 1) * 2 + (i == 2) * 3 ;
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
            % 3 stage, 2 nodes, proba 1/2
            % min p(i) x(i) s.t. x(i) + x(i-1) >= random(xi(i)), x(i) >= 0
            % Value function at stage 2 is 
            % V(x2) = 0.5*3*(max(10-x2,0) + max(1-x2,0))
            % and at stage 1 is
            % V(x1) = 0.5 * ( 2 * max(2-x,1) + 3/2 * ( max(10-max(2-x,1),0) + max(1-max(2-x,1),0) ) ) + 
            %         0.5 * ( 2 * max(7-x,1) + 3/2 * ( max(10-max(7-x,1),0) + max(1-max(7-x,1),0) ) )
            % Mini is at x1 = 5, and the cost is 20.75            
            lattice = Lattice.latticeEasy(3, 2, @demand) ;
            params = sddpSettings('stop.iterationMax',100,...
                'algo.deterministic',true,...
                'stop.pereiraCoef',0,...
                'solver','gurobi',...
                'verbose',0) ;    
            x = sddpVar(3,1) ; 
            warning('off','all') ;            
            lattice = compileLattice(lattice,@nlds,params) ;  
            output = sddp(lattice,params) ;            
            warning('on','all') ;           
            testCase.verifyEqual(output.meanCost(end),20.75) ;
            testCase.verifyEqual(output.lowerBounds(end),20.75) ;
            for s = 1:4
                testCase.verifyEqual(output.solution{1,s}.primal(1),5) ;
            end            
            % Getting solution at time 1
            nForward = 10 ;
            x1 = zeros(nForward, 1) ;
            lattice = output.lattice ;
            for i = 1:nForward                
                [~,~,~,solution] = forwardPass(lattice,'random',params) ;
                xx = lattice.getPrimalSolution(x, solution) ;
                x1(i) = xx(1) ;
            end            
            testCase.verifyEqual(mean(x1),5) ;
            testCase.verifyEqual(std(x1),0) ;            
            function out = demand(t,i)
                if t == 1
                    out = 5 ;
                elseif t == 2
                    out = (i == 1) * 2 + (i == 2) * 7 ;
                elseif t == 3
                    out = (i == 1) * 1 + (i == 2) * 10 ;
                end
            end                      
            function [cntr, obj] = nlds(scenario)
                t = scenario.time ;               
                cost = [1 2 3] ;
                obj = cost(t) * x(t) ;
                if t == 1
                    cntr = [x(t) >= 0 ; ...
                            x(t) >= scenario.data] ;
                else
                    cntr = [x(t) >= 0 ; ...
                            x(t) + x(t-1) >= scenario.data] ;
                end
            end            
        end                        
    end
end