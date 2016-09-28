classdef PurgeCutsTests < matlab.unittest.TestCase
    
    methods (Test)
        function testResults(testCase)
            for PCN = [15 36 75 125 500]
                testPCN(PCN);
            end   
            function testPCN(PCN)
                lattice = Lattice.latticeEasy(2, 2, @demand) ; 
                params = sddpSettings('stop.iterationMax',150,...
                                      'stop.iterationMin',150,...
                                      'solver','gurobi',...
                                      'algo.McCount',5,...
                                      'algo.purgeCuts',true,...
                                      'algo.purgeCutsNumber',PCN,...
                                      'algo.checkRedondance',false,...
                                      'verbose',0) ;
                var.x = sddpVar(1,1) ;
                var.s = sddpVar(1,1) ;
                warning('off','all') ;            
                lattice = compileLattice(lattice,@(scenario)nlds(scenario,var),params) ;  
                output = sddp(lattice,params) ;
                warning('on','all') ;
                lattice = output.lattice;                
                s1 = size(lattice.graph{1}{1}.model.cutCoeffs,1) ;
                s2 = size(lattice.graph{1}{1}.model.cutRHS,1) ;
                testCase.verifyEqual(s1,s2);
                testCase.verifyEqual(s1,PCN) ;  
                for n = 1:2
                    s1 = size(lattice.graph{2}{n}.model.cutCoeffs,1) ;
                    s2 = size(lattice.graph{2}{n}.model.cutRHS,1) ;
                    testCase.verifyEqual(s1,s2);
                    testCase.verifyEqual(s1,0) ;
                end
            end
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
    end
end