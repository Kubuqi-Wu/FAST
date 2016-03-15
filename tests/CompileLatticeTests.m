classdef CompileLatticeTests < matlab.unittest.TestCase
    
    methods (Test)
        
        % 2-stage test
        function test1(testCase)
            
            x = sddpVar(2,1) ;
            s = sddpVar(2,1) ;
            costProd = [2;3] ;
            sellPrice = [4;3] ;
            lattice = Lattice.latticeEasy(2, 2) ;
            
            params = sddpSettings() ;
            lattice = lattice.compileLattice(@nlds,params) ;
            
            testCase.verifyEqual(lattice.H,2) ;
            testCase.verifyEqual(size(lattice.graph),[2 1]) ;
            testCase.verifyEqual(size(lattice.graph{1}),[1 1]) ;
            testCase.verifyEqual(size(lattice.graph{2}),[2 1]) ;
            
            % Check (1,1)
            n11 = lattice.graph{1}{1} ;
            testCase.verifyEqual(n11.time,1) ;
            testCase.verifyEqual(n11.index,1) ;
            testCase.verifyEqual(n11.transitionProba,[0.5;0.5]) ;
            testCase.verifyEqual(n11.data,[]) ;
            model = n11.model ;
            testCase.verifyEqual(model.c,[2;3]) ;
            testCase.verifyEqual(model.W,speye(2)) ;
            testCase.verifyEqual(model.T,sparse(2,0)) ;
            testCase.verifyEqual(model.h,[0;0]) ;
            testCase.verifyEqual(model.k,0) ;
            testCase.verifyEqual(model.forwardIdx,[1;2]) ;
            testCase.verifyEqual(model.modelVarIdx,[1;2]) ;
            testCase.verifyEqual(model.modelCntrIdx,[1;1]) ;
            
            % Check (2,1)
            n21 = lattice.graph{2}{1} ;
            testCase.verifyEqual(n21.time,2) ;
            testCase.verifyEqual(n21.index,1) ;
            testCase.verifyEqual(n21.transitionProba,[]) ;
            testCase.verifyEqual(n21.data,[]) ;
            model = n21.model ;
            testCase.verifyEqual(model.c,-[4;3]) ;
            testCase.verifyEqual(model.W,[speye(2);-speye(2);-speye(2)]) ;
            testCase.verifyEqual(model.T,[sparse(2,2);speye(2,2);sparse(2,2)]) ;
            testCase.verifyEqual(model.h,[zeros(4,1);-[2;3]]) ;
            testCase.verifyEqual(model.k,0) ;
            testCase.verifyEqual(model.forwardIdx,double.empty(0,1)) ;
            testCase.verifyEqual(model.modelVarIdx,[3;4]) ;
            testCase.verifyEqual(model.modelCntrIdx,[2;2;3;3;4;4]) ;
            
            % Check (2,2)
            n22 = lattice.graph{2}{2} ;
            testCase.verifyEqual(n22.time,2) ;
            testCase.verifyEqual(n22.index,2) ;
            testCase.verifyEqual(n22.transitionProba,[]) ;
            testCase.verifyEqual(n22.data,[]) ;
            model = n22.model ;
            testCase.verifyEqual(model.c,-[4;3]) ;
            testCase.verifyEqual(model.W,[speye(2);-speye(2);-speye(2)]) ;
            testCase.verifyEqual(model.T,[sparse(2,2);speye(2,2);sparse(2,2)]) ;
            testCase.verifyEqual(model.h,[zeros(4,1);-[2;4]]) ;
            testCase.verifyEqual(model.k,0) ;
            testCase.verifyEqual(model.forwardIdx,double.empty(0,1)) ;
            testCase.verifyEqual(model.modelVarIdx,[3;4]) ;
            testCase.verifyEqual(model.modelCntrIdx,[5;5;6;6;7;7]) ;
            
            function [cntr, obj] = nlds(scenario)
                t = scenario.getTime() ;
                i = scenario.getIndex() ;
                if t == 1
                    cntr = [x >= 0] ;
                    obj = costProd' * x ;
                elseif t == 2
                    if i == 1
                        demand = [2;3] ;
                    elseif i == 2
                        demand = [2;4] ;
                    else
                        error('Wrong index') ;
                    end
                    cntr = [s >= 0 ;
                        s <= x ;
                        s <= demand] ;
                    obj = - sellPrice' * s ;
                else
                    error('Wrong time') ;
                end
            end
        end
        
        % Weird case
        function test2(testCase)
            % Strange model
            % (t=1) ----- (t=2) ---- (t=3)
            %               *          *
            %   *
            %               *          *
            %
            % b1,b2       a1-3        c1
            
            a = sddpVar(3,1) ;
            b = sddpVar(2,1) ;
            c = sddpVar(1,1) ;
            lattice = Lattice.latticeEasy(3, 2) ;
            
            params = sddpSettings() ;
            lattice = lattice.compileLattice(@nlds,params) ;
            
            testCase.verifyEqual(lattice.H,3) ;
            testCase.verifyEqual(size(lattice.graph),[3 1]) ;
            testCase.verifyEqual(size(lattice.graph{1}),[1 1]) ;
            testCase.verifyEqual(size(lattice.graph{2}),[2 1]) ;
            testCase.verifyEqual(size(lattice.graph{3}),[2 1]) ;
            
            % Check (1,1)
            n11 = lattice.graph{1}{1} ;
            testCase.verifyEqual(n11.time,1) ;
            testCase.verifyEqual(n11.index,1) ;
            testCase.verifyEqual(n11.transitionProba,[0.5;0.5]) ;
            testCase.verifyEqual(n11.data,[]) ;
            model = n11.model ;
            testCase.verifyEqual(model.c,[1;2]) ;
            testCase.verifyEqual(model.W,speye(2)) ;
            testCase.verifyEqual(model.T,sparse(2,0)) ;
            testCase.verifyEqual(model.h,[0;0]) ;
            testCase.verifyEqual(model.k,3) ;
            testCase.verifyEqual(model.forwardIdx,[1;2]) ;
            testCase.verifyEqual(model.modelVarIdx,[8;9]) ; % see test1
            testCase.verifyEqual(model.modelCntrIdx,[8;8]) ; % see test1
            
            % Check (2,1)
            n21 = lattice.graph{2}{1} ;
            testCase.verifyEqual(n21.time,2) ;
            testCase.verifyEqual(n21.index,1) ;
            testCase.verifyEqual(n21.transitionProba,[0.5;0.5]) ;
            testCase.verifyEqual(n21.data,[]) ;
            model = n21.model ;
            testCase.verifyEqual(model.c,[1;1;1]) ;
            testCase.verifyEqual(model.W,[speye(3);-speye(3);speye(3)]) ;
            testCase.verifyEqual(model.T,[sparse(3,2);ones(3,1) sparse(3,1);sparse(3,2)]) ;
            testCase.verifyEqual(model.h,[0;0;0;-1;-1;-1;0;0;0]) ;
            testCase.verifyEqual(model.k,0) ;
            testCase.verifyEqual(model.forwardIdx,[1;2;3]) ;
            testCase.verifyEqual(model.modelVarIdx,[5;6;7]) ; % see test1
            testCase.verifyEqual(model.modelCntrIdx,[9;9;9;10;10;10;11;11;11]) ; % see test1
            
            % Check (2,2)
            n22 = lattice.graph{2}{2} ;
            testCase.verifyEqual(n22.time,2) ;
            testCase.verifyEqual(n22.index,2) ;
            testCase.verifyEqual(n22.transitionProba,[0.5;0.5]) ;
            testCase.verifyEqual(n22.data,[]) ;
            model = n22.model ;
            testCase.verifyEqual(model.c,[1;1;1]) ;
            testCase.verifyEqual(model.W,[speye(3);-speye(3);speye(3)]) ;
            testCase.verifyEqual(model.T,[sparse(3,2);sparse(3,1) ones(3,1);sparse(3,2)]) ;
            testCase.verifyEqual(model.h,[0;0;0;2;2;2;0;0;0]) ;
            testCase.verifyEqual(model.k,0) ;
            testCase.verifyEqual(model.forwardIdx,[1;2;3]) ;
            testCase.verifyEqual(model.modelVarIdx,[5;6;7]) ; % see test1
            testCase.verifyEqual(model.modelCntrIdx,[12;12;12;13;13;13;14;14;14]) ; % see test1
            
            % Check (3,1)
            n31 = lattice.graph{3}{1} ;
            testCase.verifyEqual(n31.time,3) ;
            testCase.verifyEqual(n31.index,1) ;
            testCase.verifyEqual(n31.transitionProba,[]) ;
            testCase.verifyEqual(n31.data,[]) ;
            model = n31.model ;
            testCase.verifyEqual(model.c,0) ;
            testCase.verifyEqual(model.W,sparse([1;-1])) ;
            testCase.verifyEqual(model.T,sparse([0 0 0 ; 1 1 1])) ;
            testCase.verifyEqual(model.h,[0;0]) ;
            testCase.verifyEqual(model.k,2) ;
            testCase.verifyEqual(model.forwardIdx,double.empty(0,1)) ;
            testCase.verifyEqual(model.modelVarIdx,[10]) ; % see test1
            testCase.verifyEqual(model.modelCntrIdx,[15;16]) ; % see test1
            
            % Check (3,2)
            n32 = lattice.graph{3}{2} ;
            testCase.verifyEqual(n32.time,3) ;
            testCase.verifyEqual(n32.index,2) ;
            testCase.verifyEqual(n32.transitionProba,[]) ;
            testCase.verifyEqual(n32.data,[]) ;
            model = n32.model ;
            testCase.verifyEqual(model.c,0) ;
            testCase.verifyEqual(model.W,sparse([1;-1])) ;
            testCase.verifyEqual(model.T,sparse([0 0 0 ; -1 -1 -1])) ;
            testCase.verifyEqual(model.h,[0;-1]) ;
            testCase.verifyEqual(model.k,2) ;
            testCase.verifyEqual(model.forwardIdx,double.empty(0,1)) ;
            testCase.verifyEqual(model.modelVarIdx,[10]) ; % see test1
            testCase.verifyEqual(model.modelCntrIdx,[17;18]) ; % see test1
            
            
            function [cntr, obj] = nlds(scenario)
                t = scenario.getTime() ;
                i = scenario.getIndex() ;
                if t == 1
                    cntr = b >= 0 ;
                    obj = b(1)+2*b(2)+3 ;
                elseif t == 2
                    if i == 1
                        cntr = [a >= 0 ;
                            a <= b(1)+1 ;
                            a >= 0*b(2)] ;
                    elseif i == 2
                        cntr = [a >= 0 ;
                            a <= b(2)-2 ;
                            a >= 0*b(1)] ;
                    else
                        error('Wrong index') ;
                    end
                    obj = sum(a) ;
                elseif t == 3
                    if i == 1
                        cntr = [c >= 0 ;
                            c <= sum(a)] ;
                    elseif i == 2
                        cntr = [c >= 0 ;
                            c <= - sum(a) + 1] ;
                    else
                        error('Wrong index') ;
                    end
                    obj = sddpConst(2) ;
                else
                    error('Wrong time') ;
                end
            end
        end
        
        % forwardIdx Test
        function test3(testCase)
            
            x = sddpVar(10,1) ;
            lattice = Lattice.latticeEasy(3, 2) ;
            
            params = sddpSettings() ;
            lattice = lattice.compileLattice(@nlds,params) ;
            
            testCase.verifyEqual(lattice.H,3) ;
            testCase.verifyEqual(size(lattice.graph),[3 1]) ;
            testCase.verifyEqual(size(lattice.graph{1}),[1 1]) ;
            testCase.verifyEqual(size(lattice.graph{2}),[2 1]) ;
            
            model = lattice.graph{1}{1}.model ;
            testCase.verifyEqual(model.forwardIdx,[2]) ;
            model = lattice.graph{2}{1}.model ;
            testCase.verifyEqual(model.forwardIdx,[1;3]) ;
            model = lattice.graph{2}{2}.model ;
            testCase.verifyEqual(model.forwardIdx,[1;3]) ;
            model = lattice.graph{3}{1}.model ;
            testCase.verifyEqual(model.forwardIdx,double.empty(0,1)) ;
            model = lattice.graph{3}{2}.model ;
            testCase.verifyEqual(model.forwardIdx,double.empty(0,1)) ;            
            
            function [cntr, obj] = nlds(scenario)
                t = scenario.getTime() ;
                i = scenario.getIndex() ;
                if t == 1
                    cntr = x([1 4 8]) >= 0 ;
                    obj = sum(x([1 4 8])) ;
                elseif t == 2
                    if i == 1
                        cntr = sum(x([2 9 10])) <= x(4) ;
                    elseif i == 2
                        cntr = sum(x([2 9 10])) <= 2*x(4) ;
                    end
                    obj = sddpConst(1) ;
               elseif t == 3
                    if i == 1
                        cntr = sum(x([3 5 6 7])) <= x(2)+x(10) ;
                    elseif i == 2
                        cntr = sum(x([3 5 6 7])) <= 2*x(2)+3*x(10) ;
                    end
                    obj = sddpConst(1) ;
                end
            end
        end
        
        % Check error generation for badly formatted NLDS
        function testError(testCase)
            
            x = sddpVar() ;
            y = sddpVar() ;
            z = sddpVar() ;            
            lattice = Lattice.latticeEasy(2, 2) ;            
            params = sddpSettings() ;            
            testCase.verifyError(@()lattice.compileLattice(@nlds1,params),?MException) ;
            testCase.verifyError(@()lattice.compileLattice(@nlds2,params),?MException) ;
            testCase.verifyError(@()lattice.compileLattice(@nlds3,params),?MException) ;
                                   
            % Should crash because of unboundedness
            function [cntr, obj] = nlds1(scenario)
                t = scenario.getTime() ;
                i = scenario.getIndex() ;
                if t == 1
                    cntr = x >= 0 ;
                    obj = x ;
                elseif t == 2
                    if i == 1
                        cntr = y <= x ;
                    elseif i == 2
                        cntr = z <= x ;
                    end                    
                    obj = y+z ;
                end
            end
            % Should crash because variable not present at every node
            function [cntr, obj] = nlds2(scenario)
                t = scenario.getTime() ;
                i = scenario.getIndex() ;
                if t == 1
                    cntr = x >= 0 ;
                    obj = x ;
                elseif t == 2
                    if i == 1
                        cntr = y <= x ;
                        obj = y ;
                    elseif i == 2
                        cntr = z <= x ;
                        obj = z ;
                    end                                        
                end
            end
            % Should crash because not the same W
            function [cntr, obj] = nlds3(scenario)
                t = scenario.getTime() ;
                i = scenario.getIndex() ;
                if t == 1
                    cntr = x >= 0 ;
                    obj = x ;
                elseif t == 2
                    if i == 1
                        cntr = [y <= x ; ...
                                2*z <= x] ;
                    elseif i == 2
                        cntr = [y <= x ; ...
                                3*z <= x] ;
                    end                    
                    obj = y+z ;
                end
            end
            
        end
        
        function testError2(testCase)
        
            x = sddpVar() ;
            y = sddpVar() ;
            z = sddpVar() ;            
            lattice = Lattice.latticeEasy(3, 1) ;            
            params = sddpSettings() ;       
            lattice.compileLattice(@nlds,params) ;
            %testCase.verifyError(@()lattice.compileLattice(@nlds,params),?MException) ;
            
            function [cntr,obj] = nlds(scenario)
                t = scenario.getTime() ;
                if t == 1
                    cntr = [x >= 0] ;
                    obj = x ;
                elseif t == 2
                    cntr = [y <= x ; y >= 0] ;
                    obj = 2*y ;
                elseif t == 3
                    cntr = [z <= x ; z >= 0] ;
                    obj = 3*z ;
                end
            end
        end
        
        
        
    end
end