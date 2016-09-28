classdef LatticeTests < matlab.unittest.TestCase
    
    methods (Test)
        
        function testLatticeEasyMarkov(testCase) 
            H = 5;
            N = 4;
            P = [0.1 0.2 0.3 0.4 ;
                 0.2 0.4 0.3 0.1 ;
                 0.25 0.25 0.25 0.25 ;
                 0.9 0.05 0.04 0.01] ;
            lattice = Lattice.latticeEasyMarkov(H, P, 1) ;
            for i = 1:N
                testCase.verifyEqual(lattice.graph{1}{1}.transitionProba, P(1,:)');
            end
            for t = 2:H-1 
                for i = 1:N
                    testCase.verifyEqual(lattice.graph{t}{i}.transitionProba, P(i,:)');
                end
            end
            for i = 1:N
                testCase.verifyEqual(lattice.graph{H}{i}.transitionProba, []);
            end
        end
        
        function testLatticeEasyProbaConst(testCase)            
            lattice = Lattice.latticeEasyProbaConst(5, 3, [0.3;0.1;0.6], @(t,i) t+i)  ;
            testCase.verifyEqual(lattice.H,5) ;
            testCase.verifyEqual(size(lattice.graph),[5 1]) ;
            testCase.verifyEqual(lattice.expectedGraph,[]) ;
            testCase.verifyEqual(lattice.scenarioTable,[]) ;    
            for t = 1:5
                if t == 1
                    expectedSize = [1 1] ;
                    nNodes = 1 ;
                else
                    expectedSize = [3 1] ;
                    nNodes = 3 ;
                end
                if t < 5
                    transitionProba = [0.3;0.1;0.6] ;
                else
                    transitionProba = [] ;
                end
                testCase.verifyEqual(size(lattice.graph{t}),expectedSize) ;                                
                for n = 1:nNodes                    
                    testCase.verifyEqual(lattice.graph{t}{n}.index,n) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.time,t) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.transitionProba,transitionProba) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.data,t+n) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.model,[]) ;                    
                end                
            end            
        end
        
        function testLatticeEasyProbaConstError(testCase)            
            testCase.verifyError(@()Lattice.latticeEasyProbaConst(3, 3, [1/4;1/4;1/4;1/4]),?MException)                        
            testCase.verifyError(@()Lattice.latticeEasyProbaConst(2, 2, [1;1]),?MException)                        
        end
        
        function testLatticeEasy(testCase)
            lattice = Lattice.latticeEasy(5, 3, @(t,i) t+i) ;
            testCase.verifyEqual(lattice.H,5) ;
            testCase.verifyEqual(size(lattice.graph),[5 1]) ;
            testCase.verifyEqual(lattice.expectedGraph,[]) ;
            testCase.verifyEqual(lattice.scenarioTable,[]) ;    
            for t = 1:5
                if t == 1
                    expectedSize = [1 1] ;
                    nNodes = 1 ;
                else
                    expectedSize = [3 1] ;
                    nNodes = 3 ;
                end
                if t < 5
                    transitionProba = 1/3 * ones(3, 1) ;
                else
                    transitionProba = [] ;
                end
                testCase.verifyEqual(size(lattice.graph{t}),expectedSize) ;                                
                for n = 1:nNodes                    
                    testCase.verifyEqual(lattice.graph{t}{n}.index,n) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.time,t) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.transitionProba,transitionProba) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.data,t+n) ;
                    testCase.verifyEqual(lattice.graph{t}{n}.model,[]) ;                    
                end                
            end
        end    
        
        function testLatticeEasyErrors(testCase)
            %testCase.verifyError(@()Lattice.latticeEasy(1, 2),?MException) ;
            testCase.verifyError(@()Lattice.latticeEasy(2, 0),?MException) ;
            testCase.verifyError(@()Lattice.latticeEasy(-1, 0),?MException) ;
            testCase.verifyError(@()Lattice.latticeEasy(-1, 3),?MException) ;
            testCase.verifyError(@()Lattice.latticeEasy(2, -2),?MException) ;            
        end
    end
end