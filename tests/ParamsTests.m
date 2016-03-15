classdef ParamsTests < matlab.unittest.TestCase
    
    methods (Test)
        function testParams1(testCase)
            params = sddpSettings('verbose',1,...
                'runId','lol',...
                'sTop.iTeratIonMin',10,...
                'algo.deterMinistic',true,...
                'log.logfile','lolbis',...
                'precise.compute',false) ;
            testCase.verifyEqual(params.verbose,1) ;
            testCase.verifyEqual(params.runId,'lol') ;
            testCase.verifyEqual(params.stop.iterationMin,10) ;
            testCase.verifyEqual(params.algo.deterministic,true) ;
            testCase.verifyEqual(params.log.logFile,'lolbis') ;
            testCase.verifyEqual(params.precise.compute,false) ;
        end      
        
        function testParams2(testCase)            
            testCase.verifyError(@()sddpSettings('verbose2',1),?MException) ;
        end
        
        function testParams3(testCase)            
            testCase.verifyError(@()sddpSettings('algo.deterministicLol',1),?MException) ;
        end
        
        function testParams4(testCase)            
            testCase.verifyError(@() sddpSettings('termination.compute',1),?MException) ;
        end
        
        function testParamsTypeError(testCase)
            testCase.verifyError(@()sddpSettings('verbose',-1),?MException) ;
            testCase.verifyError(@()sddpSettings('verbose',3),?MException) ;
            testCase.verifyError(@()sddpSettings('precise.compute','lol'),?MException) ;
            testCase.verifyError(@()sddpSettings('stop.iterationMax',-3),?MException) ;
            testCase.verifyError(@()sddpSettings('runId',true),?MException) ;
            testCase.verifyError(@()sddpSettings('log.logFile',@(x)x^2),?MException) ;
        end
    end
end