classdef MexTests < matlab.unittest.TestCase
    
    methods (Test)
        
        function testMergeCells(testCase)
            
            for ite = 1:100
                m = randi(20) ;
                n = randi(20) ;
                C = cell(m, n) ;
                for i = 1:m
                    for j = 1:n
                        C{i, j} = randi(100, [randi(10) randi(10)]) ;
                    end
                end
                merged = mergeCells(C) ;
                mergedCheck = [] ;
                for i = 1:numel(C)
                    data = C{i} ;
                    mergedCheck = [mergedCheck ; data(:)] ;
                end
                testCase.verifyEqual(any(merged ~= mergedCheck),false)
            end
            
        end
        
        function testPlusScalarFast(testCase)
            for ite = 1:100
                id1 = unique(randi(1000000, [10 1])) ;
                id2 = unique(randi(1000000, [8 1])) ;
                coefs1 = rand(size(id1)) ;
                coefs2 = rand(size(id2)) ;
                [ids, coefs] = plusScalarFast(id1, id2, coefs1, coefs2) ;
                [common, i1, i2] = intersect(id1, id2, 'stable') ;
                only1 = setdiff(1:numel(id1), i1) ;
                only2 = setdiff(1:numel(id2), i2) ;
                idsCheck = [common ; id1(only1) ; id2(only2)] ;
                coefsCheck = [coefs1(i1)+coefs2(i2) ; coefs1(only1) ; coefs2(only2)] ;
                testCase.verifyEqual((all(ids == idsCheck) && all(coefs == coefsCheck)),true) ;
            end
        end
        
        function testForwardMapping(testCase)
            for ite = 1:100
                backward = (randi(100, [150, 1])) ;
                ids = (randi(100, [30, 1])) ;
                forward = forwardMapping(backward, ids) ;
                forwardCheck = zeros(size(ids)) ;
                for i = 1:length(forwardCheck)
                    res = find(backward == ids(i)) ;
                    if ~isempty(res)
                        forwardCheck(i) = res(1) ;
                    end
                end
                testCase.verifyEqual(all(forward == forwardCheck),true) ;
            end
        end
    end
end