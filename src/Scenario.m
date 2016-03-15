classdef Scenario
    
    properties
        time ;  % 1 ... H
        index ;  % 1 ... nNode(t) or -1 for expected Scenario
        transitionProba ; % [double]
        data ; % typically struct()        
        model ; % NldsModel
    end
    
    methods
        
        function scenario = Scenario(time, index, transitionProba, data)
            scenario.time = time ;
            scenario.index = index ;
            scenario.transitionProba = transitionProba ;
            scenario.data = data ;
            scenario.model = [] ;
        end
        
        function time = getTime(scenario)
            time = scenario.time ;
        end
        
        function index = getIndex(scenario)
            index = scenario.index ;
        end
        
        function proba = getTransitionProba(scenario, index)
            proba = scenario.transitionProba(index) ;
        end
        
        function scenario = storeModel(scenario, model)
            if isempty(scenario.model)
                scenario.model = model ;
            else
                scenario.model = scenario.model.updateModel(model) ;
            end
        end
        
        function [solution, diagnostics] = solve(scenario, solutionPreviousTime, withoutTheta, params)
            [solution, diagnostics] = scenario.model.solve(solutionPreviousTime, withoutTheta, params) ;
        end
        
        
        function [solution, constraints, diagnostics] = solveUnderDecision(scenario, solutionPreviousTime, withoutTheta, decision, params)
            [solution, constraints, diagnostics] = scenario.model.solveUnderDecision(solutionPreviousTime, withoutTheta, decision, params) ;
        end
        
        function scenario = addCut(scenario, cutCoeffs, cutRHS, params)
            if isempty(scenario.model)
                error('Scenario.model should have been generated first.') ;
            end
            scenario.model = scenario.model.addCut(cutCoeffs, cutRHS, params) ;
        end
        
        function assignValueToVariables(scenario,x)
            scenario.model.assignValueToVariables(x) ;
        end
        
        function handle = plotCuts(scenario,idx,lb,ub,addCost)
            
            handle = figure ;
            hold on ;
            if numel(idx) == 1
                x = linspace(lb,ub,100) ;
                z = - inf * zeros(size(x)) ;
                for k = 1:size(scenario.model.cutCoeffs,1)
                    e = scenario.model.cutRHS(k) ;
                    E = scenario.model.cutCoeffs(k, idx) ;
                    plot(x, e - E * x, '--r') ;
                    z = max(z, e - E * x) ;                    
                end
                if addCost
                    z = z + scenario.model.c(idx) * x ;
                end
                plot(x, z, '-b') ;
            elseif numel(idx) == 2
                [x1,x2] = meshgrid(linspace(lb(1),ub(1),20),linspace(lb(2),ub(2),20)) ;
                z = - inf * zeros(size(x1)) ;
                for k = 1:size(scenario.model.cutCoeffs,1)
                    e = scenario.model.cutRHS(k) ;
                    E = scenario.model.cutCoeffs(k, idx) ;
                    for i = 1:size(x1,1) 
                        for j = 1:size(x1,2)
                            x = [x1(i,j) ; x2(i,j)] ;
                            z(i,j) = max(z(i,j), e - E * x) ;
                        end
                    end
                end
                if addCost
                    for i = 1:size(x1,1) 
                        for j = 1:size(x1,2)
                            x = [x1(i,j) ; x2(i,j)] ;
                            z(i,j) = z(i,j) + scenario.model.c(idx)' * x ;
                        end
                    end
                end
                surf(x1,x2,z) ;
                xlabel('x1') ;
                ylabel('x2') ;
                zlabel('z') ;
            else
                error('idx should contains one or two indices') ;
            end
        end
                
    end
    
end