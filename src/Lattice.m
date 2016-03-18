classdef Lattice
    
    properties
        H ; % int
        graph ; % cell(cell(Scenario))
        expectedGraph ; % cell(Scenario)
        scenarioTable ; % nScenario x H array
        scenarioTableProba ; % nScenario x 1 array
    end
    
    methods
        
        % Constructor
        function lattice = Lattice(H, graph)
            lattice.H = H ;
            lattice.graph = graph ;
            lattice.expectedGraph = [] ;
        end
        
        % Return H
        function H = getH(lattice)
            H = lattice.H ;
        end
        
        % Return the total number of **scenarios** (!= nodes, basically
        % 2^(H-1) if 2 nodes/stages)
        function S = getNScenarios(lattice)
            S = 1 ;
            for t = 1:lattice.H
                S = S * length(lattice.graph{t}) ;
            end
        end
        
        % If currentScenario == [], return a random scenario at time 1
        % Otherwise, return a random scenario a time t+1
        % If t == H, return []
        function nextScenario = nextRandomScenario(lattice,currentScenario)
            if isempty(currentScenario)
                idx = randi(length(lattice.graph{1})) ;
                nextScenario = lattice.graph{1}{idx} ;
            else
                time = currentScenario.getTime() ;
                if time == lattice.H
                    nextScenario = [] ;
                else
                    real = rand() ;
                    scenarioIdx = 0 ;
                    while real >= 0
                        scenarioIdx = scenarioIdx + 1 ;
                        real = real - currentScenario.transitionProba(scenarioIdx) ;
                    end
                    nextScenario = lattice.graph{time+1}{scenarioIdx} ;
                end
            end
        end
        
        % If currentScenario == [], return the nth scenario at time 1 (n
        % increasing from call to call)
        % Otherwise, return the next unseen scenario a time t+1
        % If t == H, return []
        function nextScenario = nextDeterministicScenario(lattice,currentScenario,Mc)
            if isempty(currentScenario)
                idx = lattice.scenarioTable(Mc,1) ;
                nextScenario = lattice.graph{1}{idx} ;
            else
                time = currentScenario.getTime() ;
                if time == lattice.H
                    nextScenario = [] ;
                else
                    idx = lattice.scenarioTable(Mc,time+1) ;
                    nextScenario = lattice.graph{time+1}{idx} ;
                end
            end
        end
        
        function lattice = generateScenarioTable(lattice)
            % Compute scenario table
            nScenario = zeros(lattice.H, 1) ;
            for t = 1:lattice.H
                nScenario(t) = length(lattice.graph{t}) ;
            end
            lattice.scenarioTable = fillMatrix(nScenario) ; % nScenario x H
            % Compute probabilities
            nScenarios = size(lattice.scenarioTable, 1) ;
            lattice.scenarioTableProba = zeros(nScenarios, 1) ;
            for nSce = 1:nScenarios
                lattice.scenarioTableProba(nSce) = 1 ;
                for t = 1:lattice.H-1
                    lattice.scenarioTableProba(nSce) = lattice.scenarioTableProba(nSce) * ...
                        lattice.graph{t}{lattice.scenarioTable(nSce,t)}.transitionProba(lattice.scenarioTable(nSce,t+1)) ;
                end
            end
        end
        
        % If currentScenario == [], return the first scenario at time t
        % Else, return the next scenario, or [] if at the end
        function nextScenario = explore(lattice,currentScenario,t)
            if isempty(currentScenario)
                scenarioIdx = 1 ;
            else
                scenarioIdx = currentScenario.index + 1 ;
            end
            if scenarioIdx > length(lattice.graph{t})
                nextScenario = [] ;
            else
                nextScenario = lattice.graph{t}{scenarioIdx} ;
            end
        end
        
        % Builds all the models from the @nlds(scenario) function
        % regular = true => standard : compileLattice
        % regular = false => expected : compileExpectedLattice
        function lattice = compileAllLattices(lattice, nlds, regular)
            % Breaks (with an error) if the models are structurally different from node to
            % node at the same stage. Basically, add zero-coef variable to
            % avoid that
            H = lattice.H ;
            
            % 1) Compile all models at each node of the lattice, and return
            % A, b, c, k, backMapVar and backMapCntr
            basicModels = cell(H,1) ;
            for t = 1:H
                nNodes = length(lattice.graph{t}) ;
                basicModels{t} = cell(nNodes,1) ;
                if regular
                    nNodes = length(lattice.graph{t}) ;
                    for n = 1:nNodes
                        % We use the model provided by the user
                        [cntr,obj] = nlds(lattice.graph{t}{n}) ;
                        % We get, from there, all the numerical data
                        basicModels{t}{n} = MILPModel(cntr, obj) ;
                    end
                else
                    % We use the model provided by the user
                    [cntr,obj] = nlds(lattice.expectedGraph{t}) ;
                    % We get, from there, all the numerical data
                    basicModels{t}{1} = MILPModel(cntr, obj) ;
                end
            end
            
            % 2) We need to check that all the nodes at a given stage have
            % the same x(t) and x(t-1)
            for t = 1:H
                if regular
                    nNodes = length(lattice.graph{t}) ;
                else
                    nNodes = 1 ;
                end
                backMapVarNode1 = basicModels{t}{1}.backMapVar ;
                for n = 2:nNodes
                    backMapVarNodeN = basicModels{t}{n}.backMapVar ;
                    if any(size(backMapVarNode1) ~= size(backMapVarNodeN)) || ...
                            any(backMapVarNode1 ~= backMapVarNodeN)
                        error('The variables x(t) and x(t-1) should all be the same at a given stage') ;
                    end
                end
            end
            
            % 3) Now, we need to detect the dependant and independant
            % variables from stage to stage
            depVarForw = cell(H,1) ;
            indepVarForw = cell(H,1) ;
            depVarBack = cell(H,1) ;
            indepVarBack = cell(H,1) ;
            for t = 1:H-1
                backMapVarT = basicModels{t}{1}.backMapVar ;
                backMapVarT1 = basicModels{t+1}{1}.backMapVar ;
                [indepVarForw{t},indepVarBack{t+1},...
                    depVarForw{t},depVarBack{t+1}] = ...
                    variableDependance(backMapVarT,backMapVarT1) ;
            end
            indepVarForw{H} = (1:numel(basicModels{H}{1}.backMapVar))' ;
            depVarForw{H} = [] ;
            indepVarBack{1} = (1:numel(basicModels{1}{1}.backMapVar))' ;
            depVarBack{1} = [] ;                        
            
            % 4) From here, we can build the NLDS model
            for t = 1:H
                
                if regular
                    nNodes = length(lattice.graph{t}) ;
                else
                    nNodes = 1 ;
                end
                
                for n = 1:nNodes
                    
                    % Some annoying mapping
                    indepVarB       = indepVarBack{t} ;  % x(t)
                    depVarB         = depVarBack{t} ;    % x(t-1)
                    depVarF         = depVarForw{t} ;
                    backMapVar      = basicModels{t}{n}.backMapVar ; % all variable (in model index)
                    indepBackMapVar = backMapVar(indepVarB) ; % all variable indep from stage t-1
                    forwardMapVar   = backMapVar(depVarF) ;
                    [~,forwardIdx]  = intersect(indepBackMapVar,forwardMapVar) ;
                    
                    % Extract matrices
                    basicModel = basicModels{t}{n} ;
                    modelVarIdx = basicModel.backMapVar(indepVarB) ;
                    modelCntrIdx = basicModel.backMapCntr ;
                    
                    W = basicModel.A(:,indepVarB) ;
                    h = basicModel.b ;
                    T = basicModel.A(:,depVarB) ;
                    c = basicModel.c(indepVarB) ;
                    k = basicModel.k ;
                    
                    if norm(basicModel.c(depVarB)) ~= 0
                        error('There cannot be dependant (i.e. x(t-1)) variable in the objective') ;
                    end
                    
                    nldsModel = NldsModel(W,h,T,c,k,forwardIdx,modelVarIdx,modelCntrIdx) ;
                    if regular
                        lattice.graph{t}{n} = lattice.graph{t}{n}.storeModel(nldsModel) ;
                    else
                        lattice.expectedGraph{t} = lattice.expectedGraph{t}.storeModel(nldsModel) ;
                    end
                    
                end
            end
            
            % 5) Check that variable appear only in 1 stage
            for t = 1:H
                % Can check only node 1 since variables are the same across
                % stage
                if regular
                    vart = lattice.graph{t}{1}.model.modelVarIdx ;
                else
                    vart = lattice.expectedGraph{t}.model.modelVarIdx ;
                end
                for s = [1:t-1 t+1:H]                    
                    if regular
                        vars = lattice.graph{s}{1}.model.modelVarIdx ;
                    else
                        vars = lattice.expectedGraph{s}.model.modelVarIdx ;
                    end
                    if numel(intersect(vart,vars)) ~= 0
                        error('fast:Lattice:compileAllLattices:variableDependance',...
                            'A variable from stage %d also appears in stage %d. This is not permitted.',t,s) ;
                    end
                end
            end
            
            % 6) We finally check that the W matrix is the same across the
            % nodes at a given stage
            if regular
                for t = 1:H
                    nNodes = length(lattice.graph{t}) ;
                    W1 = lattice.graph{t}{1}.model.W ;
                    for n = 2:nNodes
                        WN = lattice.graph{t}{n}.model.W ;
                        if any(size(W1) ~= size(WN)) || ...
                                any(W1(:) ~= WN(:))
                            error('The constraints involving x(t) (i.e. the W matrix) should be the same at a given stage') ;
                        end
                    end
                end
            end
        end
        
        function lattice = compileExpectedLattice(lattice, nlds, params)
            lattice = compileAllLattices(lattice, nlds, false) ;
        end
        
        
        function lattice = compileLattice(lattice, nlds, params)
            lattice = compileAllLattices(lattice, nlds, true) ;
        end
        
        function lattice = generateExpectedModels(lattice, ndls, params)
            % Need to be rewritten anyway
            error('Not implemented yet') ;
        end
        
        % Return all nodes at a given stage
        function scenariosCell = getScenariosCells(lattice, time)
            scenariosCell = lattice.graph{time} ;
        end
        
        % cutCoeffs = {[],[],[]}
        % cutRHS    = (x,x,x)
        function lattice = addCuts(lattice, scenario, cutCoeffs, cutRHS, params)
            time = scenario.time ;
            index = scenario.index ;
            cutCoeffs = cat(1, cutCoeffs{:}) ;
            cutRHS = reshape(cutRHS,[numel(cutRHS) 1]) ;
            lattice.graph{time}{index} = lattice.graph{time}{index}.addCut(cutCoeffs, cutRHS, params) ;
        end
        
        % Display models
        function displayModels(lattice)
            for t = 1:lattice.H
                for s = 1:length(lattice.graph{t})
                    fprintf('%d -- %d\n',t,s) ;
                    disp('min c''x(t)+k s.t. Wx(t) >= h - Tx(t-1)') ;
                    disp('c') ;
                    disp(lattice.graph{t}{s}.model.c)
                    disp('k') ;
                    disp(full(lattice.graph{t}{s}.model.k))
                    disp('W') ;
                    disp(full(lattice.graph{t}{s}.model.W))
                    disp('h') ;
                    disp(full(lattice.graph{t}{s}.model.h))
                    disp('T') ;
                    disp(full(lattice.graph{t}{s}.model.T))
                    disp('modelValIdx')
                    disp(full(lattice.graph{t}{s}.model.modelVarIdx))
                    disp('modelCntrIdx')
                    disp(full(lattice.graph{t}{s}.model.modelCntrIdx))
                end
            end
        end
        
        % Plot the lattice
        function lattice = plotLattice(lattice,dataToString)
            if nargin == 1
                dataToString = @(data) '' ;
            elseif isempty(dataToString)
                dataToString = @(data) '' ;
            end
            hold on ;
            nNodesMax = 0 ;
            for t = 1:lattice.H
                nNodes = length(lattice.graph{t}) ;
                nNodesMax = max(nNodesMax,nNodes) ;
                for i = 1:nNodes % Add node
                    currentT = t ;
                    currentI = i-(1+nNodes)/2  ;
                    plot(currentT,currentI,'.','MarkerSize',20) ;
                    text(currentT+0.03,currentI+0.03,[num2str(lattice.graph{t}{i}.index) ' ' dataToString(lattice.graph{t}{i}.data)]) ;
                    if t < lattice.H % Add links with weight = proba
                        nNodesNext = length(lattice.graph{t+1}) ;
                        for j = 1:nNodesNext
                            nextT = t+1 ;
                            nextI = j-(1+nNodesNext)/2 ;
                            proba = lattice.graph{t}{i}.transitionProba(j) ;
                            
                            x = (currentT+nextT)/2 ;
                            y = (currentI+nextI)/2 ;
                            
                            plot([currentT nextT],[currentI nextI],'LineWidth',2*proba) ;
                            
                            angle = atand(nextI - currentI) ;
                            
                            text(x,y+0.03,num2str(proba),'Rotation',angle) ;
                        end
                    end
                end
            end
            xlim([0 lattice.H+1]) ;
            ylim([-nNodesMax/2 nNodesMax/2]) ;
            hold off ;
        end
        
        
        function path = randomPath(lattice)
            sce = [];
            path = zeros(lattice.H,1);
            for t=1:lattice.H
                sce = lattice.nextRandomScenario(sce);
                path(t) = sce.index;
            end
        end
        
        function lattice = clearCuts(lattice)
            for i=1:lattice.H
                for j=1:length(lattice.graph{i})
                    lattice.graph{i}{j}.model = lattice.graph{i}{j}.model.clearCuts();
                end
            end
        end
        
        function primalSolution = getPrimalSolution(lattice, variable, solutionForward)
        %GETPRIMALSOLUTION Primal solution of a forward pass
        %
        % x = lattice.GETPRIMALSOLUTION(variable, solutionForward) (where
        % solutionForward is the output of forwardPass) returns the values
        % corresponding to variable in x.
        %
        % Example: 
        % x = sddpVar()
        % (...)
        % output = sddp(...)
        % lattice = output.lattice
        % for i = 1:5
        %   [~,~,~,solution] = forwardPass(lattice,lattice.randomPath(),params) ;    
        %   xVec(i) = lattice.getPrimalSolution(x, solution) ;        
        % end
        %
        % See also FORWARDPASS
            if numel(solutionForward) ~= lattice.H || ~ iscell(solutionForward)
                error('solutionForward should be a cell of H elements') ;
            end
            primalSolution = nan(size(variable)) ;
            for t = 1:lattice.H
                model = lattice.graph{t}{1}.model ;
                primalSolution = model.updatePrimalSolution(primalSolution, variable, solutionForward{t}) ;
            end
        end
        
        function [cutsCoeffs, cutsRHS] = getCuts(lattice, time, scenarioId, variables)
        %GETCUTSCOEFFS Get cut coeffs and RHS at a given node
        %
        % [coeffs, rhs] = lattice.GETCUTSCOEFFS(time, scenarioId,
        % variables) returns the coefficients of variables at node (time,
        % scenarioId).
        % variables should contains exactly the variables involved at node
        % (time, scenarioId), in any order. The order of the coefficients
        % in the output is the same as the order of the variables provided.
        % So variables is an array of sddpVar().
        %
        % Internally, the cuts are stored as
        % coeffs^t x + theta >= rhs
        %
        % See also GETPRIMALSOLUTION
            if ~ isvector(variables)
                error('sddp:Lattice:getCutsCoeffs','Variables should be a vector (row or column) of sddpVar()') ;
            end
            variables = variables(:) ;
            for i = 1:numel(variables)
                if ~ variables(i).isVariable()
                    error('sddp:Lattice:getCutsCoeffs','Variables should only contain independant variables') ;
                end
            end
            [cutsCoeffs, cutsRHS] = lattice.graph{time}{scenarioId}.model.getCutsCoeffs(variables) ;       
        end
        
        function lattice = initExpectedLattice(lattice, data)
            
            H = lattice.getH();
            expectIdx = -1;
            
            if H < 2
                error('H should be >= 2') ;
            end
            if nargin == 1
                data = @(t,i) [] ;
            elseif isempty(data)
                data = @(t,i) [] ;
            end
            lattice.expectedGraph = cell(H,1) ;
            for t = 1:H
                if t < H
                    lattice.expectedGraph{t} = Scenario(t, expectIdx, 1, data(t,expectIdx)) ;
                else
                    lattice.expectedGraph{t} = Scenario(t, expectIdx, [], data(t,expectIdx)) ;
                end
                
            end
            
        end
    end
    
    
    methods(Static)
        
        % Create simple lattice with 1 node at the top, nNode scenario at
        % each stage and transition probabilities described by transProb.
        % We assume serial independence.
        % transitionProba is a nNodes times H-1 matrix build on the following
        % way: transprob(i,j) = probability to go to node i of stage j+1
        % from stage j. transitionProba should be of size nNodes x H-1
        % Data is optional, or can be empty.
        % If non-empty, data should be a function of two variables, time and
        % index, returning some data to be stored in the node (time,index).
        function lattice = LatticeEasyProbaNonConst(H, nNodes, transitionProba, data)
            
            if H < 1
                error('H should be >= 1') ;
            end
            if nNodes < 1
                error('nNode should be >= 1') ;
            end
            if any(size(transitionProba) ~= [nNodes H-1])
                error('probas should be of size [nNodes H-1].') ;
            end
            % Each column of transitionProba should sum to 1
            if any(abs(sum(transitionProba, 1)-1) > 1e-12)
                error('values in each column of transitionProba should sum to 1 (tol 1e-12)') ;
            end
            
            if nargin == 3
                data = @(t,i) [] ;
            elseif isempty(data)
                data = @(t,i) [] ;
            end
            
            graph = cell(H,1) ;
            graph{1} = cell(1,1) ;
            graph{1}{1} = Scenario(1, 1, transitionProba(:,1), data(1,1)) ;
            for t = 2:H
                graph{t} = cell(nNodes,1) ;
                for i = 1:nNodes
                    if t < H
                        graph{t}{i} = Scenario(t, i, transitionProba(:,t), data(t,i)) ;
                    else
                        graph{t}{i} = Scenario(t, i, [], data(t,i)) ;
                    end
                end
            end
            lattice = Lattice(H, graph) ;
            
        end
        
        % Create a lattice with H stages, 1 node at the top, and where
        % transition probability from each node of stage t to the nNodes of
        % stage t+1 are given by transitionProba
        % Data is optional, or can be empty.
        % If non-empty, data should be a function of two variables, time and
        % index, returning some data to be stored in the node (time,index).
        function lattice = latticeEasyProbaConst(H, nNodes, transitionProba, data)
            
            if H < 1
                error('H should be >= 1') ;
            end
            if nNodes < 1
                error('nNode should be >= 1') ;
            end
            if any(size(transitionProba) ~= [nNodes 1])
                error('probas should be of size [nNodes 1]') ;
            end
            if abs(sum(transitionProba)-1) > 1e-12
                error('values in transitionProba don''t sum to 1 (tol 1e-12)') ;
            end
            
            if nargin == 3
                data = @(t,i) [] ;
            elseif isempty(data)
                data = @(t,i) [] ;
            end
            
            graph = cell(H,1) ;
            graph{1} = cell(1,1) ;
            graph{1}{1} = Scenario(1, 1, transitionProba, data(1,1)) ;
            for t = 2:H
                graph{t} = cell(nNodes,1) ;
                for i = 1:nNodes
                    if t < H
                        graph{t}{i} = Scenario(t, i, transitionProba, data(t,i)) ;
                    else
                        graph{t}{i} = Scenario(t, i, [], data(t,i)) ;
                    end
                end
            end
            lattice = Lattice(H, graph) ;
            
        end
        
        
        % Create simple lattice with 1 node at the top, nNode scenario at
        % each stage and equals transition probabilities.
        % Data is optional, or can be empty.
        % If non-empty, data should be a function of two variables, time and
        % index, returning some data to be stored in the node (time,index).
        function lattice = latticeEasy(H, nNodes, data)
            
            if H < 1
                error('H should be >= 1') ;
            end
            if nNodes < 1
                error('nNode should be >= 1') ;
            end
            
            if nargin == 2
                data = @(t,i) [] ;
            elseif isempty(data)
                data = @(t,i) [] ;
            end
            
            transitionProba = 1/nNodes * ones(nNodes,1) ;
            graph = cell(H,1) ;
            graph{1} = cell(1,1) ;
            graph{1}{1} = Scenario(1, 1, transitionProba, data(1,1)) ;
            for t = 2:H
                graph{t} = cell(nNodes,1) ;
                for i = 1:nNodes
                    if t < H
                        graph{t}{i} = Scenario(t, i, transitionProba, data(t,i)) ;
                    else
                        graph{t}{i} = Scenario(t, i, [], data(t,i)) ;
                    end
                end
            end
            lattice = Lattice(H, graph) ;
            
        end
        
    end
end