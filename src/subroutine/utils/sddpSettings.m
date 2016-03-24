function params = sddpSettings(varargin)
%% SDDPSettings Constructs option structure for the SDDP algorithm
% 
% params = SDDPSETTINGS(name1,value1,name2,value2,...)
% params = SDDPSETTINGS(paramsOld,name1,value1,...)
%
% Takes (Name,Value) pairs of argument
% The first argument may be a struct. In this case, parameters are
% overwritten by the one passed in argument.
% Case doesn't matter.
% Basic checks are performed on the input.
%
% Arguments (default value) :
%   * runId (date-hour) : the prefix of the log file
%   * verbose (1) : the level of verbose. 0 = nothing, 1 = standard, 2 =
%   lot of information printed
%   * solver ('linprog') : the solver to be used : 'linprog', 'mosek', 'gurobi', 'glpk'
%
%   * log.useDiary (false) : true if you want to write the log file on disk, false otherwise
%   * log.logFile ('sddp.log') : the name of the log file
%   * log.saveTempResults (false) : true if you want to save temps results
%   in the folder log.resultsFolder
%   * log.resultsFolder ('results') : the result folder where temp
%   results are stored
%
%   * algo.McCount (5) : the number of monte-carlo trials to use at each
%   iteration
%   * algo.purgeCuts (false) : wether the algorithm should remove too old
%   cuts
%   * algo.purgeCutsNumber (500) : the number of cuts to keep. The algorithm will 
%   keep the last algo.purgeCutsNumber cuts while removing the older ones and the duplicates
%   * algo.minTheta (-1e3) : the lower bound on theta. Necessary to avoid
%   unboundedness at the beginning.
%   * algo.deterministic (false) : true if you want to have a deterministic
%   algorithm. In this case, the algorithm will explore ALL possible
%   scenarios. Use with caution, since this number of scenario can be quite
%   high. This will discard the value of algo.McCount.
%   * algo.checkRedondance (true) : check redondance when adding a cut, in
%   order to avoid duplicate. Usefull, but take some times (~10/15%).
%   * algo.redondanceTol (1e-6) : the relative tolerance, to decide when
%   cuts are equals or not
%
%   * stop.timeMin (0) : the minimum time the algorithm should run. SDPP will not
%   stop before timeMin
%   * stop.timeMax (inf) : the maximum time. SDDP will stop immediately 
%   after timeMax
%   * stop.iterationMin (0) : the minimum number of iteration
%   * stop.iterationMax (20) : the maximum number of iteration
%   Note that the algorithm never stop before both timeMin and iterationMin
%   are met, but stops immediately after either timeMax or iterationMax is
%   met.
%   * stop.stopWhen ('pereira') : the stopping criterion. Can be 'pereira',
%   'std', 'pereira and std' or 'never'.
%   * stop.pereiraCoef (2) : the coefficient in Pereira's stopping
%   criterion. 1 corresponds to a 68% confidence interval, 2 to a 95%, etc.
%   The smaller the more precise.
%   * stop.stdMcCoef (0) : the coefficient in the minimum std criterion.
%   * stop.regular (true) : if set to true, the stopping criterion is
%   check at each iteration using the algo.McCount trials
%   * stop.precise (false) : if set to true, the stopping critetion is check 
%   during the "precise" computations (see below).
%
%   * precise.compute (false) : wether the algorithm should do "big batch"
%   of forward pass every precise.iterationStep iteration to compute a 
%   precise mean cost
%   * precise.iterationStep (inf) : the number of iterations between two
%   precise computation
%   * precise.count (0) : the number of Monte-Carlo trials at each precise
%   computation.
%   * precise.computeEnd (false) : wether SDDP should run a "big batch" of
%   forward pass after termination, to get a precise final mean cost
%
%   Other parameters should not be modified by the user
%
%   See also: SDDP


% Inspired (but not copied :-)) by the code sdpsettings of Johan Lofberg from the Yalmip toolbox
% Argument may be string, double, boolean

%% NAMES
names = {
    % General   
    'runId'
    'verbose'
    'debug'
    'solver'
    
    % Version
    'version.id'
    'version.type' 
    'version.warning'
    'version.copyright'
    
    % Logging    
    'log.useDiary'
    'log.saveTempResults'
    'log.logFile'
    'log.resultsFolder'
    
    % Algorithm
    'algo.McCount'
    'algo.purgeCuts'
    'algo.purgeCutsNumber'
    'algo.minTheta'
    'algo.deterministic'
    'algo.checkRedondance'
    'algo.redondanceTol'
    
    % Termination
    'stop.timeMin'      % Min time in [s]
    'stop.timeMax'      % Max time in [s]
    'stop.iterationMin' % Min # of iterations
    'stop.iterationMax' % Max # of iterations
    'stop.stopWhen'     % Need to all be met in order to stop    
    'stop.pereiraCoef'  % Coef in Pereira stoping criterion
    'stop.stdMcCoef'    % Requires an std <= abs(lb) to stop
    'stop.regular'      % Stop if met in regular case 
    'stop.precise'      % Stop if met in precise case
    
    % Precise Computation
    'precise.compute'        % do it or not
    'precise.computeEnd'     % do it or not at the end (after convergence)
    'precise.iterationStep'  % every # iterations
    'precise.count'          % # of MC
    
} ;
lowerNames = lower(names) ;

%% TYPES CHECK
% General    
paramsType.runId = 'string' ;
paramsType.verbose = [0 2] ;
paramsType.debug = 'boolean' ;
paramsType.solver = {'gurobi','linprog','mosek','glpk'} ;

% Version
paramsType.version.id = 'string' ;
paramsType.version.type = 'string' ;
paramsType.version.warning = 'string' ;
paramsType.version.copyright = 'string' ;

% Logging
paramsType.log.useDiary = 'boolean' ;    
paramsType.log.saveTempResults = 'boolean' ;
paramsType.log.logFile = 'string' ;
paramsType.log.resultsFolder = 'string' ;

% Algorithm
paramsType.algo.McCount = [1 inf] ; 
paramsType.algo.purgeCuts = 'boolean' ; 
paramsType.algo.purgeCutsNumber = [1 inf] ; 
paramsType.algo.minTheta = [-inf inf] ;
paramsType.algo.deterministic = 'boolean' ;
paramsType.algo.checkRedondance = 'boolean' ;
paramsType.algo.redondanceTol = [0 inf] ;

% Termination
paramsType.stop.timeMin = [0 inf] ;
paramsType.stop.timeMax = [0 inf] ;
paramsType.stop.iterationMin = [0 inf] ;
paramsType.stop.iterationMax = [0 inf] ;
paramsType.stop.stopWhen = {'pereira','std','pereira and std','never'} ;
paramsType.stop.pereiraCoef = [0 inf] ;
paramsType.stop.stdMcCoef = [0 1] ; 
paramsType.stop.regular = 'boolean' ;
paramsType.stop.precise = 'boolean' ;

% Precise
paramsType.precise.compute = 'boolean' ;
paramsType.precise.computeEnd = 'boolean' ;
paramsType.precise.iterationStep = [0 inf] ;
paramsType.precise.count = [0 inf] ;

%% LOAD & BUILD
if nargin > 0 && isstruct(varargin{1})
    % TODO : Should add a check here
    params = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General    
    params.runId = getRunId() ;
    params.verbose = 1 ;
    params.debug = false ;
    params.solver = 'linprog' ;
    
    % Version
    params.version.id = '0.9.1b' ;
    params.version.type = 'beta' ; % dev, alpha, beta or release 
    params.version.warning = 'This is still a beta version. Use with caution :-)' ;
    params.version.copyright = 'FAST  Copyright (C) 2015-2016, Cambier Léopold and Scieur Damien.\nThis program comes with ABSOLUTELY NO WARRANTY; for details see the LICENSE file.\nThis is free software, and you are welcome to redistribute it under the conditions of the GPLv3 license. See the LICENSE file.' ;
    
    % Logging
    params.log.useDiary = false ;
    params.log.saveTempResults = false ;
    params.log.logFile = 'sddp.log' ;
    params.log.resultsFolder = 'results' ;
    
    % Algorithm
    params.algo.McCount = 5 ; 
    params.algo.purgeCuts = true ; 
    params.algo.purgeCutsNumber = 500 ; 
    params.algo.minTheta = -1e3 ;
    params.algo.deterministic = false ;
    params.algo.checkRedondance = true ;
    params.algo.redondanceTol = 1e-6 ;
    
    % Termination
    params.stop.timeMin = 0 ;
    params.stop.timeMax = inf ;    
    params.stop.iterationMin = 0 ;
    params.stop.iterationMax = 20 ;
    params.stop.stopWhen = 'pereira' ; % 'pereira', 'std', 'pereira and std', 'never'
    params.stop.pereiraCoef = 2 ;
    params.stop.stdMcCoef = 0 ; 
    params.stop.regular = true ;
    params.stop.precise = false ;

    % Precise
    params.precise.compute = false ; 
    params.precise.computeEnd = false ; 
    params.precise.iterationStep = inf ;
    params.precise.count = 0 ;
end

expectString = true ;
for idArg = paramstart:nargin    
    arg = varargin{idArg} ;
    % We're expecting a string
    if expectString
        if ~ ischar(arg)
            error('Argument name %d should be a string.',idArg) ;
        end
        % It is indeed a string !
        % What is is real name in params-name table ?
        idName = strmatch(lower(arg),lowerNames,'exact') ;         
        if isempty(idName)
            error('Argument %s does''nt match any parameters name.',arg) ;
        elseif numel(idName) > 1
            error('Argument %s match multiple parameters name.',arg) ;
        end
        expectString = false ;
    % We're expecting the value
    else
        % Check name type
        name = names{idName} ;
        if isa(arg, 'function_handle')
            argStr = func2str(arg) ;
        elseif ismatrix(arg) || ischar(arg)
            argStr = num2str(arg) ;
        else
            error('Argument type should be string, double or function handle') ;
        end
        type = eval(['paramsType.' name]) ;
        if strcmp(type,'string')
            if ~ ischar(arg)
                error('Argument %s (param %s) should be a string.',argStr, name) ;                
            end
        elseif strcmp(type,'boolean')
            if ~ (arg == 0 || arg == 1)
                error('Argument %s (param %s) should be a boolean (0 or 1, true or false)',argStr, name) ;
            end
        elseif strcmp(type,'handle')
            if ~ isa(arg, 'function_handle')               
                error('Argument %s (param %s) should be a handle function', argStr, name) ;
            end            
        elseif isnumeric(type)
            if any(size(type) ~= [1 2])                
                error('Internal error : Wrong size matrix type') ;
            end
            if ~ isnumeric(arg) || ~ isscalar(arg)
                error('Argument %s (params %s) should be a double scalar',argStr, name) ;
            end
            if arg < type(1) || arg > type(2)
                error('Argument %s (params %s) should lie in the interval [%e,%e]',argStr, name, type(1), type(2)) ;
            end
        elseif iscell(type)
            if isempty(strmatch(arg, type))
                % TODO : display accepted values
                error('Argument %s (params %s) is not a supported string value', argStr, name) ;
            end
        else
            error('Internal error : Non supported type.') ;
        end        
        eval(['params.' name ' = arg ;']); 
        expectString = true ;
    end
end

if ~ expectString
    error('(Name,Value) should appears in pairs, but it looks like the value corresponding to argument %s is missing.',arg) ;
end

% Particular checks
if params.algo.deterministic && (params.precise.compute || params.precise.computeEnd)
    error('algo.deterministic is incompatible with precise.compute or precise.computeEnd') ;
end
if params.stop.precise 
    if ~ params.precise.compute
        error('using stop.precise requires precise.compute') ;
    end
end


end


