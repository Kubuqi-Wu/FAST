function interactive_demo_optimize_stock()
warning('Not updated. Will crash.') ;

clc;
clear all;
close all;

idx = 1;
tot = 14;

disp(' ')
disp(' ')
disp('--------------------- DEMO1 : Optimizing stock  ---------------------')
fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('This demo shows how to modelize a simple problem in order to use the ')
disp('SDDP algorithm.')
disp(' ')
disp(' ')
disp('If you want to stop the demo before the end type CTRL-C')
disp(' ')
disp(' ')
disp('Press any key to start...')
disp(' ')
disp(' ')
pause
idx=idx+1;



fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('Assume we want to solve the following problem :')
disp(' ')
disp('We work on a two-step timeline.')
disp('  - At time 1, we can decide to produce some amount of a product.')
disp('  - At time 2, we can sell this product, up to some demand.')
disp('But the problem is that demand at time 2 is unknown and random.')
disp(' ')
disp('The goal is to minimize the expectation of the costs.')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;



fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('Let''s define "C" the cost of production, "P" the sale price.')
disp('Let''s denote by "x" the amount of product we produce at time 1, and')
disp('by "s" the amount of product we sell at time 2.')
disp('Finally, let "d" be the (random) demand at time 2, function of ')
disp('the scenario xi.')
disp(' ')
disp('Formally, the problem at time 1 can be stated as')
disp('          min C x + E(V(x))')
disp('    s.t.  x >= 0')
disp('where E(V(x)) is the expectation of the cost of future stages.')
disp(' ')
disp('At stage 2, the problem is')
disp('          min - P s')
disp('    s.t.  s >= 0')
disp('          s <= d(xi)')
disp('          s <= x')
disp(' ')
disp('We will now see how to solve this problem using FAST.')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;




fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('   Step (1) : Build the lattice')
disp('   ------------------------------')
disp(' ')
disp('The lattice must describe all possible future scenarios. For example,')
disp('in our case we have a two-stage problem. Then, the horizon "H" of the')
disp('lattice will be equal to H=2.')
disp(' ')
disp('Also, we will assume that we can have two scenarios : low or high')
disp('demand (say d = 2 or d = 3) with probability 1/2 each. Then, the ')
disp('lattice will have only two nodes by stages.')
disp(' ')
disp('There exist a simple constructor for the class lattice : ')
disp('>>  Lattice.latticeEasy(H,nNodes,@scenarioParameters)')
disp(' ')
disp('In our case, H=2, nNodes=2 and the last function handle must')
disp('describe our demand in function of (t,xi).')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;





fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp(' ')
disp(' We will write a function demand(t,xi) :')
disp(' ')
disp('>>    function out = demand(t,xi)')
disp('>>        if t == 1  % We do not have demand at stage 1')
disp('>>            out = [] ; ')
disp('>>        else % At time 2, the demand is random, either 2 or 3.')
disp('>>            out = (xi == 1) * 2 + (xi == 2) * 3 ;')
disp('>>        end')
disp('>>    end')
disp(' ')
disp('And call the constructor :')
disp('>> H=2 ; nNodes=2')
disp('>> lattice = Lattice.latticeEasy(H, nNodes, @demand) ;')
disp(' ')
disp('Press any key to launch this call...')
disp(' ')
disp(' ')
pause
idx=idx+1;

    function out = demand(t,xi)
        if t == 1 % We don't have any randomnes at stage 1, no we don't store anything
            out = [] ;
        else % At time 2, the demand is random, either 2 or 3.
            out = (xi == 1) * 2 + (xi == 2) * 3 ;
        end
    end
H = 2; 
nNodes = 2 ;
lattice = Lattice.latticeEasy(H, nNodes, @demand) ;




fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('Let''s see what the lattice looks like.')
disp('The following 2 lines will display the lattice, as well as transition')
disp('probabilities (1/2 by default) and the data stored, using the @(data) num2str(data)')
disp('to convert the output of the demand function into string. This is')
disp('optional, and the following two lines will work too')
disp(' ')
disp('>> figure ;')
disp('>> lattice.plotLattice(@(data) num2str(data)) ;')
disp(' ')
disp('Press any key to launch this call...')
disp(' ')
disp(' ')
pause
idx=idx+1;
figure ;
lattice.plotLattice(@(data) num2str(data)) ;




fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('   Step (2) : Variables declaration')
disp('   ------------------------------')
disp(' ')
disp('Now, we need to define the variables we will use in our problem.')
disp('Here, we only need x (of size [1 x 1]) and s (size [1 x 1]).')
disp('So we build a function returning a structure with these two variables.')
disp('Note the ''full'' argument : it is a good habit to always put that ')
disp('(see Yalmip documentation for more info).')
disp('>>    function var = declareVariables()')
disp('>>        var.x = sdpvar(1,1,''full'') ;')
disp('>>        var.s = sdpvar(1,1,''full'') ;')
disp('>>    end')
disp(' ')
disp('Press any key to launch this call...')
disp(' ')
disp(' ')
pause
idx=idx+1;


    function var = declareVariables()
        var.x = sdpvar(1,1,'full') ;
        var.s = sdpvar(1,1,'full') ;
    end




fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('   Step (3) : NDLS declaration')
disp('   ------------------------------')
disp(' ')
disp('Finally, we need to create the function that, at each node,') 
disp('return the NDLS problem')
disp(' ')
disp('We have to write a function like:')
disp('>> function [constraints, objective] = ndls(scenario, var)')
disp(' ')
disp('Where var are the variables declared before (in our example, we ')
disp('have var.x and vax.s) and scenario contains information about where')
disp('we are in the lattice, like scenario.time (the stage), scenario.idx ')
disp('(the index of the node) and in scenario.data (in our case, equal to ')
disp('the demand d(scenario.idx)).')
disp(' ')
disp('The function must return Yalmip constraints and Yalmip objective for')
disp('the node (scenario.time,scenario.idx).')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;



fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('We will now see what we must put in this function.')
disp('First, we get all required data:')
disp('>> C = 1;')
disp('>> S = 2;')
disp('>> x = var.x;')
disp('>> s = var.s;')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;





fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('We need to describe the NDLS at stage 1 (when we produce the goods)')
disp(' ')
disp('At stage 1, the problem is')
disp('           min C x + E(V(x))')
disp('     s.t.  x >= 0')
disp(' So we write the following constraints :')
disp('           x >= 0')
disp(' and the objective function, without the E(V(x)) : the toolbox')
disp(' automatically takes care of that')
disp('           min C * x')
disp('Then we must write :')
disp('>>        if(scenario.time == 1)')
disp('>>            % Constraint')
disp('>>            cntr = x >= 0 ;')
disp('>>            % Objective')
disp('>>            obj  = C * x ;')
disp('>>        end')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;


fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('Now we are at stage 2 (we sell the goods).')
disp(' ')
disp('At stage 2, the problem is')
disp('           min - P s')
disp('    s.t.  s >= 0')
disp('          s <= d(xi)')
disp('          s <= x')
disp('(Remind that we minimize costs, so the ''-'' is needed)')
disp('>>        if(scenario.time == 2)')
disp('>>            % Constraints')
disp('>>            % data is the random demand store into each node (scenario)')
disp('>>            notExceedDemand = s <= scenario.data ; ')
disp('>>            notExceedStock  = s <= x ;')
disp('>>            pos             = s >= 0 ;')
disp('>>            % Concatenate all constraints together')
disp('>>            cntr = [notExceedDemand ; notExceedStock ; pos];')
disp('>>            % The objective')
disp('>>            obj = -S * s;')
disp('>>        end')
disp(' ')
disp('Press any key to continue...')
disp(' ')
disp(' ')
pause
idx=idx+1;





fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('Here is the complete function.')
disp(' ')
disp('>>    function [cntr, obj] = ndls(scenario, var)')
disp('>>        C 	= 1 ;')
disp('>>        S 	= 2 ;')
disp('>>        x   = var.x ;')
disp('>>        s   = var.s ;')
disp('>>        if(scenario.time == 1)')
disp('>>            cntr = x >= 0 ;')
disp('>>            obj  = C * x ;')
disp('>>        end')
disp('>>        if(scenario.time == 2)')
disp('>>            pos             = s >= 0 ;')
disp('>>            notExceedDemand = s <= scenario.data ; ')
disp('>>            notExceedStock  = s <= x ;')
disp('>>            cntr = [notExceedDemand ; notExceedStock ; pos];')
disp('>>            obj = -S * s;')
disp('>>        end')
disp('>>    end')
disp(' ')
disp('Press any key to launch this call...')
pause
idx=idx+1;

    function [cntr, obj] = ndls(scenario, var)
        C 	= 1 ;
        S 	= 2 ;
        x   = var.x ;
        s   = var.s ;
        
        % It at stage 1, the problem is
        %           min C x + E(V(x))
        %           s.t.  x >= 0
        % So we write the following constraints :
        %           x >= 0
        % and the objective function, without the E(V(x)) : the toolbox
        % automatically takes care of that
        %           min C * x
        if(scenario.time == 1)
            % Constraint
            cntr = x >= 0 ;
            % Objective
            obj  = C * x ;
        end
        
        % If at stage 2, the problem is
        %           min - P s
        %           s.t.  s >= 0
        %                 s <= d(xi)
        %                 s <= x
        % (dy default we minimize, so the '-' is needed)
        if(scenario.time == 2)
            % Constraints
            pos             = s >= 0 ;
            notExceedDemand = s <= scenario.data ; % data is the random 
                                                   % demand store into each 
                                                   % node (scenario)
            notExceedStock  = s <= x ;
            % Concatenate all constraints together
            cntr = [notExceedDemand ; notExceedStock ; pos];
            % The objective
            obj = -S * s;
        end
    end



fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('   Step (4) : Run SDDP!')
disp('   ------------------------------')
disp(' ')
disp('We are now abble to solve our problem with SDDP algorithm.')
disp('The followings are the settings.')
disp('params = sddpSettings(''yalmip.declareVariables'',@declareVariables,...')
disp('    ''algo.McCount'',2, ...')
disp('    ''algo.iterationMax'',10,...')
disp('    ''solver'',''linprog'') ;')
disp(' ')
disp('You can leave the default one, but it is a good idea to take a look')
disp('at the documentation to see what all of them have to do with')
disp('>> help sddpsettings')
disp(' ')
disp('Now, just type')
disp('>> output = sddp(@ndls,lattice,params) ;')
disp(' ')
disp('Press any key to launch this call...')
disp(' ')
disp(' ')
pause
idx=idx+1;


params = sddpSettings('yalmip.declareVariables',@declareVariables,...
    'algo.McCount',2, ...
    'algo.iterationMax',10,...
    'solver','linprog') ;
output = sddp(@ndls,lattice,params) ;



fprintf('\n\n%g / %g \n',idx,tot)
disp(' ')
disp('Your problem is now solved!. You can visualise the evolution of ')
disp('the mean cost and the lower bound. They should meet at the end.')
disp(' ')
disp('>> plotOutput(output) ; ')
disp(' ')
disp('This is the end of the demo. You can find more information with')
disp('the documention.')
disp('>> help sdpsettings')
disp('>> help sddp')
disp(' ')
disp('Pres any key to end the demo.')
pause
plotOutput(output);
end






