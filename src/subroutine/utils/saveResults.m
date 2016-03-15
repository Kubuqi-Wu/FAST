function saveResults(solutionForward, lowerBounds, meanCosts, stds, iteration, lattice, params, filename)

% A little be redundant, no ?
output.lattice     = lattice ;
output.lowerBounds = lowerBounds ;
output.meanCost    = meanCosts ;
output.stds        = stds ;
output.params      = params ;

save(filename,'output','solutionForward','lowerBounds','meanCosts', 'stds','iteration','lattice','params') ;

end