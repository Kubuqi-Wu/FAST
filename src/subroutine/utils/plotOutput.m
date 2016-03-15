function fig = plotOutput(output)

Mc = output.params.algo.McCount ;
fig = figure;
hold on;
plot(output.lowerBounds,'linewidth',1.5) ;
plot(output.meanCost,'r','linewidth',1.5) ;
plot(output.meanCost-output.stds/sqrt(Mc),'--r','linewidth',1) ;
plot(output.meanCost-2*output.stds/sqrt(Mc),'-.r','linewidth',0.5) ;
plot(output.meanCost+output.stds/sqrt(Mc),'--r','linewidth',1) ;
plot(output.meanCost+2*output.stds/sqrt(Mc),'-.r','linewidth',0.5) ;
legend('Lower bound','Mean Costs','67% confidence interval','95% confidence interval')
xlabel('Iteration');
ylabel('Value')