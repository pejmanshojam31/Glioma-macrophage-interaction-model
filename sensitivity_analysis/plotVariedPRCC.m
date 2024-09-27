function prccPlot = plotVariedPRCC(M,N,x,labelstring,parameters,prcc)
    figure()
    hold on
    box on
    for mm=1:M
        plot(x, prcc(mm,:),'LineWidth',3.0);
    end
    prccPlot = gca; 
    xlabel('time (days)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('PRCC value', 'FontSize', 14, 'FontWeight', 'bold');
    title(['PRCC Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples'], 'FontSize', 14, 'FontWeight', 'bold');
    legend(parameters.name,'Location','EastOutside', 'FontSize', 14, 'FontWeight', 'bold')
end