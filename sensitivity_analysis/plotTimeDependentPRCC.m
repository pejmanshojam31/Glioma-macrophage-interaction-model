function prccPlot = plotTimeDependentPRCC(timePoints, prcc, labelstring, paramNames)
    % Create a new figure
    figure();
    hold on;
    box on;

    % Plot PRCC for each parameter
    for mm = 1:size(prcc, 1)
        plot(timePoints, prcc(mm, :), 'LineWidth', 2.0);
    end

    % Set plot properties
    prccPlot = gca;
    xlabel('Time');
    ylabel('PRCC value');
    title(['Time-dependent Correlation Plot of ', labelstring, ', ', num2str(size(prcc, 2)), ' time points']);

    % Create a legend with parameter names
    if exist('paramNames', 'var')
        legend(paramNames, 'Location', 'EastOutside');
    else
        legend(arrayfun(@(x) ['Param ', num2str(x)], 1:size(prcc, 1), 'UniformOutput', false), 'Location', 'EastOutside');
    end

    hold off;
end
