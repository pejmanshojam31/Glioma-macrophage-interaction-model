% Example of how to use the function:
plotGliomaDensityAtTime(1080, Simdata, time, xmesh, N); % Plots data for time = 150 days
function plotGliomaDensityAtTime(targetTime, Simdata, time, xmesh, N)
    % Define simulation parameters
    Tf = 3 * 3.60e+2;  % Final time of the simulation

    % Check if the target time is within the simulation time range
    if targetTime < 0 || targetTime > Tf
        error('Target time is outside the simulation time range.');
    end

    % Find the index of the closest time point to the target time
    [~, timeIndex] = min(abs(time - targetTime));

    % Initialize arrays to store the values
    numSpatialPoints = length(xmesh);
    allDensities = zeros(N, numSpatialPoints);
    
    % Extract data for the specified time step
    for i = 1:N
        allDensities(i, :) = Simdata(i).p(timeIndex, :);
    end
    
    % Calculate mean, minimum, and maximum
    meanDensity = mean(allDensities, 1);
    minDensity = min(allDensities, [], 1);
    maxDensity = max(allDensities, [], 1);

    % Plotting
    figure;
    hold on;
    plot(xmesh, meanDensity, 'b', 'LineWidth', 2); % Plot mean density
    fill([xmesh, fliplr(xmesh)], [minDensity, fliplr(maxDensity)], 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0); % Shadowed region
    xlabel('Space (mm)');
    ylabel('Glioma Density');
    title(['Glioma Density vs. Space at Time = ', num2str(targetTime), ' days']);
    hold off;
end

