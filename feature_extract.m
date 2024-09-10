function feature_extract(Np)
    L     = 300;    % (mm) length of the domain
    xh    = 7e-1; % (mm) spatial discretization
    xmesh = 0:xh:L; % Discretized spatial domain
    NSiml = Np;     % Number of patients
    all_parameters = cell(NSiml, 1);
    % The dimensions are [NSiml, 1] initially, since the time step j is chosen randomly for each patient
    m1_values_80_60 = cell(NSiml, 1);
    m2_values_80_60 = cell(NSiml, 1);
    m1_values_40_20 = cell(NSiml, 1);
    m2_values_40_20 = cell(NSiml, 1);
    m1_values_10_02 = cell(NSiml, 1);
    m2_values_10_02 = cell(NSiml, 1);
% Preallocate for results storage
    resultsTable = table([], [], [], [], [], [], [], [], [], [], [], [], ...
        'VariableNames', {
            'SimulationID', 'TimeStep', 'Parameters', 'M1_Values_80_60', 'M2_Values_80_60',... 
            'M1_Values_40_20', 'M2_Values_40_20', 'M1_Values_10_02', 'M2_Values_10_02',... 
            'IW_ec', 'Mass_80', 'Mass_16'
        });

j = randi([4, 16]);
% Loop through each simulation (patient)
for k = 1:NSiml
        % Randomly select a time step between 4 and 16 for this patient

        
        % Load the simulation results for this patient
        simul_pdesys = load(strcat('C:/Users/pejma/Nextcloud/PDE model/Final2/simulationResults_SimulID_', num2str(k), '.mat')); % Assuming 'k' indexes patients/simulations
        simulationResults = simul_pdesys.simulationResults;
        results = simulationResults.Results;
        parameters = simulationResults.Parameters;
        parameters_str = jsonencode(parameters);
        % Perform calculations at the selected time step
        % Note: Adjust the func_InfiltrationWidth and func_Mass function calls as necessary based on their definitions
        [IW_ec(k,j), random_xmesh80_60, m1_80_60, m2_80_60, random_xmesh40_20, m1_40_20, m2_40_20, random_xmesh10_02, m1_10_02, m2_10_02] = func_InfiltrationWidth(results.p, results.m1, results.m2, xmesh);
        [mass_80(k),mass_16(k)] = func_Mass(results.p(j,:), xmesh);
        
        % Store m1 and m2 values for each range for this patient and time step
        m1_values_80_60{k, j} = m1_80_60;
        m2_values_80_60{k, j} = m2_80_60;
        m1_values_40_20{k, j} = m1_40_20;
        m2_values_40_20{k, j} = m2_40_20;
        m1_values_10_02{k, j} = m1_10_02;
        m2_values_10_02{k, j} = m2_10_02;
                % Flatten or serialize m1, m2 values for CSV
        m1_80_60_str = mat2str(m1_80_60);
        m2_80_60_str = mat2str(m2_80_60);
        m1_40_20_str = mat2str(m1_40_20);
        m2_40_20_str = mat2str(m2_40_20);
        m1_10_02_str = mat2str(m1_10_02);
        m2_10_02_str = mat2str(m2_10_02);
        % Store or process other necessary parameters or results as needed
        all_parameters{k} = parameters;
               % Append results to the table
        newRow = {k, j, parameters_str, m1_80_60_str, m2_80_60_str, m1_40_20_str, m2_40_20_str, m1_10_02_str, m2_10_02_str, IW_ec, mass_80, mass_16};
        tempTable = cell2table(newRow, 'VariableNames', {
    'SimulationID', 'TimeStep', 'Parameters', 'M1_Values_80_60', 'M2_Values_80_60',... 
    'M1_Values_40_20', 'M2_Values_40_20', 'M1_Values_10_02', 'M2_Values_10_02',... 
    'IW_ec', 'Mass_80', 'Mass_16'
});
        resultsTable = [resultsTable; tempTable];
            % Write results to CSV
        writetable(resultsTable, 'patient_simulation_results.csv');
end

function [tfl, random_xmesh80_60, m1_80_60, m2_80_60, random_xmesh40_20, m1_40_20, m2_40_20, random_xmesh10_02, m1_10_02, m2_10_02] = func_InfiltrationWidth(p, m1, m2, xmesh)
    % Calculate the maximum value of p for threshold calculations
    maxVal = max(p);
    
    % Initialize indices to NaN; they'll be overwritten when found.
    thresholds = [0.80, 0.60, 0.40, 0.20, 0.10, 0.02];
    inds = NaN(1, length(thresholds));
    
    % Find indices for each threshold
    for i = 1:length(thresholds)
        tempInd = find(p >= thresholds(i) * maxVal, 1, 'first');
        if ~isempty(tempInd)
            inds(i) = tempInd;
        end
    end
    
    % Extract indices for clarity
    ind80 = inds(1);
    ind60 = inds(2);
    ind40 = inds(3);
    ind20 = inds(4);
    ind10 = inds(5);
    ind02 = inds(6);
    
    % Calculate tfl based on the found indices
    tfl = xmesh(ind02) - xmesh(ind80); % Use xmesh values for tfl calculation
    
    % Helper function to select random values
function [rxmesh, rm1, rm2] = selectRandom(startInd, endInd, n)
    if ~isnan(startInd) && ~isnan(endInd)
        range = endInd - startInd + 1;
        selectedInds = randperm(range, min(n, range)) + startInd - 1;
        rxmesh = xmesh(selectedInds);
        rm1 = m1(selectedInds);
        rm2 = m2(selectedInds);
    else
        rxmesh = [];
        rm1 = [];
        rm2 = [];
    end
end

    % Select random xmesh values and corresponding m1, m2 values within each specified range
    [random_xmesh80_60, m1_80_60, m2_80_60] = selectRandom(ind80, ind60, 5);
    [random_xmesh40_20, m1_40_20, m2_40_20] = selectRandom(ind40, ind20, 5);
    [random_xmesh10_02, m1_10_02, m2_10_02] = selectRandom(ind10, ind02, 5);
end
function [mass_80, mass_16] = func_Mass(p, xmesh)
    % Calculate the maximum value of p to find 80% and 16% thresholds
    max_p = max(p);
    threshold_80 = 0.80 * max_p;
    threshold_16 = 0.16 * max_p;
    
    for t = 1:length(p)
        if(p(1,t) >= (threshold_80))
            ind_80 = t;
        end
    end
    for t = 1:length(p)
        if(p(1,t) >= (threshold_16))
            ind_16 = t;
        end
    end
    start_index = find(xmesh >= 0, 1); % Find the index where xmesh is zero or the first available point
    if isempty(start_index)
        start_index = 1; % Default to the first index if no exact zero point is found
    end
    
    % Trapz from start_index to ind_80 and ind_16
    mass_80 = trapz(xmesh(start_index:ind_80), p(start_index:ind_80) / max_p);
    mass_16 = trapz(xmesh(start_index:ind_16), p(start_index:ind_16) / max_p);
end
end