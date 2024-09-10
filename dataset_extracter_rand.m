L     = 300;    % (mm) length of the domain
xh    = 7e-1; % (mm) spatial discretization
xmesh = 0:xh:L; % Discretized spatial domain
NSiml = 20000;     % Number of patients
all_parameters = cell(NSiml, 1);
IW_at1 = zeros(NSiml, 1); 
IW_at2 = zeros(NSiml, 1); 
IW_at3 = zeros(NSiml, 1); 
random_values_m1 = zeros(NSiml,1);
random_values_m2 = zeros(NSiml,1);
random_values_p = zeros(NSiml,1);
mass_80 =zeros(NSiml,1);
mass_16 =zeros(NSiml,1);
ki67 = zeros(NSiml,1);
column_headers = {'param_b', 'param_d', 'param_h2','param_r', 'param_h3','param_Dm','param_ts','param_S', 'IW_at1', 'IW_at2','IW_at3','ki67',...
                  'm1_1','m2_1','p1','T1_MRI','FLAIR_MRI'};
all_data_for_csv = column_headers; % This will be your first row
for k = 1:NSiml
        simul_pdesys = load(strcat('E:/scene5/simulationResults_SimulID_', num2str(k), '.mat')); % Assuming 'k' indexes patients/simulations
        simul_pdesys2 = load(strcat('E:/scene5/simulationResults2_SimulID_', num2str(k), '.mat'));
        simulationResults = simul_pdesys.simulationResults;
        results = simulationResults.Results;
        simulationResults2 = simul_pdesys2.simulationResults2;
        results2 = simulationResults2.Results;
        parameters = simulationResults.Parameters;
        [ki67(k,:)] = solveForPI(parameters.param_b);
        Xmesh2 = simulationResults2.Parameters.Xmesh2;
        [mass_80(k,:), mass_16(k,:)] = func_Mass(results.p(end,:), xmesh);
        [IW_at1(k,:)] = func_IW(results2.pr(3,:), Xmesh2);
        [IW_at2(k,:)] = func_IW(results2.pr(6,:), Xmesh2);
        [IW_at3(k,:)] = func_IW(results2.pr(12,:), Xmesh2);
        [random_values_m1(k,:),random_values_m2(k,:), random_values_p(k,:)] = func_InfiltrationWidth(results.p(end,:),results.m1(end,:),results.m2(end,:), xmesh);
        param_values = [parameters.param_b, parameters.param_d, parameters.param_h2,parameters.param_r  ,parameters.param_h3  ,parameters.param_Dm, parameters.param_ts, parameters.param_S  ];
        numeric_data_row = num2cell([IW_at1(k,:), IW_at2(k,:), IW_at3(k,:), ki67(k,:), random_values_m1(k,:), random_values_m2(k,:), random_values_p(k,:), mass_80(k,:), mass_16(k,:)]);
        combined_row = [num2cell(param_values), numeric_data_row];
        all_data_for_csv = [all_data_for_csv; combined_row];
end
writecell(all_data_for_csv, 'random_dataset_8Apr.csv');
function [random_values_m1, random_values_m2, random_values_p] = func_InfiltrationWidth(vals, vals1, vals2, xmesh)
    % Initialize variables
    ind99 = [];
    ind01 = [];
    maxVal = max(vals); % Compute once to optimize

    for k = 1:length(vals)
        if isempty(ind99) && vals(1,k) >= (0.99 * maxVal)
            ind99 = k; % First index where condition is met
        end
        if vals(1,k) >= (0.05 * maxVal)
            ind01 = k; % Continuously update to last index meeting condition
        end
    end

    % Ensure ind80 and ind01 are found and correct
    if isempty(ind99) || isempty(ind01)
        error('Ind99 or Ind01 not found. Check the values array for adequate range.');
    end

    % Set random seed for reproducibility
    rng(123);

    % Generate a random index between ind80 and ind01
    random_index = randi([ind99, ind01], 1);

    % Get random values based on the generated index
    % Ensuring index is within bounds for xmesh
    xmesh_idx = max(1, min(length(xmesh), round(xmesh(1,random_index))));
    random_values_m1 = vals1(end, xmesh_idx); % Random values from m1
    random_values_m2 = vals2(end, xmesh_idx); % Random values from m2
    random_values_p = vals(end, xmesh_idx); % Random values for p
end


function [IW_at] = func_IW(vals, xmesh)
    for k = 1:length(vals)
        if(vals(1,k) >= (0.80 * max(vals)))
            ind80 = k;
        end
    end
    for k = length(vals):-1:1
        if(vals(1,k) <= (0.02 * max(vals)))
            ind02 = k;
        end
    end
    IW_at = xmesh(1,ind02) - xmesh(1,ind80);
end



function [ki67] = solveForPI(a_val)
    syms a t ta AI PI
    equation = a_val * (1 - AI - PI) - ((1/t) * (PI + PI^2) - (1/ta) * AI * PI) == 0;
    t_val = 24/24;
    ta_val = 8/24;
    AI_val = 0.7/100;
    equation_substituted = subs(equation, [a, t, ta, AI], [a_val, t_val, ta_val, AI_val]);
    solutions = solve(equation_substituted, PI);
    positive_solution = max(double(solutions));
    if positive_solution < 0
        error('No positive solution found.');
    else
        ki67 = positive_solution*100;
    end
end

function [mass_80, mass_16] = func_Mass(p, xmesh)
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
    start_index = find(xmesh >= 0, 1);
    if isempty(start_index)
        start_index = 1; 
    end
    mass_80 = trapz(xmesh(start_index:ind_80), p(start_index:ind_80) / max_p);
    mass_16 = trapz(xmesh(start_index:ind_16), p(start_index:ind_16) / max_p);
end
