% Load the parameters from the file
all_parameters = load('all_parameters_LHS_20000.txt');

% Get the number of patients/simulations from the number of rows in the loaded data
N = size(all_parameters, 1);

% Iterate over each patient/simulation
for i = 1:N
    % Extract parameters for the current patient/simulation
    D = all_parameters(i, 1);
    b = all_parameters(i, 2);
    k12 = all_parameters(i, 3);
    ts = all_parameters(i, 4);
    Dm = all_parameters(i, 5);
    r = all_parameters(i, 6);
    h2 = all_parameters(i, 7);
    h3 = all_parameters(i, 8);
    delta = all_parameters(i, 9);
    S = all_parameters(i, 10);
    
    % Call the PDE2 function with the extracted parameters
    PDE4(D, b, k12, ts, Dm, r, h2, h3, delta, S, i); % Assuming PDE2 function signature matches the order and number of parameters
end