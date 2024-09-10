function gen_Pat2(Np)

N = Np; % Define the number of patients

rng(10); % Fix the seed for rand (to make results reproducible)

% Set the max and min values of the parameters for the full model
Dmin = 0.0273000; Dmax = 0.2730000;
bmin = 0.00273000; bmax = 0.02730000;
% Generate the parameters for the full model from a uniform distribution
% Generate the list of motility D
D_list = Dmin + (Dmax-Dmin).*rand(1,N);
% % Gen. list of proliferation b
b_list = bmin + (bmax-bmin).*rand(1,N);
k12 = 2.5e-5+(2.5e-2-2.5e-5).*rand(1,N);
ts = 0.1 +(0.6-0.1).*rand(1,N);
Dm = 1e-3 +(1e-1-1e-3).*rand(1,N);
r = 0.4+ (0.8-0.4).*rand(1,N);
h2 = 1e-4+(1e-2-1e-4).*rand(1,N);
h3 = 1e-4+(1e-2-1e-4).*rand(1,N);
delta = 1e-6 +(1e-4-1e-6).*rand(1,N);
S = 0.1 +(0.4-0.1).*rand(1,N);

all_parameters = [D_list; b_list; k12; ts; Dm; r; h2; h3; delta; S]';
% Save the matrix to a file
filename = 'all_parameters_test.txt'; % Saves the file in MATLAB's current working directory
% Simplified file path for testing
fileID = fopen(filename,'w');
if fileID == -1
    error('Cannot open file for writing. Check permissions and path.');
end
for i = 1:size(all_parameters,1)
    fprintf(fileID, '%1.4e\t', all_parameters(i,:));
    fprintf(fileID, '\n');
end
fclose(fileID);