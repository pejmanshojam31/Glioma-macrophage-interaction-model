function LHS2(Np)

N = Np; % Define the number of patients
rng(10);
Dmin = 0.00273000; Dmax = 0.2730000;
bmin = 0.00273000; bmax = 0.02730000;
% Use Latin Hypercube Sampling to generate parameter values
samples = lhsdesign(N, 10); % 10 parameters, N samples
D_list = Dmin + (Dmax-Dmin).*samples(:,1)';
b_list = bmin + (bmax-bmin).*samples(:,2)';
k12 = 2.5e-5+(2.5e-2-2.5e-5).*samples(:,3)';
ts = 0.1 +(0.6-0.1).*samples(:,4)';
Dm = 1e-3 +(1e-1-1e-3).*samples(:,5)';
r = 0.4+ (0.8-0.4).*samples(:,6)';
h2 = 5.73e-3+(1.14e-1-5.73e-3).*samples(:,7)';
h3 = 1e-3+(1e-2-1e-3).*samples(:,8)';
delta = 1e-6 +(1e-4-1e-6).*samples(:,9)';
S = 0.1 +(0.4-0.1).*samples(:,10)';

all_parameters = [D_list; b_list; k12; ts; Dm; r; h2; h3; delta; S]';

% Save the matrix to a file
filename = 'all_parameters_LHS_20000.txt'; % Saves the file in MATLAB's current working directory
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
