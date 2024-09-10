function[] = PDE4(param_d,param_b,param_k12,param_ts,param_Dm,param_r,param_h2,param_h3,param_delta,param_S,Simul_id)

Dc = param_d;
bc = param_b;
det=5;
delta = param_delta;  % Baba et al. (2008) #0.03983 
K = 100;  % Tumor Carrying capacity => Normal value is 100 cells/mm
Kmm= 30;
%Switching terms
k12 = param_k12;   
k21 = 2.5e-19; 
gamman = 1.01; %Pietro 2019
tn =0.5;     %Regularization term
ts = param_ts;
d1 = 0.87; %Deathrate wang et al. (2012)
d2 = 0.1;
Dm = param_Dm;      %N/A Macrophages diffusion rate
r = param_r;  %Macrophage division rate  r<1 causing instability
%Oxygen
n0 = 1;   %(nmol/mm) i.c. for oxygen concentration
h1 = 3.37e-1; %(1/days), scientific reports Juan Carlos from Bam3 
%---------------------CHANGE_THIS--------------------------------------
h2 = param_h2;%1.14e-1;%+)/2;  %mm/cells.days scientific reports Juan Carlos%[5.73e-3,1.14e-1]
Do = 1.51e2; %mm^2/day, scientific reports Juan Carlos
ncr = 2.5e-1; %Oxygen concentration threshold for hypoxia
Kp = 0.5*K;
S= param_S*Kmm;% based on the carrying capacity of microvessels in tumor J.C.L 2016
h3 = param_h3;

th = 30; % (days) time discretization
Tf = 3 * 3.60e+2; % (days) Final time of the simulation
L = 300; % (mm) length of the domain
xh = 5e-1; % (mm) spatial discretization
time = 0:th:Tf; % Time-steps at which the solutions are returned
xmesh = 0:xh:L; % Discretized spatial domain
m = 0;
opt = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
% Generate a random resection time between 12 and 18 months
resectionTimeMonths = 12 + (18-12) * rand();
resectionTimeDays = resectionTimeMonths * 30; % Convert months to days

% Find the closest time index for resection
[~, resectionIndex] = min(abs(time - resectionTimeDays));
sol1 = pdepe(m, @(x, t, u, DuDx) pdesyscoef(x, t, u, DuDx), @pdesysic, @pdesysbc, xmesh, time(1:resectionIndex), opt);
solAtResection = sol1(end, :, :); % Solution at resection time
solAtResection(:, :, 1) = 0; % Apply resection by setting tumor density to zero
% Assuming the first component of sol1(end, :, 1) is tumor density
tumorDensityBeforeResection = sol1(end, :, 1);
resectionThreshold = det;
exceedingIndices = find(tumorDensityBeforeResection > resectionThreshold);
if isempty(exceedingIndices)
    % Handle the case where no resection is needed based on the criteria
    x_resect = []; % This could be a default value or lead to skipping resection
else
    % For simplicity, take the first index where resection criteria are met
    % You might want to adjust this based on your specific needs
    x_resect = xmesh(exceedingIndices(end));
end
Xmesh2 = x_resect:xh:L;
sol2 = pdepe(m, @(x, t, u, DuDx) pdesyscoef(x, t, u, DuDx), @pdesysic2, @pdesysbc, Xmesh2, time(resectionIndex:end), opt); 
p = sol1(:,:,1); % Extract the p variable
m1 =sol1(:,:,2);
m2 =sol1(:,:,3);
o2 = sol1(:,:,4);

pr = sol2(:,:,1); % Extract the p variable
m1r =sol2(:,:,2);
m2r =sol2(:,:,3);
o2r = sol2(:,:,4);

% Create a structure to hold both parameters and results
simulationResults.Parameters = struct('param_d', param_d, 'param_b', param_b, 'param_k12', param_k12, 'param_ts', param_ts, 'param_Dm', param_Dm, 'param_r', param_r, 'param_h2', param_h2, 'param_h3', param_h3, 'param_delta', param_delta, 'param_S', param_S);
simulationResults.Results.p = p;
simulationResults.Results.m1 = m1;
simulationResults.Results.m2 = m2;
simulationResults.Results.o2 = o2;

simulationResults2.Parameters = struct('param_d', param_d, 'param_b', param_b,  'param_k12', param_k12, 'param_ts', param_ts, 'param_Dm', param_Dm, 'param_r', param_r, 'param_h2', param_h2, 'param_h3', param_h3, 'param_delta', param_delta, 'param_S', param_S, 'Xmesh2', Xmesh2);
simulationResults2.Results.pr = pr;
simulationResults2.Results.m1r = m1r;
simulationResults2.Results.m2r = m2r;
simulationResults2.Results.o2r = o2r;



% Save the structure to a .mat file
filename = sprintf('E:/scene5/SimulationResults_SimulID_%d.mat', Simul_id);
save(filename, 'simulationResults');
filename = sprintf('E:/scene5/SimulationResults2_SimulID_%d.mat', Simul_id);

save(filename, 'simulationResults2');

function [c,f,s] = pdesyscoef(xmesh,time,u,DuDx) % Equation to solve
   % macrophage Oxygen consumption
c = [1; 1; 1; 1];
f1 = Dc.*DuDx(1).*((ts.*u(3)./Kmm)+tn.*(gamman-u(4)))./((ts.*u(3)/Kmm)+tn.*gamman)+...
     Dc.*u(1).*((ts.*DuDx(3)./Kmm) - tn.*DuDx(4))./((ts.*u(3)/Kmm)+tn.*gamman)-...
     Dc.*u(1).*(ts.*u(3)./Kmm + tn.*(gamman-u(4))).*ts.*DuDx(3)./Kmm.*((ts.*u(3)/Kmm)+tn.*gamman).^2;
f2 = Dm.*DuDx(2);
f3 = Dm.*DuDx(3);
f4 = Do.*DuDx(4);
s1 = (bc.*u(1).*(1-(u(1)/K)).*(tn.*u(4)./((ts.*u(3)/Kmm)+tn.*gamman))-delta.*u(2).*u(1));
s2 = (u(1)/(u(1)+Kp))*S+(-k12.*u(2).*u(1) + k21.*u(3))-d1.*u(2);
s3 = -(-k12.*u(2).*u(1)+ k21.*u(3))+(r.*u(3).*(1-(u(3)+u(2))/Kmm)./(1+exp(2*1.0e+1.*(u(4)-ncr))))-d2.*u(3);
s4 = h1.*(n0-u(4))-h2.*u(1).*u(4)-h3.*(u(2)+u(3)).*u(4);
% .*(u(1)./(u(1)+Kp))
f = [f1; f2; f3; f4];    
s = [s1; s2; s3; s4]; 
end
% ---------------------------------------------
function u0 = pdesysic(xmesh)
   
        pm0 = 40;  %=> (cells/mm) initial condition for cell density
        m10 = 0; %=> I.C for M1 macrophages
        m200 = 0 ; %=> I.C for M2 macrophages
        n0 = 1; %=> (nmol/mm) i.c. for oxygen concentration
        u0 = [pm0./(1+exp(2*1.0e+1.*(xmesh-0.5))); 
              m10;
              m200;
              n0];
end
function u0 = pdesysic2(xmesh)
        m10 = 0; %=> I.C for M1 macrophages
        m200 = 0 ; %=> I.C for M2 macrophages
        n0 = 1; %=> (nmol/mm) i.c. for oxygen concentration
        u0 = [det * exp(-sqrt(bc/Dc) * (xmesh-x_resect)); 
              m10;
              m200;
              n0];
end

% ---------------------------------------------
function [pl, ql, pr, qr] = pdesysbc(xl, ul, xr, ur, time)
        % No-flux Boundary condition
        pl = [0;0;0;0];
        ql = [1;1;1;1];
        pr = [0;0;0;0];
        qr = [1;1;1;1];
end
end