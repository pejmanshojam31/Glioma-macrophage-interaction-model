function[] = LHSPRCCIW2(M,N)
%global d1 Dc bc NN k12 k21 Kpm gamman tn ts d Dm r K n0 h1 h2 Do ncr 
%%_________________________________________________________________________
%
% FIRST STEP: Latin Hypercube Sampling of Model Parameters 
%__________________________________________________________________________
% This requests user input via the command line to get started:
M = input('Number of parameters to sample?: '); % Don't include add'l parameters you would like to leave fixed
while rem(M,1)~=0 || M<=0
    M = input('Number of parameters should be an integer > 0. Please re-enter the number of parameters: ');
end
N = input('Number of samples to draw? (recommend 100 to 1000): '); 
while rem(N,1)~=0 || N<=0
    N = input('Number of samples should be an integer > 0. Please re-enter the number of samples: ');
end
% This code will prompt for sample distribution specifics and return drawn and randomly paired parameter samples
[parameters,A] = DrawSamples(M,N);
outFileStr = 'LHS-testlinear'; % give workspace an appropriate unique name
outFileName1 = [outFileStr,'_samples.mat'];
outFileName2 = [outFileStr,'_results.mat'];
save(outFileName1, 'parameters', 'A')
binCount = ceil(N/10); 
histPlot = plotSampleHists(M,parameters,binCount);
pause(1) %Time to dock/maximize the figure before it saves
% Name and save histograms visualizing the sample distributions:
figurelabel1=([outFileStr,'-N',num2str(N),'-histograms.pdf']);
% figurelabel2=([outFileStr,'-N',num2str(N),'-histograms.png']);
saveas(histPlot, figurelabel1);
% saveas(histPlot, figurelabel2);
%%_________________________________________________________________________
%
% SECOND STEP: Solve Model Function with Sampled Parameters
%__________________________________________________________________________
m = 0;
th    = 30;           %=> (days) time discretization 
Tf    = 3 * 3.60e+2;  %=> (days) Final time of the simulation
L     = 300;          %=> (mm) length of the domain
xh    = 7e-1;       %=> (mm) spatial discretization
time  = 0:th:Tf;      %=> Time-steps at which the solutions are returned
xmesh = 0:xh:L;       %=> Discretized spatial domain
% delta = 2e-5;  % Baba et al. (2008) #0.03983 
% Dc = 0.373;
% bc = 0.2;
K = 1e2;  % Tumor Carrying capacity => Normal value is 100 cells/mm
Kmm = 30;
% k12 = 2.5e-1;   
k21 = 2.5e-19; 
gamman = 1.01; %Pietro 2019
% tn = 0.5;     %Regularization term
% %Macrophages
d1 = 0.87; %Deathrate wang et al. (2012)
d2 = 0.1;
%Oxygen
Kp = 0.1*K;
% S = 5;
n0 = 1.0;   %(nmol/mm) i.c. for oxygen concentration
h1 = 3.37e-1; %(1/days), scientific reports Juan Carlos from Bam3 
% h2 = 1.14e-1;%1.14e-1;%+)/2;  %mm/cells.days scientific reports Juan Carlos%[5.73e-3,1.14e-1]
Do = 1.51e2; %mm^2/day, scientific reports Juan Carlos
ncr = 2.5e-1; %Oxygen concentration threshold for hypoxia
% S = 0.2;% based on the carrying capacity of microvessels in tumor J.C.L 2016
% h3 = 3e-2;   % macrophage Oxygen consumption
%==>> effective diffusion and proliferation parameters
% alpha = (ts*(m20.^2/(m20.^2+Kpm.^2))+tn *((gamman-n0)/(tn*gamman+ts)));
% betta = (ts*(1-(m20.^2/(m20.^2+Kpm.^2))+tn *(n0)/(tn*gamman+ts)));
% q.Dc = Dc/alpha;
% q.bc = bc/betta;
Simdata = struct; % initialize struct to store outputs for computing PRCCs
tic % start measuring time to solve equations to monitor progress
% Loop over the parameter sample pairs for Monte Carlo simulations
for j=1:N   
    sampledparams=A(j,1:M);
    fprintf('Parameters passed for sample: ');fprintf('%u',j);fprintf(' of ');fprintf('%u\n',N);
%          Simdata(j).p = ThirdPDE(xmesh, time, sampledparams, d1,Dc, bc, NN, k12, k21, Kpm, gamman, tn, ts, d, m20, Dm,...
%              r,K, n0, h1, h2, Do, ncr);
%        Simdata(j).p = lastPDE(xmesh, time, sampledparams, d1,Dc, bc, NN, k12, k21, Kpm, gamman, tn, ts, d, Dm,...
%      r,K, n0, h1, h2, Do, ncr);
%     lastPDE(xmesh, time,sampledparams, d1, Dc, bc, NN, k12, k21, Kpm, gamman, tn, ts, d, Dm,...
%                r, K, n0, h1, h2, Do, ncr)
    delta = sampledparams(1);%2e-3;  % (Normal) Baba et al. (2008) #0.03983 sampledparams, N, Kpm, m20, n0, Do, ncr, v  
    Dc = sampledparams(2);%Dmax;
    bc = sampledparams(3);
    k12 = sampledparams(4);%2.5;   
%     k21 = sampledparams(5);%2.5e-1;
%     Kpm =sampledparams(6);
% %     gamman = sampledparams(6);%1.01; %Pietro 2019
    tn = sampledparams(5);%2.0;     %Regularization term
    ts = sampledparams(6);%1.0;   % regularization term
%     d1 = sampledparams(6);%0.75; %Deathrate wang et al. (2012)
%     d2 = sampledparams(7);%0.75; %Deathrate wang et al. (2012)
    Dm = sampledparams(7);%5e-2;      %N/A Macrophages diffusion rate
    r = sampledparams(8);%1.5;  %Macrophage division rate  
%     h1 = sampledparams(9);%3.37e-1; %(1/days), scientific reports Juan Carlos from Bam3 
    h2 = sampledparams(9);%(5.73e-3+1.14e-1)/2;  %mm/cells.days scientific reports Juan Carlos%[5.73e-3,1.14e-1]
%     Do = sampledparams(7);%1.51e2; %mm^2/day, scientific reports Juan Carlos
% %     v = sampledparams(16) ; % based on the carrying capacity of microvessels in tumor J.C.L 2016
    h3 = sampledparams(10);   % macrophage Oxygen consumption
%     S = sampledparams(10); 
    S = sampledparams(11)*Kmm;
%  EDIT THE FOLLOWING FUNCTION CALL FOR YOUR OWN MODEL:
        opt = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
        sol = pdepe(m, @pdesyscoef, @pdesysic, @pdesysbc, xmesh, time,opt);
        Simdata(j).p = sol(:,:,1);
        Simdata(j).m1 = sol(:,:,2);
        Simdata(j).m2 = sol(:,:,3);
        Simdata(j).n = sol(:,:,4);
    toc
end
save(outFileName2, 'Simdata') % add Simdata to the saved .mat file
OutputOfInterest = zeros(N,length(time),length(xmesh));
InfiltrationWidths = zeros(N, length(time));
TumorMass = zeros(N, length(time));
% tumor_concentration = zeros(N, length(time));
for si = 1:N
    for ti = 1:length(time)
        vals = Simdata(si).p(ti, :); % Extract glioma density at a specific time
        InfiltrationWidths(si, ti) = func_InfiltrationWidth(vals, xmesh);
        TumorMass(si, ti) = func_Mass(vals, xmesh);

    end
end
save(outFileName2, 'Simdata', 'OutputOfInterest', 'InfiltrationWidths','TumorMass');
labelstring = 'Glioma concentration';
    prccIW = VariedPRCC(M, N, A, InfiltrationWidths, time, []);
    prccTM = VariedPRCC(M, N, A, TumorMass, time, []);
    prccPlotIW = plotVariedPRCC(M, N, time, 'IW', parameters, prccIW);
    prccPlotTM = plotVariedPRCC(M, N, time, 'TS', parameters, prccTM);

figurelabel1 = ([outFileStr, '-N', num2str(N), '-InfiltrationWidth-PRCC.fig']);
saveas(prccPlotIW, figurelabel1);
figurelabel3 = ([outFileStr, '-N', num2str(N), '-TumorMass-PRCC.fig']);
saveas(prccPlotTM, figurelabel3);


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
% ---------------------------------------------
    function [pl, ql, pr, qr] = pdesysbc(xl, ul, xr, ur, time)
            % No-flux Boundary condition
            pl = [0;0;0;0];
            ql = [1;1;1;1];
            pr = [0;0;0;0];
            qr = [1;1;1;1];
    end
end