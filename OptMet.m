%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         OptMet - Optimal design of resonant metamaterial beams
%                           3-D PRINTED BEAMS
% Sergio Cantero Chinchilla
% V01 - 06/07/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code initialisation 
restoredefaultpath
clearvars; close all; clc

% Frequency range for the vibration attenuation design of the metamaterial
% beam in [Hz]:
fmin=280;
fmax=380;

% Discretisation of the mass percentage variable:
m=0.05:0.05:1;

% Load or create a new set of samples from the prior PDF p(\th)
%
% Definition of the prior information of \theta -(Comment if already created)
% Young's modulus - Gaussian distributed -(Experimentally obtained)
E_mean=1.6217*1e+3; 
E_std=49.8990; 
% Density - Gaussian distributed -(Experimentally obtained)
rho_mean=948.9624;
rho_std=7.3896;
% Generate n_samples samples
n_samples=100;
genSamples(n_samples,E_mean,E_std,rho_mean,rho_std); 

% Load the samples:
load('./dat/th.mat')

% Trade-off variable between expectation and variance in [0,1]:
A = 0.5;

% Maximum number of resonators considered:
MaxRes=15;
ResRange=1:MaxRes;

% Cost function
C=pchip([1,8,12,15],[0,.3,.65,.9],ResRange);


% Initialise variables:
U_p=zeros(length(m),MaxRes);
U_mean=zeros(length(m),MaxRes);
U_var=zeros(length(m),MaxRes);
U=zeros(length(m),MaxRes);
idxMIN_Up=zeros(1,MaxRes);
idxMIN_Um=zeros(1,MaxRes);
idxMIN_Uv=zeros(1,MaxRes);

% Exhaustive search for the optimal design:
for n=ResRange
    
    % Parallel for to ease the computational burden
    parfevalOnAll(gcp(), @warning, 0, 'off', 'MATLAB:nearlySingularMatrix');
    parfor i=1:length(m)
        [U_p(i,n), U_mean(i,n), U_var(i,n)] = ObjFun(m(i),...
            n, th, A, fmin, fmax);
        fprintf('Number of resonators: %d; m=%f\n',n,m(i))
    end
    
    % Trade-off variable (performance - cost)
    B=abs(U_p(:,n));
    
    % Objective function evaluation:
    U(:,n)=U_p(:,n)+B.*C(n);
    
    % Plot results
    [~,idxMIN_Up(n)]=min(U_p(:,n));
    [~,idxMIN_Um(n)]=min(U_mean(:,n));
    [~,idxMIN_Uv(n)]=min(U_var(:,n));
    figure; 
    plot(m,U_p(:,n),'-k',m,U_mean(:,n),'--k',m,U_var(:,n),':k')
    hold on
    plot(m(idxMIN_Up(n)),U_p(idxMIN_Up(n),n),'ok')
    plot(m(idxMIN_Um(n)),U_mean(idxMIN_Um(n),n),'sk')
    plot(m(idxMIN_Uv(n)),U_var(idxMIN_Uv(n),n),'dk')
    xlim([0, max(m)])
    xlabel(strcat({'Mass percentage of '}, num2str(n),{' resonators'}),...
        'interpreter','latex','fontsize',10)
    ylabel(strcat('Sum of FRF in [',num2str(fmin),',',num2str(fmax),'] Hz'),...
        'interpreter','latex','fontsize',10)
    hold off
    set(gcf, 'Units', 'centimeters', 'OuterPosition', [12, 10.3, 12, 10]);
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    set(gca,'TickLabelInterpreter','latex','fontsize',10)
    legend({'U-p','U-mean','U-var'},'location','best','interpreter','latex',...
        'fontsize',10)
    print(gcf,strcat('./res/N_res',num2str(n),'.pdf'),'-dpdf')
    close all
    % Save results
    save(strcat('./res/NRes_',num2str(n),'_Rng.mat'),'U','U_p','U_mean',...
        'U_var','idxMIN_Up')
end

% Optimal design:
[row,col]=find(U==min(U(:)));
OptMass=m(row);
OptNum=ResRange(col);
ObjFunVal=U(row,col);
fprintf('\n\nThe robust optimal design is:\nMass=%f;\nNumber of resonators=%d;\nObjective function value=%f\n',OptMass,OptNum,ObjFunVal)

% Plot the objective function and identify the optimal value:
figure; imagesc(ResRange,m,U);set(gca,'YDir','normal'); hold on;
xlabel('Number of resonators','interpreter','latex','fontsize',10)
ylabel('Mass percentage','interpreter','latex','fontsize',10)
colormap(flipud(gray))
scatter(OptNum,OptMass,30,'o','filled')
h=colorbar; h.TickLabelInterpreter='latex';
set(gcf, 'Units', 'centimeters', 'OuterPosition', [12, 10.3, 12, 10]);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gca,'TickLabelInterpreter','latex','fontsize',10)
print(gcf,'OptDesign','-dpdf')

% Save the properties of the optimal design
save('OptDesign.mat','OptMass','OptNum','ObjFunVal')

% Close the parallel pool
poolobj=gcp('nocreate');
delete(poolobj);