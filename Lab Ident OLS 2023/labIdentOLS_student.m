% Aerodynamic coefficients Identification with OLS
% This labwork is derived from:
% [1] Flight Vehicle System Identification - A Time Domain Methodology, Second Edition, 
% Ravindra V. Jategaonkar, Published by AIAA, Reston, VA 20191, USA
% 
% Companion software and data available online https://arc.aiaa.org/doi/suppl/10.2514/4.102790
%
% Original data have been pre-processes for educational purpose
%
%                  Yves Briere ISAE 2023
%
%%  Load data and create variables in the workspace
clear,clc;
load('data_ident_AC.mat');
time = (1:max(size(CD)))'*0.04;
% List of data created:
% 
% Outputs (dependant variables) :
% CD    drag coefficient
% CL    lift coefficient
% CY     side force coefficient
% Cm_ac pitching moment coefficient at AC center
% Cl_ac rolling moment coefficient at AC center
% Cn_ac yawing moment coefficient at AC center 

% Inputs (independant variables) :
% de elevator deflection (rad)
% da aileron deflection (rad)
% dr rudder deflection (rad)
% p normalized angular rate (p*b/V) cf p208
% q normalized angular rate (q*l/V) cf p208
% r normalized angular rate (r*b/V) cf p208
% alpha angle of attack at CG (rad)
% beta sidesleep angle at CG (rad)
% V true air speed at CG (rad)
% rho density of air
% Fe left engine thrust
% Fer right engine thrust
% ald alphadot (rad/s)
% M Mach number

% Data have pre-processed with a-priori knowledge about the aircraft (mass, 
% inertia, dimensions, etc.). For instance the CD coefficient was obtained from:
% 
% CD = -Cx cos(alpha)-Cz sin(alpha)
% 
% and Cx = (m acc -Feng)/q/S
% 
% where m is the mass of the aircraft, acc the acceleration of the center of 
% gravity along the x axis, Feng the force of the engine, q the measured dynamic 
% pressure and S the reference area. More details in [1]

%% Data visualization

figure(1);
tiledlayout(6,1)
ax21 = nexttile;plot(time,CD);ylabel('CD');
ax22 = nexttile;plot(time,CL);ylabel('CL');
ax23 = nexttile;plot(time,CY);ylabel('CY');
ax24 = nexttile;plot(time,Cm_ac);ylabel('Cm_ac');
ax25 = nexttile;plot(time,Cl_ac);ylabel('Cl_ac');
ax26 = nexttile;plot(time,Cn_ac);ylabel('Cn_ac');
figure(2);
tiledlayout(3,1)
ax11 = nexttile;plot(time,de,time,da,time,dr);ylabel('de,da,dr (rad)');legend('de','da','dr');
ax12 = nexttile;plot(time,p,time,q,time,r);ylabel('p,q,r (norm)');legend('p','q','r');
ax13 = nexttile;plot(time,alpha,time,beta);ylabel('alpha,beta (rad)');legend('alpha','beta');
xlabel('time (s)');
xlabel('time (s)');
linkaxes([ax11 ax12 ax13 ax21 ax22 ax23 ax24 ax25 ax26 ],'x');

%% Identification of aerodynamic coefficients
% 
% Drag :  CD = CD0 + CDal*alpha
% Lift :  CL = CL0 + Clalpha*alpha
% Pitch : Cm_ac = Cm0 + etc.
% Side :  CY = CY0 + CY_beta*beta
% Roll :  Cl_ac = Cl0 + Cl_beta*beta + Cl_p*p + Cl_r*r + Cl_da*da
% Yaw :   Cn_ac = Cn0 + Cn_beta*beta + Cn_p*p + Cn_r*r + Cn_dr*dr
%
% Modify section 1 (Organize data) only. Sections 2, 3 and 4 automatically
% computes LS estimate and displays relevant results
%

% 1 Organize data
Y_name = 'Cm_ac'; % name of dependant variable
theta_name = {'Cm0','Cm_alpha','Cm_de', 'Cm_da'}; % name of parameters
PHI_name = {'one','alpha','de','da'}; % name of independant variables

T = length(time); %Number of samples
one = ones(T,1);
[~,N]=size(PHI_name); % Number of independant variables

% 2 Create and concatenate variables
Y = eval(Y_name);
PHI = [];
for k = 1:N
    PHI = [PHI eval(PHI_name{k})];
end

% 3 Estimate OLS
thetahat = (PHI'*PHI)^-1*PHI'*Y; % OLS solution
% Compute statistics
Yhat = PHI*thetahat; % Predicted output
sigmahat2 = (Y-Yhat)'*(Y-Yhat)/(T-N); % Estimated standard deviation
R2 = 1-((Y-Yhat)'*(Y-Yhat))/((Y-mean(Y))'*(Y-mean(Y))); % Fit
D = (PHI'*PHI)^-1;
covtheta = sigmahat2*D; % Covariance(theta)

% 4 Extract and display results
equation = [Y_name ' = ' theta_name{1}];
for k = 2:N
    equation = [equation ' + ' theta_name{k} '*' PHI_name{k}];
end
disp(equation)
for k= 1:N
    disp([theta_name{k} ' = ' num2str(thetahat(k)) ' +/- ' num2str(abs(2*sqrt(covtheta(k,k)))) ' at 95% (t = theta/sqrt(cov(theta)) = ' num2str(thetahat(k)/abs(sqrt(covtheta(k,k))))]);
end
disp(['fit: R2 = ' num2str(100*R2) '%']);
