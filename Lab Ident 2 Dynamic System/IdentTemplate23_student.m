%% Identification of linear system from sampled data

%% Load data for identification

clear all;
clc;
load("dataident.mat");
figure;
subplot(211);
ax1 = plot(t,u);ylabel('u PRBS');
subplot(212);
ax2 = plot(t,y);ylabel('y');
xlabel('time');
title('Response to PRBS excitation')

figure;
stairs(tstep,ystep);
xlabel('time');ylabel('y');title('step response');

% Data coming from simulation from a pure linear system:
% y(t ) =  F(s) u(t) + e(t)
%
% First set of data: 120s simulation. The system F(s) is constant during
% the first half of the simulation. There is a change in some parameter(s)
% at the half of the simulation. System switches from F1(s) to F2(s)
%
% Second set of data is the step response of F1(s)
%

%% Data preparation

uident = u(1:300);
yident = y(1:300);
uvalid = u(300:600);
yvalid = y(300:600);

%% Identification of a discrete time model

na = 1;
nb = 1;
nk = 2;
nc = 2;
%% 

na = 0;
nb = 100;
nk = 0;
[A,B] = RLS(uident,yident,na,nb,nk);
B
stem(B)
%% RLS
%intitalized
nk = 0;
na_opt = na;
nb_opt = nb;
J_min = zeros(10,10);
fit_mat = zeros(10,10);
for na = 1:10
    for nb = 1:10
        [A,B] = RLS(uident,yident,na,nb,nk);
        Frls = tf(B,A,dt,'variable','z^-1');
        [err, J, AIC,yhat] = computerrors(uvalid,yvalid,A,B,[]);
        yhat = lsim(Frls,uvalid,0:dt:(length(uvalid)-1)*dt);
        fit_mat(na,nb) = 100*(1-norm(yvalid(100:end)-yhat(100:end))/norm(yvalid(100:end)-mean(yvalid(100:end))));
        na_opt = na;
        nb_opt = nb;
        J_min(na_opt,nb_opt) = J;
    end
end
figure;
surf(J_min)
%optimal value na = 6; nb = 5 from surf plot
na = 7;
nb = 6;
[A,B] = RLS(uident,yident,na,nb,nk);
Frls = tf(B,A,dt,'variable','z^-1');
[err, J, AIC,yhat] = computerrors(uvalid,yvalid,A,B,[]);
yhat = lsim(Frls,uvalid,0:dt:(length(uvalid)-1)*dt);
figure;
plot([yvalid yhat]);
fit = 100*(1-norm(yvalid(100:end)-yhat(100:end))/norm(yvalid(100:end)-mean(yvalid(100:end))));
mse = norm(yvalid(100:end)-yhat(100:end))/length(yvalid(100:end));
legend(['fit = ' num2str(fit) '% mse = ' num2str(mse)]);
title('RLS fit');%Analyse residuals
N = length(uvalid);
x = xcorr(err,err);
x = x/x(N);
figure;
stem(-40:1:40,x(N-40:N+40));
line([-40 40],[2.17/sqrt(N) 2.17/sqrt(N)],'color','r');
line([-40 40],[-2.17/sqrt(N) -2.17/sqrt(N)],'color','r');
title('Discrete model: rls');

%model validation - RLS
figure;
[A,B,C] = ELS(uident,yident,na,nb,nc,nk);
sysid = tf(B,A,dt,'variable','z^-1');
ysim = lsim(sysid,ustep,tstep);
plot(tstep,ystep,tstep,ysim);

%% ELS
nk = 0;
na_opt = na;
nb_opt = nb;
J_min = zeros(10,10);
fit_mat = zeros(10,10);
for na = 1:10
    for nb = 1:10
    [A,B] = ELS(uident,yident,na,nb,nc,nk);
    Fels = tf(B,A,dt,'variable','z^-1');
    [err, J, AIC,yhat] = computerrors(uvalid,yvalid,A,B,[]);
    yhat = lsim(Fels,uvalid,0:dt:(length(uvalid)-1)*dt);
    fit_mat(na,nb) = 100*(1-norm(yvalid(100:end)-yhat(100:end))/norm(yvalid(100:end)-mean(yvalid(100:end))));
    end
end 
figure;
surf(fit_mat)
figure;
plot([yvalid yhat]);
na = 7;
nb = 4;
fit = 100*(1-norm(yvalid(100:end)-yhat(100:end))/norm(yvalid(100:end)-mean(yvalid(100:end))));
mse = norm(yvalid(100:end)-yhat(100:end))/length(yvalid(100:end));
legend(['fit = ' num2str(fit) '% mse = ' num2str(mse)]);
title('ELS fit');

% Analyse residuals
N = length(uvalid);
x = xcorr(err,err);
x = x/x(N);
figure;
stem(-40:1:40,x(N-40:N+40));
line([-40 40],[2.17/sqrt(N) 2.17/sqrt(N)],'color','r');
line([-40 40],[-2.17/sqrt(N) -2.17/sqrt(N)],'color','r');
title('Discrete model: els');


%model validation - ELS
figure;
[A,B,C] = ELS(uident,yident,na,nb,nc,nk);
sysid = tf(B,A,dt,'variable','z^-1');
ysim = lsim(sysid,ustep,tstep);
plot(tstep,ystep,tstep,ysim);
%% Identification of a continuous time model (svf)

lambda = 5;
na = 6;
nb = 4;
nk = 4;
Fsvf = svf(uident,yident,dt,na,nb,nk,lambda);
yhat = lsim(Fsvf,uvalid,0:dt:(length(uvalid)-1)*dt);
yvalid2 = yvalid(100:end); %forget first 100 samples to absorb initial condition mismatch
yhat2 = yhat(100:end);
err = yhat2-yvalid2;
N = length(yvalid2);
figure;
plot([yvalid2 yhat2]);
fit = 100*(1-norm(yvalid2-yhat2)/norm(yvalid2-mean(yvalid2)));
mse = norm(yvalid2-yhat2)/length(yvalid2);
legend(['fit = ' num2str(fit) '% mse = ' num2str(mse)]);
title('SVF fit')
figure;
x = xcorr(err,err);
x = x/x(N);
stem(-40:1:40,x(N-40:N+40));
line([-40 40],[2.17/sqrt(N) 2.17/sqrt(N)],'color','r');
line([-40 40],[-2.17/sqrt(N) -2.17/sqrt(N)],'color','r');
title('SVF error auto-correlation')


%% Identification of a continuous time model (ivsvf)

lambda = 1;
na = 6;
nb = 4;
nk = 1;
Fiv = ivsvf(uident,yident,dt,na,nb,Fsvf,nk,lambda);
yhat = lsim(Fiv,uvalid,0:dt:(length(uvalid)-1)*dt);
yvalid2 = yvalid(100:end); %forget first 100 samples to absorb initial condition mismatch
yhat2 = yhat(100:end);
err = yhat2-yvalid2;
N = length(yvalid2);
figure;
plot([yvalid2 yhat2]);
fit = 100*(1-norm(yvalid2-yhat2)/norm(yvalid2-mean(yvalid2)));
mse = norm(yvalid2-yhat2)/length(yvalid2);
legend(['fit = ' num2str(fit) '% mse = ' num2str(mse)]);
title('IV fit');
figure;
x = xcorr(err,err);
x = x/x(N);
stem(-40:1:40,x(N-40:N+40));
line([-40 40],[2.17/sqrt(N) 2.17/sqrt(N)],'color','r');
line([-40 40],[-2.17/sqrt(N) -2.17/sqrt(N)],'color','r');
title('IV error auto-correlation');


%% Step responses comparisons

y1 = lsim([Frls;Fels],ustep,tstep);
y2 = lsim([Fsvf;Fiv],ustep,tstep);
figure;
plot(tstep,ystep,tstep,y1,tstep,y2);
legend('real','discrete: rls','discrete: els','continuous: svf', 'continuous: iv');
title('Step response comparison');

%% Algorithms

%% Least Square

function [A,B,yhat,sigmahat,covtheta] = LS(uid,yid,na,nb,nk)
uid = reshape(uid,[],1);
yid = reshape(yid,[],1);

n = max(na,nb);
Y = yid(n+1:end);
PHI = [];
for i = 1:na
    PHI = [PHI -yid(n-i+1:end-i)];
end
for i = 1:nb
    PHI = [PHI uid(n-i+1:end-i)];
end
theta = PHI\Y;
A = [1;theta(1:na)]';
B = [zeros(nk+1,1);theta(na+1:na+nb)]';
yhat = PHI*theta;
sigmahat = (Y-yhat)'*(Y-yhat)/(length(yhat)-na-nb);
covtheta = sigmahat^2*inv(PHI'*PHI);
end

%% Recursive Least Square

function [A,B] = RLS(uid,yid,na,nb,nk)
uid = reshape(uid,[],1);
yid = reshape(yid,[],1);
lambda = 1; %forgetting factor
%initialization
theta = zeros(na+nb,1);
D = 10^8*eye(na+nb);
%Recursive algorithm
for t = max([na nb+nk]):(length(yid)-1)
    phi = [-yid(t:-1:t-na+1);uid(t-nk:-1:t-nb-nk+1)];
    yhatprior = theta'*phi;
    eps0 = yid(t+1)-yhatprior;
    epsposterior = eps0/(1+phi'*D*phi);
    theta = theta + D*phi*epsposterior;
%     D = D - D*phi*phi'*D/(1+phi'*D*phi);
    D = (D - D*phi*phi'*D/(lambda+phi'*D*phi))/lambda; %including forgetting factor
end
A = [1;theta(1:na)]';
B = [zeros(nk+1,1);theta(na+1:na+nb)]';
end

%% Recursive Extended Least Square

function [A,B,C] = ELS(uid,yid,na,nb,nc,nk)
uid = reshape(uid,[],1);
yid = reshape(yid,[],1);
lambda = 1;
%initialization
theta = zeros(na+nb+nc,1);
D = 10^6*eye(na+nb+nc);
epsposterior = zeros(length(yid),1);
eps0 = zeros(length(yid),1);
%Recursive algorithm
for t = max([na nb+nk nc]):(length(yid)-1)
    phi = [-yid(t:-1:t-na+1);
        uid(t-nk:-1:t-nb-nk+1);
        epsposterior(t:-1:t-nc+1)];
    yhatprior = theta'*phi;
    eps0(t+1) = yid(t+1)-yhatprior;
    epsposterior(t+1) = eps0(t+1)/(1+phi'*D*phi);
    theta = theta + D*phi*epsposterior(t+1);
    D = (D - D*phi*phi'*D/(lambda+phi'*D*phi))/lambda;
end

A = [1;theta(1:na)]';
B = [zeros(nk+1,1);theta(na+1:na+nb)]';
C = theta(na+nb+1:end)';
end

%% Prediction Error and AIC criteria

function [err, J, AIC,yhat] = computerrors(uvalid,yvalid,A,B,C)
N = length(uvalid);
if nargin == 4
    C = [];
end
uvalid = reshape(uvalid,[],1);
yvalid = reshape(yvalid,[],1);
A = A(2:end);
B = B(2:end);
na=length(A);
nb=length(B);
nc=length(C);
theta = [A';B';C'];
err = zeros(N,1);
yhat = zeros(N,1);
%Recursive algorithm
for t = max([na nb nc]):(length(yvalid)-1)
    phi = [-yvalid(t:-1:t-na+1);uvalid(t:-1:t-nb+1);err(t:-1:t-nc+1)];
    yhat(t+1) = theta'*phi;
    err(t+1) = yvalid(t+1)-yhat(t+1);
end
J = 1/N*(err'*err);
AIC = log(J)+2*(na+nb+nc)/N;
end

%% Identification of continuous time model

function F = svf(uid,yid,Ts,na,nb,nk,lambda)
% Identification of a continuous-time model from sampled data
% uid: input data
% yid: output data
% Ts: sampling time (s)
% na: order of the denominator
% nb: order of the denominator
% nk: delay in samples (if known), default 0
% lambda: time constant of the svf, default lambda = 5*Ts
% Algorithm based on "Identification of Continuous-time Models from Sampled Data
% Garnier, Wang. Springer
narginchk(5,7);
switch nargin
    case 6
        lambda = Ts*5;
    case 5
        lambda = Ts*5;
        nk = 0;
end
if isempty(nk), nk = 0; end
if isempty(lambda), lambda = Ts*5;end
if nb>na error('nb must be less or equal to na (proper transfer function)');end

N = length(uid);
uid = reshape(uid,[N 1]);
yid = reshape(yid,[N 1]);

T = 0:Ts:(N-1)*Ts;
s = tf('s');
f = 1/(s+lambda)^(na+1);
PHIu = lsim(f,uid,T);
for k = 1:nb
    f = s*f;
    PHIu = [lsim(f,uid,T) PHIu];
end
s = tf('s');
f = 1/(s+lambda)^(na+1);
PHIy = -lsim(f,yid,T);
for k = 1:na-1
    f = s*f;
    PHIy = [-lsim(f,yid,T) PHIy];
end
PHI = [PHIy PHIu];
f = s*f;
Y = lsim(f,yid,T);
theta = pinv(PHI)*Y;
A = [1 theta(1:na)'];
B = theta(na+1:end)';
F = tf(B,A);
end

function err = svf_errors(uvalid,yvalid,Ts,F,na,nb,lambda)
% Identification of a continuous-time model from sampled data
N = length(uvalid);
uvalid = reshape(uvalid,[N 1]);
yvalid = reshape(yvalid,[N 1]);
T = (0:Ts:(N-1)*Ts)';
yhat = lsim(F,uvalid,T);
err = yvalid(100:end) - yhat(100:end);
plot([yvalid yhat])
% 
end

%%
function F = ivsvf(uid,yid,Ts,na,nb,F0,nk,lambda)
% Identification of a continuous-time model from sampled data
% uid: input data
% yid: output data
% Ts: sampling time (s)
% na: order of the denominator
% nb: order of the denominator
% nk: delay in samples (if known), default 0
% lambda: time constant of the svf, default lambda = 5*Ts
% F0: initial guess for the model (lti model)
% Algorithm based on "Identification of Continuous-time Models from Sampled Data
% Garnier, Wang. Springer
narginchk(5,8);
switch nargin
    case 7
        lambda = Ts*5;
    case 6
        lambda = Ts*5;
        nk = 0;
end
if isempty(nk), nk = 0; end
if isempty(lambda), lambda = Ts*5;end
if nb>na error('nb must be less or equal to na (proper transfer function)');end

N = length(uid);
uid = reshape(uid,[N 1]);
yid = reshape(yid,[N 1]);

T = 0:Ts:(N-1)*Ts;
s = tf('s');
f = 1/(s+lambda)^(na+1);
PHIu = lsim(f,uid,T);
for k = 1:nb
    f = s*f;
    PHIu = [lsim(f,uid,T) PHIu];
end
s = tf('s');
f = 1/(s+lambda)^(na+1);
PHIy = -lsim(f,yid,T);
for k = 1:na-1
    f = s*f;
    PHIy = [-lsim(f,yid,T) PHIy];
end
z = lsim(F0,uid,T);
s = tf('s');
f = 1/(s+lambda)^(na+1);
PHIz = -lsim(f,z,T);
for k = 1:na-1
    f = s*f;
    PHIz = [-lsim(f,z,T) PHIz];
end

PHI = [PHIy PHIu];
Z = [PHIz PHIu];
f = s*f;
Y = lsim(f,yid,T);
theta = (Z'*PHI)\Z'*Y; % (better than (Z'*PHIu)^-1*Z'*Y)
A = [1 theta(1:na)'];
B = theta(na+1:end)';
F = tf(B,A);
end