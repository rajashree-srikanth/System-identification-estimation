%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D orientation estimation template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONFIGURATION and INITIALISATION
dt=0.01;
COM='COM7';
useSimulation = false;
useSimulation_visuals = true; %useSimulation;
%% CALIBRATION
% To re-run calibration, just clear MAG0 or run the calibration script
% separately:
if ~exist('MAG0', 'var')
    runcalibration;
end
%% affichage du parallélépipède
launch3D;

%% estimation attitude
qa = 0.0001; %state orientation noise
ra = 100*ACCVAR(2:3); %Acceleration measurement noise
rg = GYRVAR(1); %Gyro measurement noise
Q = diag([qa]); %state noise matrix - denotes the trust in system model
R = diag([ra]); %output noise matrix - denotes the trust in sesnor model
X = [0]'; %initial state theta
P = deg2rad(10)^2; %inital variance 

tp=0;
ii=0;
obs = [];
xtrue = [];
xhat = [];
if useSimulation
    imu411('open',COM, 'simimu_2Dmaneuver');
else
    imu411('open',COM);
end
while(1)
    ii=ii+1;
    % Read sensor
    [d, xt] = imu411('read');
    %cette lecture est bloquante et cadence la boucle a 0.01s
    obs(:, end+1) = d;
    xtrue(:, end+1) = xt;
    % Rappel :
    % d(1)    = Time                (s)
    % d(2:4)  = Gyroscope     X,Y,Z (°/s)
    % d(5:7)  = Accelerometer X,Y,Z (g)
    
    t=d(1);
    % Predict
    X = X;
    F = [1];
    P = F*P*F' + Q;
    % Update
    if ~isnan(t)
        Y = [d(7) d(6)]';
        Yhat = [cos(X) -sin(X)]';
        V = Y - Yhat;
        H = [-sin(X) -cos(X)]';
        S = H*P*H' + R;
        K = P*H'*inv(S);
        X = X + K*V;
        P = P - K*S*K';
    end
    
    % Update Visualisation:
    xhat(:, end+1) = [angle2quat(0, 0, -X)'; zeros(3,1)];
    DCM_k = angle2dcm(0, 0, X, 'ZYX')';
    update3D;
end