%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D orientation estimation template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONFIGURATION and INITIALISATION
dt=0.01;
COM='COM7';
useSimulation = true;
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
qa = 0.0001; %Bruit d'état orientation
ra = 100*ACCVAR(2:3); %Bruit de mesure accéléro
rg = GYRVAR(1); %Bruit de mesure gyro
Q = diag([qa]);
R = diag([ra]);
X=[0]'; %Etat : rotation selon x (rad)
P = deg2rad(10)^2;

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
    [d, xt] = imu411('read'); %cette lecture est bloquante et cadence la boucle a 0.01s
    obs(:, end+1) = d;
    xtrue(:, end+1) = xt;
    % Rappel :
    % d(1)    = Time                (s)
    % d(2:4)  = Gyroscope     X,Y,Z (°/s)
    % d(5:7)  = Accelerometer X,Y,Z (g)
    
    t=d(1);
    % Predict
    X_temp = [X(1) X(2)] ;
    F = eye(2);
    F(1,2) = dt;
    P = F*P*F' + Q;

    % Use gyro measurement as control input
    u = [0; d(2)*pi/180]; % Assuming the gyro measurement is in degrees per second
    X_temp = X_temp + u*dt;

    % Update
    if ~isnan(t)
        Y = [d(6); d(7)];
        Yhat = [-sin(X_temp(1)); cos(X_temp(1))];
        v = Y - Yhat;
        H = [-cos(X_temp(1)) 0; -sin(X_temp(1)) 0];
        S = H*P*H' + R;
        K = P*H'*pinv(S);
        X_temp = X_temp + K*v;
        P = P - K*S*K'
    end
    
    % Update Visualisation:
    xhat(:, end+1) = [angle2quat(0, 0, -X)'; zeros(3,1)];
    DCM_k = angle2dcm(0, 0, X, 'ZYX')';
    update3D;

end