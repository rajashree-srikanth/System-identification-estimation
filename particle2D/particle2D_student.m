% %Particle Filter example
clear 
%% Parameters
Tmax = 20; %Number of time steps
Q = [0.001 0
    0 0.001]; %State noise
R = [0.04 0
    0 0.04]; %Measurement noise

N = 100000;    %Number of particles
minx = -2;
maxx = 2;
miny = -2;
maxy = 2;

%initial position
xtrue = [1
    0];

%Initial guess
P = [0.02 0
    0 0.02];             %Initial variance of estimated position
M = xtrue;   %Initial Guess for the position

Wi = 1/N*ones(1,N); %initial weight of the particles
xi = M + chol(P')*randn(2,N); %Initial particles


%Plots initialization
figure(1); %Particle, real positions, measurement, ellipsoid of confidence
axis equal;
axis([minx maxx miny maxy]);
hx = line(xtrue(1),xtrue(2),1,'Marker','o','MarkerFaceColor','k');
hy = line(xtrue(1),xtrue(2),1,'Marker','x','MarkerSize',12);
hold on;
hz = scatter(xi(1,:),xi(2,:),40,[1 0 0]);
hz.MarkerFaceColor = 'r';
hz.MarkerFaceAlpha = 'flat';
hz.AlphaDataMapping = 'none';
hz.AlphaData = Wi/max(Wi);
%Create and plot Ellipsoid
X0=sum(xi.*[Wi;Wi],2);
M=(xi.*Wi-X0*ones(1,N).*Wi)*(xi-X0*ones(1,N))';
[U,S,V]=svd(M);
a=0:0.1:2*pi;
ell0=[sqrt(S(1,1))*cos(a);sqrt(S(2,2))*sin(a)]*3; %3-sigma interval,
                    % corresponds to 98.8% of confidence interval (well explained in J.Sola Phd Thesis)
ell=X0*ones(1,length(a))+V'*ell0;
hxell = line(ell(1,:),ell(2,:),'LineWidth',2,'Color','r');
figure(2);  %Particles weights (sorted)
hstem = stem(sort(Wi));


%% Loop
variance = ones(1,Tmax);
rms_arr = ones(1,Tmax);
for t = 1:Tmax
    % True position
    [v] = chol(Q)*randn(2,1);
    xtrue = xtrue...
        + [-2*pi/Tmax*sin(2*pi*t/Tmax) 2*pi/Tmax*cos(2*pi*t/Tmax)]'...
        +[v(1) v(2)]'; %trajectory is a circle plus perturbation

    %Measurement
    w = randn(1,2)*chol(R);
    ytrue = xtrue + [w(1) w(2)]'; %Measurement is true position plus noise

    %update plot
    hx.XData = xtrue(1);
    hx.YData = xtrue(2);
    hy.XData = ytrue(1);
    hy.YData = ytrue(2);
    figure(1);title(['step ' num2str(t)]);
    %pause;
    title(['step ' num2str(t) ' : Predicting']);
    %pause(0.01);

    %Predict
   
        v = chol(Q)*randn(2,N);
         xi = xi + [-2*pi/Tmax*sin(2*pi*t/Tmax); 2*pi/Tmax*cos(2*pi*t/Tmax)] + v;
   
         
    %update particle positions after prediction
    hz.XData = xi(1,:);
    hz.YData = xi(2,:);
    %Update Ellipsoid after prediction
    X0=sum(xi.*[Wi;Wi],2);
    M=(xi.*Wi-X0*ones(1,N).*Wi)*(xi.*Wi-X0*ones(1,N).*Wi)'*N;
    [U,S,V]=svd(M);
    a=0:0.1:2*pi;
    ell0=[sqrt(S(1,1))*cos(a);sqrt(S(2,2))*sin(a)]*sqrt(3);
    ell=X0*ones(1,length(a))+V'*ell0;
    set(hxell,'Xdata',ell(1,:),'Ydata',ell(2,:));
    figure(1);title(['step ' num2str(t) ' ; After Prediction']);
    %pause;

    %update particles weights
    for i = 1:N
        %Wi(i) = 1/N;
        particle_measurement = xi(:, i);
        likelihood = 1 /(2 * pi * sqrt(det(R))) * exp(-0.5 * (ytrue - particle_measurement)' * inv(R) * (ytrue - particle_measurement));

        % Update the particle weight
        Wi(i) = Wi(i) * likelihood;
    end
    Wi = Wi/sum(Wi); %Normalize
   


    %update particle weights and ellipsoid of confidence
    hz.XData = xi(1,:);
    hz.YData = xi(2,:);
    hz.AlphaData = Wi/max(Wi);
    hz.MarkerFaceAlpha = 'flat';
    hstem.YData = sort(Wi);
    %Update Ellipsoid
    X0=sum(xi.*[Wi;Wi],2);
    M=(xi.*Wi-X0*ones(1,N).*Wi)*(xi.*Wi-X0*ones(1,N).*Wi)'*N;
    [U,S,V]=svd(M);
    a=0:0.1:2*pi;
    ell0=[sqrt(S(1,1))*cos(a);sqrt(S(2,2))*sin(a)]*sqrt(3);
    ell=X0*ones(1,length(a))+V'*ell0;
    set(hxell,'Xdata',ell(1,:),'Ydata',ell(2,:));
    figure(1);title(['step ' num2str(t) ' ; After Update']);
     %figure(2);title(['step ' num2str(t) ' ; After Update ; Neff = ' num2str(Neff)]);
    %pause;
   

     %Resampling
     Neff= 1/sum(Wi.^2); %resampling criteria
     if Neff <= 0.5*N
        [Wi,xi] = resample(Wi,xi);
     end

    a = xtrue-xi;
    b = sqrt(a(1,:).^2 + a(2,:).^2);
    rms_arr(t) = max((b));
    variance(t) = var((b));
end
maxvar = max(variance)
maxrms = max(rms_arr)
save('N100000_Neff_0.5_QR.mat',"b",'a', 'variance')

%% Resampling algorithm
function [W,X] = resample(W,X)
% Resampling algorithm
% Multinomial algorithm
% W : weights of the particles
% indx : index of the particles to be resampled
W = W/sum(W); % just in case the weights are not normalized
N = length(W);
Q = cumsum(W);
i=1;
indx = 1:N; % Preallocate indx (save time)
while (i<=N)
    sampl = rand(1,1);
    j=1;
    while (Q(j)<sampl)
        j=j+1;
    end
    indx(i)=j;
    i=i+1;
end
X = X(:,indx);
W = ones(1,N)/N;
end