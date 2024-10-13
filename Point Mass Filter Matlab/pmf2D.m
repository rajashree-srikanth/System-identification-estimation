%Point Mass Filter example
clear

%% Parameters
Tmax = 20; %Number of time steps
Q = [0.001 0
    0 0.001]; %State noise
R = [0.04 0
    0 0.04]; %Measurement noise

I = 30;     %size of the grid along x
J = 30;     %size of the grid along y
N = I*J;    %Number of points in the grid
minx = -2;
maxx = 2;
miny = -2;
maxy = 2;

%initial position
xtrue = [1
    0];

%Initial guess
P = [1 0
    0 1];             %Initial variance of estimated position
M = xtrue+[-0.5 0.5]';  %Initial Guess for the position

W = zeros(I,J);
for i = 1:I
    for j = 1:J
        x = [i*(maxx-minx)/(I)+minx
            j*(maxy-miny)/(J)+miny];
        W(i,j) = 1/sqrt((2*pi)^2*det(P))*exp(-1/2*(x-M)'*P^-1*(x-M));
    end
end
W = W/max(W,[],'all'); %Normalize
[X,Y] = meshgrid(minx:(maxx-minx)/(I-1):maxx,miny:(maxy-miny)/(J-1):maxy);

%Plots initialization
figure(1);axis equal;
axis([minx maxx miny maxy]);
hx = line(xtrue(1),xtrue(2),1,'Marker','o','MarkerFaceColor','k');
hy = line(xtrue(1),xtrue(2),1,'Marker','x','MarkerSize',12);
hold on;
hz = scatter(reshape(Y,[N 1]),reshape(X,[N 1]),40,[1 0 0],'filled');
hz.AlphaData = reshape(W,[N 1]);
hz.MarkerFaceAlpha = 'flat';


%% Loop
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
    title(['step ' num2str(t)]);
    pause;
    title(['step ' num2str(t) ' : Predicting']);
    pause(0.01);
    %Predict
    for i = 1:I
        for j = 1:J
            x = [i*(maxx-minx)/(I)+minx j*(maxy-miny)/(J)+miny]';
            W(i,j) = 0;
            for ii = 1:I
                for jj = 1:J
                    xx = [ii*(maxx-minx)/(I)+minx jj*(maxy-miny)/(J)+miny]';
                    % probability to jump from (ii,jj) to (i,j)
                    dd = x - xx - [-2*pi/Tmax*sin(2*pi*t/Tmax) 2*pi/Tmax*cos(2*pi*t/Tmax)]';
                    pp = exp(-1/2*(dd)'*(Q)^-1*(dd));
                    W(i,j) = W(i,j) + W(ii,jj)*pp;
                end
            end
        end
    end
    W = W/max(W,[],'all'); %Normalize
    %update plot
    hz.AlphaData = reshape(W,[N 1]);
    hz.MarkerFaceAlpha = 'flat';
    title(['step ' num2str(t) ' ; After Prediction']);
    pause;

    %update
    for i = 1:I
        for j = 1:J
            x = [i*(maxx-minx)/(I)+minx j*(maxy-miny)/(J)+miny]';
            W(i,j) = W(i,j)*exp(-1/2*(x-ytrue)'*(R/1)^-1*(x-ytrue));
        end
    end
    W = W/max(W,[],'all'); %Normalize

    %update plot
    hz.AlphaData = reshape(W,[N 1]);
    hz.MarkerFaceAlpha = 'flat';
    title(['step ' num2str(t) ' ; After Update']);
    pause;
end