clc, clear all, close all
loadImage=imread('loadImage.jpg');
loadImage= rgb2gray(loadImage);
load('measuredData.mat');


tic
no_of_sensors = 2; % choose a number whose square root is an integer
% Particle filter parameters
M = 20; % samples of signal recieved
beta = 40;
snr = 10^4;
mu = beta / snr; % threshold
endK = 169; % final time index
niter = 200; % number of Monte Carlo iterations
% Location of the sensors
xy_max = 250;
x_loc = xy_max * rand(1,no_of_sensors);
y_loc = xy_max * rand(1,no_of_sensors);

dt = 1;
F=1;
qtrue = .1;
Q = qtrue*dt;
  
% Generate state and observation sequences
XTrue(:,1) = [x(1); y(1)]; % Initial state
for k=1:endK
    XTrue(:,k+1) = XTrue(:,k)+dt*[vx(k);vy(k)] + Q*randn(2,1);
    distance_vector(:,k) = distance( XTrue(1,k+1),XTrue(2,k+1),x_loc,y_loc)';
    Z(:,k) = Z_k(M,distance_vector(:,k),mu)';

end

N=100; % Total number of particles 
ErrSum = zeros(2,endK+1);
SqErrSum = zeros(2,endK+1);
red_ErrSum = zeros(2,endK+1);
red_SqErrSum = zeros(2,endK+1);
% Parameters for keeping track of reduced set of sensors
C_r_k = 1;
r_k = 500;
avg_sensors = zeros(1,niter);


XTrue = squeeze( XTrue ); % Retreiving values generated
Z = squeeze( Z); % in data generation script

% Particle filter implementation
% Generate the initial particles
init_vector = XTrue(:,1);
X = repmat(init_vector,1,N) + ...
sqrtm([10 0 ; 0 10]) * randn(2,N);

red_X = X;
w = ones(1,N)/N;
red_w = w;
xhat(:,1) = X*w';
red_xhat(:,1) = xhat(:,1);
for k=1:endK
    V=repmat([vx(k);vy(k)],1,N) + ...
    sqrtm([1 0 ; 0 1]) * randn(2,N);
    rand_matrix = randn(2,N);
    xExt = F*X +dt*V+ Q*rand_matrix;
    red_xExt = F*red_X + Q*rand_matrix;
    for n = 1:N
        dist = distance( xExt(1,n),xExt(2,n),x_loc,y_loc)';
        P_det = P_d(dist,M,beta,snr);
        exponent = Z(1,k);
        P_Zk_Xk = P_det(1)^exponent * ( 1 - P_det(1) )^(1-exponent);
        for l=2:no_of_sensors
            exponent = Z(l,k);
            P_Zk_Xk = P_Zk_Xk * ( P_det(l)^exponent * ( 1 - P_det(l) )^(1-exponent) );
        end
        w(n) = P_Zk_Xk;
    end
    if sum(w) > 0
    w = w/sum(w);
    else
    w = ones(1,N)/N;
    end
    if sum(red_w) > 0
    red_w = red_w/sum(red_w);
    else
    red_w = ones(1,N)/N;
    end
    % Compute the estimate of the state
    xhat(:,k+1) = xExt * w';
    red_xhat(:,k+1) = red_xExt * red_w';
    % Finding the covariance matrix
    P_k_k = red_xExt - repmat( red_xhat(:,k+1),1,N );
    P_k_k = ( repmat(red_w,2,1) .* P_k_k ) * P_k_k';
    e_k = P_k_k(1,1) + P_k_k(2,2);
    r_k = C_r_k * e_k;

    % Resample the particles
    len = length(w);
    % Cumulative Distributive Function
    cumpr = cumsum(w(1,1:len))';
    cumpr = cumpr/max(cumpr);
    red_cumpr = cumsum(red_w(1,1:len))';
    red_cumpr = red_cumpr/max(red_cumpr);
    rand_no = rand(1,1);
    u(1,1) = (1/N)*rand_no;
    a=1;
    for b = 1:N
    u(b,1)= u(1,1) + (1/N)*(b-1);
    while (u(b,1) > cumpr(a,1))
    a = a+1;
    if a > N
    break
    end
    end
    if a <= N
    x_update(:,b) = xExt(:,a);
    end
    end
    X = x_update;
    w = 1/N*ones(1,N);
    red_u(1,1) = (1/N)*rand_no;
    red_a=1;
    for red_b = 1:N
        red_u(red_b,1)= red_u(1,1) + (1/N)*(red_b-1);
        while (red_u(red_b,1) > red_cumpr(red_a,1))
            red_a = red_a+1;
            if red_a > N
            break
            end
        end
        if red_a <= N
            red_x_update(:,red_b) = red_xExt(:,red_a);
        end
    end
    red_X = red_x_update;
    red_w = w;
    figure(1)
    imshow(loadImage)
    hold on 
    plot(X(1,:),286-X(2,:),'+')
    hold off
    pause(0.2)
end

Err=XTrue-xhat;
SqErr=(XTrue-xhat).^2;


p = SqErrSum-ErrSum.^2;
figure(2)
plot(0:endK,XTrue(1,:),'b',0:endK,xhat(1,:),'r');
xlabel('time step'); ylabel('car x position');
legend('True value', 'predicted value');

figure(3)
plot(0:endK,XTrue(2,:),'b',0:endK,xhat(2,:),'r');
xlabel('time step'); ylabel('car y position');
legend('True value', 'predicted value');



