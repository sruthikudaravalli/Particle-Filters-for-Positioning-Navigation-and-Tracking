clear;clc;close all;
load('height.mat')

no_of_sensors = 2; % choose a number whose square root is an integer
% Particle filter parameters
M = 20; % samples of signal recieved
beta = 40;
snr = 10^4;
mu = beta / snr; % threshold
endK = 25; % final time index
niter = 1; % number of Monte Carlo iterations
% Location of the sensors
xy_max = 600;
x_loc = xy_max * rand(1,no_of_sensors);
y_loc = xy_max * rand(1,no_of_sensors);
z_loc=0;
dt = 1;
% % F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
F = [1 0 0 dt 0 0 dt^2/2 0 0;...
     0 1 0 0 dt 0 0 dt^2/2 0;...
     0 0 1 0 0 dt 0 0 dt^2/2;...
     0 0 0 1 0 0 dt 0 0;...
     0 0 0 0 1 0 0 0 dt;...
     0 0 0 0 0 1 0 0 dt;...
     0 0 0 0 0 0 1 0 0;...
     0 0 0 0 0 0 0 1 0;...
     0 0 0 0 0 0 0 0 1;...
    ];
A=[dt^2/2 0 0;0 dt^2/2 0;0 0 dt^2/2;...
    dt 0 0;0 dt 0;0 0 dt;...
    0 0 0;0 0 0;0 0 0];
qtrue = .1;
Q = qtrue*[dt^3/6 0 0 dt^2/2 0 0 dt 0 0;...
    0 dt^3/6 0 0 dt^2/2 0 0 dt 0;...
    0 0 dt^3/6 0 0 dt^2/2 0 0 dt;...
    0 0 0 dt^2/2 0 0 dt 0 0;...
    0 0 0 0 dt^2/2 0 0 dt 0;...
    0 0 0 0 0 dt^2/2 0 0 dt;...
    0 0 0  0 0 0 dt 0 0;...
    0 0 0  0 0 0 0 dt 0;...
    0 0 0  0 0 0 0 0 dt;...
    ];
Q=qtrue*Q;
sqrtQ=Q;
% sqrtQ = sqrtm(Q);
iter = 0;
while (iter < niter)
    iter = iter+1;
    flag = 0;
    % Generate state and observation sequences

     XTrue(:,1) = [0;50;3000; 4;2;2; 1;2;1]; % Initial state
    
    for k=1:endK
        XTrue(:,k+1) = F*XTrue(:,k) + sqrtQ*randn(9,1);
               
        distance_vector(:,k) = distance( XTrue(1,k+1),XTrue(2,k+1),x_loc,y_loc)';
        Z(:,k) = Z_k(M,distance_vector(:,k),mu)';
       
    end
   
    saved_Z(iter,:,:) = Z;
    saved_XTrue(iter,:,:) = XTrue;
    
    
end

N=100; % Total number of particles (temporary)
ErrSum = zeros(9,endK+1);
SqErrSum = zeros(9,endK+1);
red_ErrSum = zeros(9,endK+1);
red_SqErrSum = zeros(9,endK+1);
% Parameters for keeping track of reduced set of sensors
C_r_k = 1;
r_k = 500;
avg_sensors = zeros(1,niter);
for iter = 1:niter
    XTrue = squeeze( saved_XTrue(iter,:,:) ); % Retreiving values generated
    Z = squeeze( saved_Z(iter,:,:) ); % in data generation script

    % Particle filter implementation
    % Generate the initial particles
    % init_vector = [200;300;15;8];
    init_vector = XTrue(:,1);
    X = repmat(init_vector,1,N) + ...
    0.1*sqrtm([1 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0;...
    0 0 0 1 0 0 0 0 0;0 0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0 0;...
    0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 1]) * randn(9,N);
    red_X = X;
    w = ones(1,N)/N;
    red_w = w;
    xhat(:,1) = X*w';
    red_xhat(:,1) = xhat(:,1);
    for k=1:endK
        rand_matrix = randn(9,N);
        xExt = F*X + sqrtQ*rand_matrix;
        red_xExt = F*red_X + sqrtQ*rand_matrix;
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
        P_k_k = ( repmat(red_w,9,1) .* P_k_k ) * P_k_k';
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
        contour(zz)
        hold on 
        plot(red_X(1,:)/5,red_X(2,:)/5,'+')
        hold off
        pause(0.2)
    end
    avg_sensors(1,iter) = avg_sensors(1,iter) / endK;
    ErrSum = ErrSum + (XTrue-xhat);
    SqErrSum = SqErrSum + (XTrue-xhat).^2;
    red_ErrSum = red_ErrSum + (XTrue-red_xhat);
    red_SqErrSum = red_SqErrSum + (XTrue-red_xhat).^2;
    
end
p = SqErrSum/niter-ErrSum.^2/niter^2;
red_p = red_SqErrSum/niter-red_ErrSum.^2/niter^2;

for ii=1:endK+1
    h_ture(ii)=XTrue(3,ii)-zz(fix(XTrue(1,ii)/5)+1,fix(XTrue(1,ii)/5)+1);
    h_hat(ii)=xhat(3,ii)-zz(fix(xhat(1,ii)/5)+1,fix(xhat(1,ii)/5)+1);
end
figure(2);
clf;
subplot(3,1,1);
plot(0:endK,XTrue(1,:),'b--',0:endK,xhat(1,:),'r-');
title('Estimated x position sequence');
legend('True','Estimated');
subplot(3,1,2);
plot(0:endK,XTrue(2,:),'b--',0:endK,xhat(2,:),'r-');
title('Estimated y position sequence');
legend('True','Estimated');
subplot(3,1,3);
plot(0:endK,h_ture,'b--',0:endK,h_hat,'r-');
title('Estimated z position sequence');
legend('True','Estimated');


figure(3);
clf;
subplot(3,1,1);
plot(0:endK,XTrue(4,:),'b--',0:endK,xhat(4,:),'r-');
title('Estimated x velocity sequence');
legend('True','Estimated');
subplot(3,1,2);
plot(0:endK,XTrue(5,:),'b--',0:endK,xhat(5,:),'r-');
title('Estimated y velocity sequence');
legend('True','Estimated');
subplot(3,1,3);
plot(0:endK,XTrue(6,:),'b--',0:endK,xhat(6,:),'r-');
title('Estimated z velocity sequence');
legend('True','Estimated');