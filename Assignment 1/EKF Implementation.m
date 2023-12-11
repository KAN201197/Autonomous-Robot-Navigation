%% Clear past plots, variables and commands
close all; clear all; clc;

% Load data 
load 'EE5114_CA1.mat';

% acx = x-axis accelerometer reading
% acy = y-axis accelerometer reading
% acz = z-axis accelerometer reading
% 
% phi = Roll angle computed by the drone's on-board computer
% tht = Pitch angle computed by the drone's on-board computer
% psi = Yaw angle computed by the drone's on-board computer 
% 
% fix = GPS position fix signal 
% eph = GPS horizontal variance 
% epv = GPS vertical variance 
% lat = GPS Latitude
% lon = GPS Longitude
% alt = GPS altitude
% gps_nSat = Number of GPS satellites
% 
% out1 = Motor 1 signal
% out2 = Motor 2 signal
% out3 = Motor 3 signal
% out4 = Motor 4 signal

%% Accelerometer plot
figure; set(gcf,'numbertitle','off','name','Acceleration');  
subplot(3,1,1); plot(t, acx, 'b'); ylim([-2 2]); ylabel('acx (m/s^2)'); grid on; 
subplot(3,1,2); plot(t, acy, 'b'); ylim([-2 2]); ylabel('acy (m/s^2)'); grid on; 
subplot(3,1,3); plot(t, acz, 'b'); ylabel('acz (m/s^2)'); grid on; 

%% Euler angles plot
figure; set(gcf,'numbertitle','off','name','Euler Angles');  
subplot(3,1,1); plot(t, rad2deg(phi), 'b'); ylabel('Roll (degree)'); grid on; 
subplot(3,1,2); plot(t, rad2deg(tht), 'b'); ylabel('Pitch (degree)'); grid on; 
subplot(3,1,3); plot(t, rad2deg(psi), 'b'); ylabel('Yaw (degree)'); grid on; 

%% GPS plot
figure; set(gcf,'numbertitle','off','name','GPS');  
subplot(3,2,1); plot(t, lon); ylabel('Longitude'); grid on;
subplot(3,2,3); plot(t, lat); ylabel('Latitude'); grid on;
subplot(3,2,5); plot(t, alt); ylabel('Altitude'); grid on; xlabel('time (s)');

subplot(3,2,2); plot(t, gps_nSat, '.'); ylabel('Sat'); grid on;
subplot(3,2,4); plot(t, eph); ylabel('Eph'); grid on; ylim([0 5]);
subplot(3,2,6); plot(t, epv); ylabel('Epv'); grid on; ylim([0 5]);

%% Motor signal plot
figure; set(gcf,'numbertitle','off','name','Motor Signal');  
hold on;
plot(t,out1,'r');
plot(t,out2,'g');
plot(t,out3,'b');
plot(t,out4,'y');
legend('Motor1','Motor2','Motor3','Motor4'); 
ylabel('Motor inputs'); xlabel('time (s)'); ylim([1000 2000]); grid on;


%%%%%%%%%%%%%%%%%%%%%% Your own coding work start from here %%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert GPS raw measurements to local NED position values

% Define parameter to calculate prime vertical radius of curvature
a = 6378137;
b = 6356752.3142;
e2 = 1 - (b^2/a^2);
N = (a./sqrt(1 - e2 * sin(deg2rad(lat)).^2));

% Define coordinates reference point of GPS data
lat_ref = lat(1); alt_ref = alt(1); lon_ref = lon(1);

% Conversion from GPS coordinate into ECEF coordinates
xe = (N + alt) .* cos(deg2rad(lat)) .* cos(deg2rad(lon));
ye = (N + alt) .* cos(deg2rad(lat)) .* sin(deg2rad(lon));
ze = ((b^2/a^2)*N + alt) .* sin(deg2rad(lat));

% Conversion reference point from GPS data into ECEF coordinates
xe_ref = (N(1) + alt_ref) * cos(deg2rad(lat_ref)) * cos(deg2rad(lon_ref));
ye_ref = (N(1) + alt_ref) * cos(deg2rad(lat_ref)) * sin(deg2rad(lon_ref));
ze_ref = ((b^2/a^2)*N(1) + alt_ref) * sin(deg2rad(lat_ref));

% Conversion from ECEF into NED coordinates
xn = zeros(size(xe));
yn = zeros(size(ye));
zn = zeros(size(ze));

for i = 1:length(xe)
    lat_i = deg2rad(lat(i));
    lon_i = deg2rad(lon(i));
    R = [-sin(lat_i) * cos(lon_i), -sin(lon_i), -cos(lat_i) * cos(lon_i);
         -sin(lat_i) * sin(lon_i), cos(lon_i), -cos(lat_i) * sin(lon_i);
         cos(lat_i), 0, -sin(lat_i)];
    ned_dif = [xe(i); ye(i); ze(i)] - [xe_ref; ye_ref; ze_ref];
    ned_i = R' * ned_dif;
    xn(i) = ned_i(1);
    yn(i) = ned_i(2);
    zn(i) = ned_i(3);
end

%% Implement EKF to estimate NED position and velocity

% Specify time range for EKF
t_min = 1445 - 600; % 10 minutes before second take off
t_max = 1690 + 300; % 5 minutes after second landing 

% Find indices for specified time range
t_min_index = find(t >= t_min, 1);
t_max_index = find(t >= t_max, 1);
t_range = t(t_min_index:t_max_index);

% Extended Kalman Filter (EKF) Implementation

% Define measurement noise covariances
R = diag([0.1, 0.1, 0.1]); % measurement noise covariance

% Initialize state and covariance matrix
x = [xn(t_min_index); yn(t_min_index); zn(t_min_index);...
    0; 0; 0; ...
    0; 0; 0]; % [NED pos x, y, z, NED vel x, y, z, acc bias x, y, z]
P = diag([1, 1, 1, 1, 1, 1, 1, 1, 1]); % covariance matrix

% EKF loop
ned_estimates = zeros(9, length(t_range)); % this variable to store the value for each loop of EKF implementation

for i = t_min_index:t_max_index

    dt = t(i+1) - t(i);
    y = [xn(i+1); yn(i+1); zn(i+1)]; % GPS measurement model

    % Prediction step
    R_btog = [cos(psi(i))*cos(tht(i)), (cos(psi(i))*sin(tht(i))*sin(phi(i)))-(sin(psi(i))*cos(phi(i))), (cos(psi(i))*sin(tht(i))*cos(phi(i)))+(sin(psi(i))*sin(phi(i)));
              sin(psi(i))*cos(tht(i)), (sin(psi(i))*sin(tht(i))*sin(phi(i)))+(cos(psi(i))*cos(phi(i))), (sin(psi(i))*sin(tht(i))*cos(phi(i)))-(cos(psi(i))*sin(phi(i)));
              -sin(tht(i)), cos(tht(i))*sin(phi(i)), cos(tht(i))*cos(phi(i))]; % rotation matrix from body frame into NED frame
    F_upper = [eye(3), diag([dt, dt, dt]), diag([-(dt^2)/2,-(dt^2)/2,-(dt^2)/2]); ...
               zeros(3,3), eye(3), diag([-dt,-dt,-dt])] * [eye(6), zeros(6,3); zeros(3,6), R_btog];
    F = [F_upper; [zeros(3,3), zeros(3,3), eye(3)]]; % state transition matrix
    G = [diag([(dt^2)/2,(dt^2)/2,(dt^2)/2]); diag([dt,dt,dt]); zeros(3,3)] * [R_btog,[0;0;9.8]]; % control input matrix
    Q = G * diag([10, 10, 10, 0]) * G'; % Process noise covariance matrix

    x_predicted = F * x + G * [acx(i); acy(i); acz(i); 1];
    P_predicted = F * P * F' + Q;

    if rem(i,5) == 0 % Condition whether do correction step when there is a new measurement data from GPS
        % Correction step
        H = [eye(3), zeros(3,3), zeros(3,3)]; % Measurement Matrix
        K = P_predicted * H' * inv(H * P_predicted * H' + R); % Kalman Gain

        % Update Step using new GPS data measurement
        x = x_predicted + K * (y - (H * x_predicted));
        P = (eye(9) - K * H) * P_predicted;
    else
        % Update step without using GPS data measurement
        x = x_predicted;
        P = P_predicted;
    end

    ned_estimates(:,i) = x;
end

ned_estimates = ned_estimates(:,t_min_index:t_max_index); % Retrieving estimated data from t_min index into t_max index
%% Result plots

% NED coordinates Plot againt time
figure; set(gcf,'numbertitle','off','name','NED positions');  
subplot(3,1,1); plot(t, xn, 'b'); ylabel('xn (m)'); grid on; 
subplot(3,1,2); plot(t, yn, 'b'); ylabel('yn (m)'); grid on; 
subplot(3,1,3); plot(t, zn, 'b'); ylabel('zn (m)'); grid on; 
xlabel('Time (s)');

% Estimated NED Position based on EKF Implementation againt time
figure; set(gcf,'numbertitle','off','name','Estimated NED positions');  
subplot(3,1,1); plot(t_range, xn(t_min_index:t_max_index), 'b', t_range, ned_estimates(1, :), 'r'); ylabel('NED x (m)'); legend('Ground Truth', 'Estimated'); grid on;
subplot(3,1,2); plot(t_range, yn(t_min_index:t_max_index), 'b', t_range, ned_estimates(2, :), 'r'); ylabel('NED y (m)'); legend('Ground Truth', 'Estimated'); grid on;
subplot(3,1,3); plot(t_range, zn(t_min_index:t_max_index), 'b', t_range, ned_estimates(3, :), 'r'); ylabel('NED z (m)'); legend('Ground Truth', 'Estimated'); grid on;
xlabel('Time (s)');

% Estimated NED velocities based on EKF Implementation againt time

% Calculate the difference in NED and time to determine NED velocity
ned_x_diff = diff(xn(t_min_index:t_max_index));
ned_y_diff = diff(yn(t_min_index:t_max_index));
ned_z_diff = diff(zn(t_min_index:t_max_index));
t_diff = diff(t_range);

% plot estimated NED velocity vs ground truth NED velocity
figure; set(gcf,'numbertitle','off','name','Estimated NED velocities');
subplot(3,1,1); plot(t_range(2:end), ned_x_diff./t_diff, 'b', t_range, ned_estimates(4, :), 'r'); ylabel('NED v_x (m/s)'); legend('Ground Truth', 'Estimated'); grid on;
subplot(3,1,2); plot(t_range(2:end), ned_y_diff./t_diff, 'b', t_range, ned_estimates(5, :), 'r'); ylabel('NED v_y (m/s)'); legend('Ground Truth', 'Estimated'); grid on;
subplot(3,1,3); plot(t_range(2:end), ned_z_diff./t_diff, 'b', t_range, ned_estimates(6, :), 'r'); ylabel('NED v_z (m/s)'); legend('Ground Truth', 'Estimated'); grid on;
xlabel('Time (s)');

% plot accelerometer bias
figure; set(gcf,'numbertitle','off','name','Accelerometer Bias');  
subplot(3,1,1); plot(t_range, ned_estimates(7, :), 'r'); ylabel('ACC_BIAS x (m/s^2)'); grid on;
subplot(3,1,2); plot(t_range, ned_estimates(8, :), 'r'); ylabel('ACC_BIAS y (m/s^2)'); grid on;
subplot(3,1,3); plot(t_range, ned_estimates(9, :), 'r'); ylabel('ACC_BIAS z (m/s^2)'); grid on;
xlabel('Time (s)');
