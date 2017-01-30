function C_KALMAN_filter
    % C_KALMAN_filter.m
    clc;
    close all;
    clear;

    Fs = 100; % Sampling freq, Hz
    w_ref = 8; % rpm
    f = w_ref/60; % frequency, Hz
    W = 2*pi*f;
    Np = 12; % pulses per revolution
    time_per_revolution = 1/f; % sec
    time_per_pulse = time_per_revolution/Np;
    time_recorded = 30; % sec
    n = 1 : ceil(time_recorded*Fs); % Number of samples recorded;
    t = n/Fs; % time vector


    t_step = floor(t/time_per_pulse)*time_per_pulse;

    x = cos(W*t_step);
    y = sin(W*t_step);

    % Simulation of fixed point
    % x = ((x+1)/2)*(2-2^(-16))-1;
    % x = floor(x*2^16)/2^16;
    % y = floor(y*2^15)/2^15;

    figure;
    subplot(2,1,1);
    plot(t, x);
    hold on;
    plot(t, y, 'color', 'red');
    grid on;


    Tsc = 16e-3; % sec (16 ms)
    wb = 1; % rad/sec
    q = 2^(-15);
    r = 1;

    H = [1,0,0;...
         0,1,0];

    Q = q*eye(3);

    R = r*eye(2);

    P_hat = zeros(3);
    x_tilda = zeros(3,1);
    x_hat   = [0,0,0]';
    y_tilda = zeros(2,1);
    x_out = zeros(1, numel(x));
    y_out = zeros(1, numel(x));
    w_out = zeros(1, numel(x));


    for i=1:numel(x)

        x1 = x_hat(1);
        x2 = x_hat(2);
        x3 = x_hat(3);

        f_init = [cos(wb*Tsc*x3), -sin(wb*Tsc*x3), 0;...
                  sin(wb*Tsc*x3),  cos(wb*Tsc*x3), 0;...
                  0             ,  0             , 1];

    % Step I: df/dx
        F = [cos(wb*Tsc*x3), -sin(wb*Tsc*x3), -x1*wb*Tsc*sin(wb*Tsc*x3)-x2*wb*Tsc*cos(wb*Tsc*x3);...
             sin(wb*Tsc*x3),  cos(wb*Tsc*x3),  x1*wb*Tsc*cos(wb*Tsc*x3)-x2*wb*Tsc*sin(wb*Tsc*x3);...
             0             ,  0             , 1];

        x_tilda = f_init*x_hat;
        y_tilda = H*x_tilda;

    % Step II
        P_tilda = F*P_hat*F' + Q;

    % Step III
        K = P_tilda*H'/((H*P_tilda*H' + R));

        y_observed = [x(i); y(i)];

    % Step IV
        x_hat = x_tilda + K*(y_observed - y_tilda);

    % Step V
        P_hat = P_tilda - K*H*P_tilda;

        
        x_out(i) = y_tilda(1);
        y_out(i) = y_tilda(2);
        w_out(i) = wb*Tsc*x3*Fs; % angle by which the encoder is moving per sample divided by time taken

    end

    plot(t, x_out);
    plot(t, y_out, 'color', 'red');
    title ('x and y postion of encoder');
    xlabel('Time (sec)')

    subplot(2,1,2);
    plot(t, w_out*60/(2*pi));
    ylabel ('Encoder frequency (RPM)');
    xlabel('Time (sec)')
    grid minor;


    disp('V-test done!');
end