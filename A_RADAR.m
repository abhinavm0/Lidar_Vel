function A_RADAR

    % A_RADAR.m
    clc;
    close all;
    clear;

%% Design Specifications
    
    pd = 0.9;            % Probability of detection
    pfa = 1e-6;          % Probability of false alarm
    max_range = 5000;    % Maximum unambiguous range
    range_res = 50;      % Required range resolution
    tgt_rcs = 1;         % Required target radar cross section

%% Waveform
    
    prop_speed = 2.99792458e8; % physconst('LightSpeed');   % Propagation speed
    pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
    pulse_width = 1/pulse_bw;               % Pulse width
    prf = prop_speed/(2*max_range);         % Pulse repetition frequency
    fs = 2*pulse_bw;                        % Sampling rate
%     waveform = phased.RectangularWaveform(...
%         'PulseWidth',1/pulse_bw,...
%         'PRF',prf,...
%         'SampleRate',fs);
    waveform = rect_pulse_train(pulse_width, prf, fs); % Using basic Matlab functions, recreated the functionality as described in ref to phased.RectangularWaveform

%% Receiver Noise Characteristics
    
    noise_bw = pulse_bw;
    gain = 20;
    noisefigure = 0;
    samplerate = fs;
    en_port = 1;

%     receiver = phased.ReceiverPreamp(...
%         'Gain',20,...
%         'NoiseFigure',0,...
%         'SampleRate',fs,...
%         'EnableInputPort',true);

    receiver = phased_ReceiverPreamp(gain, noisefigure, samplerate, en_port); % Using basic Matlab functions, recreated the functionality as described in ref to phased.ReceiverPreamp

%% Transmitter
    % snr_min = albersheim(pd, pfa, num_pulse_int)
    snr_min = 4.9904;
    tx_gain_db = 20;

    fc = 10e9; % operating frequency 
    lambda = prop_speed/fc;

%     peak_power = radareqpow(lambda,max_range,snr_min,pulse_width,...
%         'RCS',tgt_rcs,'Gain',tx_gain);
        
    peak_power = get_radareqpow(lambda,max_range,snr_min,pulse_width,...
        tgt_rcs,tx_gain_db); % Using basic Matlab functions, recreated the functionality as described in ref to radareqpow

    
%     transmitter = phased.Transmitter(...
%         'Gain',tx_gain,...
%         'PeakPower',peak_power,...
%         'InUseOutputPort',true);

    transmitter = phased_Transmitter(tx_gain_db,peak_power,1); % Using basic Matlab functions, recreated the functionality as described in ref to phased.Transmitter

    % Antenna design to follow.
    
    disp('Radar done!');

end

%% Helper functions

function w = rect_pulse_train(pulse_width, prf, fs)
    time = 0.0001; % sec
    N_samples = time*fs;
    n = 1 : N_samples;
    w = zeros(1, N_samples);
%     f = find(mod(n/fs, 1/prf) <= pulse_width);
    w(mod(n/fs, 1/prf) <= pulse_width) = 1;
%     plot(w);    
    disp('rect_pulse_train done!');
end

function receiver = phased_ReceiverPreamp(x, gain, noisefigure, fs, en_port)
    Boltzmann_const = 1.38e-23;
    noise_temp = 290; % Kelvin (assumption)
    NoiseBandwidth = fs;
    G = 10^(gain/20);
    noisepow = Boltzmann_const*noise_temp*NoiseBandwidth;
    receiver = G*x + sqrt(noisepow/2)*(randn(size(x))+1j*randn(size(x)));
    disp('phased_ReceiverPreamp done!');
end

function peak_power = get_radareqpow(lambda,max_range,snr_min_db,pulse_width,...
        tgt_rcs,tx_gain_db)
    
    noise_temp = 290; % Kelvin (assumption)
    Boltzmann_const = 1.38064852e-23;
    tx_gain = db2pow(tx_gain_db);
    
    snr_min = db2pow(snr_min_db);
    
    peak_power = (snr_min)*((4*pi)^3)*Boltzmann_const*noise_temp*(max_range^4)/...
        (pulse_width*(tx_gain^2)*(lambda^2)*tgt_rcs);
    disp(['peak_power = ', num2str(peak_power, '%10.4e')]); 
    disp('get_radareqpow done!');
end

function transmitter = phased_Transmitter(tx_gain_db, PeakPower, InUseOutputPort)
    tx_gain = 10^(tx_gain_db/10);    
    LossFactor = 0;
    transmitter = sqrt(PeakPower*db2pow(tx_gain - LossFactor));
    disp(['transmitter = ', num2str(transmitter, '%10.4e')]); 
    disp('phased_Transmitter done!');

end