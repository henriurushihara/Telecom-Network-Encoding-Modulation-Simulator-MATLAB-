
clc;
close all;

% NOTE : YOU MUST REVIEW + ADJUST PARAMETERS BELOW
% for ALL ENCODING SCHEMES!!!

% SAMPLE BIT PATTERN 
bit_pattern = [0 1 0 0 1 1 1 0];

% OTHER SAMPLE BIT PATTERNS YOU CAN TRY. 

% Define the bit pattern
%bit_pattern = [0 0 1 1 1 1 0 1 0 0 1 1 0 0 1 1];

% bit_pattern = [0 1 0 0 1 1 0 0 0 1 1];

% bit_pattern = [0 0 1 1 0 1 0 0 0 1 0];

% bit_pattern = [0 0 1 1 1 0 0 1 0 1];

% bit_pattern = [1 0 1 1 0];

% Length of bit pattern
N = length(bit_pattern);

% NRZ-Level Encoding
nrz_signal = zeros(1, N);
for i = 1:N
    if bit_pattern(i) == 0
        nrz_signal(i) = 1;
    else
        nrz_signal(i) = 0;
    end
end

% NRZI Encoding
nrzi_signal = zeros(1, N);

%prev_state = 0;  % ASSUMPTION : Prev. voltage is low; 

prev_state = 1;  % ASSUMPTION : Prev. voltage is HIGH; 

for i = 1:N
    if bit_pattern(i) == 1
        nrzi_signal(i) = ~prev_state;
        prev_state = nrzi_signal(i);
    else
        nrzi_signal(i) = prev_state;
    end
end

% Bipolar AMI Encoding
ami_signal = zeros(1, N);

prev_ami = -1;  % Most recent preceding 1 bit has a negative voltage

% prev_ami = 1;  % Most recent preceding 1 bit has a POSITIVE voltage

for i = 1:N
    if bit_pattern(i) == 1
        prev_ami = -prev_ami;  % Alternate the polarity for next 1
        ami_signal(i) = prev_ami;
    end
end

% Pseudoternary Encoding
pseudo_signal = zeros(1, N);

prev_pseudo = -1;  % Most recent preceding 0 bit has a negative voltage

% prev_pseudo = 1;  % Most recent preceding 0 bit has a POSITIVE voltage

for i = 1:N
    if bit_pattern(i) == 0
        prev_pseudo = -prev_pseudo;  % Alternate the polarity for next 0
        pseudo_signal(i) = prev_pseudo;
    end
end

% Manchester Encoding
manchester_signal = zeros(1, 2*N);
for i = 1:N
    if bit_pattern(i) == 0
        manchester_signal(2*i-1) = 1;
        manchester_signal(2*i) = 0;
    else
        manchester_signal(2*i-1) = 0;
        manchester_signal(2*i) = 1;
    end
end

% Differential Manchester Encoding
diff_manchester_signal = zeros(1, 2*N);

prev_state = 0;  % Start with an assumed high state 

% ( THIS DEPENDS ON THE FIRST BIT OF THE SEQ!!!; 
% if 1st bit is 0 => prev_state = 1
% if 1st bit is 1 => prev_state = 0 )

for i = 1:N
    if bit_pattern(i) == 0
        diff_manchester_signal(2*i-1) = prev_state;
        diff_manchester_signal(2*i) = ~prev_state;
    else
        diff_manchester_signal(2*i-1) = ~prev_state;
        diff_manchester_signal(2*i) = prev_state;
        prev_state = ~prev_state;  % Change the state after a 1 bit
    end
end

% Plotting
figure;
subplot(6,1,1);
stairs([0.5:1:(N-0.5)], nrz_signal);
title('NRZ-L');
xlim([0.5,N+1])
ylim([-0.5 1.5]);
xlabel('Bit');
ylabel('Amplitude');
grid on;

subplot(6,1,2);
stairs([0.5:1:(N-0.5)], nrzi_signal);
title('NRZI');
xlim([0.5,N+1])
ylim([-0.5 1.5]);
xlabel('Bit');
ylabel('Amplitude');
grid on;

subplot(6,1,3);
stairs([0.5:1:(N-0.5)], ami_signal);
title('Bipolar AMI');
xlim([0.5,N+1])
ylim([-1.5 1.5]);
xlabel('Bit');
ylabel('Amplitude');
grid on;

subplot(6,1,4);
stairs([0.5:1:(N-0.5)], pseudo_signal);
title('Pseudoternary');
xlim([0.5,N+1])
ylim([-1.5 1.5]);
xlabel('Bit');
ylabel('Amplitude');
grid on;

subplot(6,1,5);
x = 0.5:0.5:N;

y = manchester_signal(1:16);

stairs(x, y);
title('Manchester');
xlim([0.5,N+1])
ylim([-0.5 1.5]);
xlabel('Bit');
ylabel('Amplitude');
grid on;

subplot(6,1,6);
x = 0.5:0.5:N;
y = diff_manchester_signal(1:2*N);
stairs(x, y);
title('Differential Manchester');
xlim([0.5,N+1])
ylim([-0.5 1.5]);
xlabel('Bit');
ylabel('Amplitude');
grid on;

% PSK AND DPSK

% Parameters for PSK and DPSK
f = 1;  % periods/cycles per bit
fs = 100;  % Sampling frequency of 100 samples/second
t = linspace(0, N, N*fs);  % Time vector
z = pi;  % Phase shift for 1 (can be adjusted)

% PSK Modulation
psk_signal = zeros(1, N*fs);
for i = 1:N
    if bit_pattern(i) == 0
        psk_signal((i-1)*fs+1 : i*fs) = sin(2*pi*f*t((i-1)*fs+1 : i*fs));
    else
        psk_signal((i-1)*fs+1 : i*fs) = sin(2*pi*f*t((i-1)*fs+1 : i*fs) + z);
    end
end

% DPSK Modulation
dpsk_signal = zeros(1, N*fs);
prev_phase = 0;
for i = 1:N
    if bit_pattern(i) == 0
        dpsk_signal((i-1)*fs+1 : i*fs) = sin(2*pi*f*t((i-1)*fs+1 : i*fs) + prev_phase);
    else
        prev_phase = mod(prev_phase + z, 2*pi);
        dpsk_signal((i-1)*fs+1 : i*fs) = sin(2*pi*f*t((i-1)*fs+1 : i*fs) + prev_phase);
    end
end

% Adding plots for PSK and DPSK
figure;  % Create a new figure

subplot(2,1,1);
plot(t, psk_signal);
title('PSK');
xlim([0 N])
ylim([-1.5 1.5]);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, dpsk_signal);
title('Differential PSK');
xlim([0 N])
ylim([-1.5 1.5]);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
