clc;
clear;

% Parameters
N = 1024;                  % Number of subcarriers
cp_len = 128;              % Length of cyclic prefix
SNR_dB = 0:2:30;           % Extended SNR range in dB (from 0 to 30 dB)

% Modulation schemes for MPSK and MQAM
psk_mod_orders = [4, 8, 16, 32, 64];   % MPSK modulation orders
qam_mod_orders = [4, 16, 64, 256];     % MQAM modulation orders

% Preallocate error rate arrays
ser_psk = zeros(length(psk_mod_orders), length(SNR_dB));
ser_qam = zeros(length(qam_mod_orders), length(SNR_dB));

% Simulate for MPSK
for m = 1:length(psk_mod_orders)
    M = psk_mod_orders(m);
    
    for snr_idx = 1:length(SNR_dB)
        % Number of errors and symbols
        num_errors = 0;
        num_symbols = 0;
        
        % Run simulation for each SNR value
        for blk = 1:100
            % Transmitter
            data = randi([0 M-1], N, 1);                  % Random data symbols
            mod_data = pskmod(data, M, 0, 'gray');        % Modulate data using PSK
            
            % OFDM Modulation
            ifft_data = ifft(mod_data, N);                % IFFT operation
            cp_data = [ifft_data(end-cp_len+1:end); ifft_data]; % Add cyclic prefix
            
            % Channel (AWGN)
            noise = (1/sqrt(2)) * (randn(size(cp_data)) + 1i*randn(size(cp_data)));
            rx_signal = cp_data + 10^(-SNR_dB(snr_idx)/20) * noise;
            
            % Receiver
            rx_signal = rx_signal(cp_len+1:end);          % Remove cyclic prefix
            fft_data = fft(rx_signal, N);                 % FFT operation
            
            % Demodulation
            demod_data = pskdemod(fft_data, M, 0, 'gray');% Demodulate data using PSK
            
            % Calculate SER
            num_errors = num_errors + sum(data ~= demod_data);
            num_symbols = num_symbols + length(data);
        end
        
        % Symbol Error Rate (SER)
        ser_psk(m, snr_idx) = num_errors / num_symbols;
    end
end

% Simulate for MQAM
for m = 1:length(qam_mod_orders)
    M = qam_mod_orders(m);
    
    for snr_idx = 1:length(SNR_dB)
        % Number of errors and symbols
        num_errors = 0;
        num_symbols = 0;
        
        % Run simulation for each SNR value
        for blk = 1:100
            % Transmitter
            data = randi([0 M-1], N, 1);                 % Random data symbols
            mod_data = qammod(data, M, 'gray');          % Modulate data using QAM
            
            % OFDM Modulation
            ifft_data = ifft(mod_data, N);               % IFFT operation
            cp_data = [ifft_data(end-cp_len+1:end); ifft_data]; % Add cyclic prefix
            
            % Channel (AWGN)
            noise = (1/sqrt(2)) * (randn(size(cp_data)) + 1i*randn(size(cp_data)));
            rx_signal = cp_data + 10^(-SNR_dB(snr_idx)/20) * noise;
            
            % Receiver
            rx_signal = rx_signal(cp_len+1:end);         % Remove cyclic prefix
            fft_data = fft(rx_signal, N);                % FFT operation
            
            % Demodulation
            demod_data = qamdemod(fft_data, M, 'gray');  % Demodulate data using QAM
            
            % Calculate SER
            num_errors = num_errors + sum(data ~= demod_data);
            num_symbols = num_symbols + length(data);
        end
        
        % Symbol Error Rate (SER)
        ser_qam(m, snr_idx) = num_errors / num_symbols;
    end
end

% Plotting the results for MPSK
figure;
for m = 1:length(psk_mod_orders)
    semilogy(SNR_dB, ser_psk(m, :), '-o', 'LineWidth', 2);
    hold on;
end
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Symbol Error Rate (P_s)');
legend(arrayfun(@(x) [num2str(x) '-PSK'], psk_mod_orders, 'UniformOutput', false), 'Location', 'southwest');
title('Symbol Error Rate Performance of MPSK-CP-OFDM over AWGN Channel');

% Plotting the results for MQAM
figure;
for m = 1:length(qam_mod_orders)
    semilogy(SNR_dB, ser_qam(m, :), '-o', 'LineWidth', 2);
    hold on;
end
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Symbol Error Rate (P_s)');
legend(arrayfun(@(x) [num2str(x) '-QAM'], qam_mod_orders, 'UniformOutput', false), 'Location', 'southwest');
title('Symbol Error Rate Performance of MQAM-CP-OFDM over AWGN Channel');