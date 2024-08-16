% Define different values of M (number of symbols)
M_values = [2, 4, 8, 16, 32];

% Eb/N0 range in dB
Eb_N0_dB = 0:20;

% Convert Eb/N0 from dB to linear scale
Eb_N0 = 10.^(Eb_N0_dB / 10);

% Initialize arrays to store results
Ps_mpsk = zeros(length(M_values), length(Eb_N0));
Ps_mpam = zeros(length(M_values), length(Eb_N0));
Ps_mqam = zeros(length(M_values), length(Eb_N0));
Ps_mfsk = zeros(length(M_values), length(Eb_N0));
Ps_mfsk_noncoherent = zeros(length(M_values), length(Eb_N0));

% Loop over each value of M for each modulation scheme
for i = 1:length(M_values)
    M = M_values(i);
    
    % MPSK
    Pe_mpsk = qfunc(sqrt(2 * Eb_N0 * sin(pi/M)^2));
    Ps_mpsk(i, :) = 1 - (1 - Pe_mpsk).^log2(M);
    
    % MPAM
    Pe_mpam = qfunc(sqrt(3 * log2(M) * Eb_N0 / (M^2 - 1)));
    Ps_mpam(i, :) = 1 - (1 - Pe_mpam).^log2(M);
    
    % MQAM
    Pe_mqam = 2 * (1 - 1/sqrt(M)) * qfunc(sqrt(3 * log2(M) * Eb_N0 / (2 * (M - 1))));
    Ps_mqam(i, :) = 1 - (1 - Pe_mqam).^log2(M);
    
    % MFSK
    Ps_mfsk(i, :) = 2 * qfunc(sqrt((6 * log2(M)) / (M^2 - 1) * Eb_N0));
    
    % Non-coherent MFSK
    Ps_mfsk_noncoherent(i, :) = 2 * qfunc(sqrt((3 * log2(M)) / (2 * (M^2 - 1)) * Eb_N0));
end

% Plotting
figHandles = gobjects(5, 1);

% MPSK
figHandles(1) = figure;
semilogy(Eb_N0_dB, Ps_mpsk','-o');
title('Probability of Symbol Error for MPSK');
legend('M=2', 'M=4', 'M=8', 'M=16', 'M=32');
grid on;

% MPAM
figHandles(2) = figure;
semilogy(Eb_N0_dB, Ps_mpam');
title('Probability of Symbol Error for MPAM');
legend('M=2', 'M=4', 'M=8', 'M=16', 'M=32');
grid on;

% MQAM
figHandles(3) = figure;
semilogy(Eb_N0_dB, Ps_mqam','-o');
title('Probability of Symbol Error for MQAM');
legend('M=2', 'M=4', 'M=8', 'M=16', 'M=32');
grid on;

% MFSK
figHandles(4) = figure;
semilogy(Eb_N0_dB, Ps_mfsk');
title('Probability of Symbol Error for MFSK');
legend('M=2', 'M=4', 'M=8', 'M=16', 'M=32');
grid on;

% Non-coherent MFSK
figHandles(5) = figure;
semilogy(Eb_N0_dB, Ps_mfsk_noncoherent');
title('Probability of Symbol Error for Non-coherent MFSK');
legend('M=2', 'M=4', 'M=8', 'M=16', 'M=32');
grid on;

% Set common xlabel and ylabel
for i = 1:length(figHandles)
    figure(figHandles(i));
    xlabel('Eb/N0 (dB)');
    ylabel('Probability of Symbol Error (Ps)');
end