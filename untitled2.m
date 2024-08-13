close all
clear all
clc

% Inicialização dos parâmetros
N_max = 16; % Número máximo de elementos da RIS
Es = 1; % Energia média do sinal
N0 = 1; % Densidade espectral do ruído
MC = 1e4; % Iterações de Monte Carlo
N_s = 1e3; % Quantidade de simbolos transmitidos
SNR_teo_dB = -10:2:20;
SNR_teo_lin = 10.^(SNR_teo_dB/10); 

% Inicialização das variáveis
BER_sim = zeros(1, length(SNR_teo_dB));
BER_sim_q = zeros(1, length(SNR_teo_dB));
eta_ideal = zeros(1, length(SNR_teo_dB));
eta_quantizado = zeros(1, length(SNR_teo_dB));

% Loop de Monte Carlo
for k = 1:length(SNR_teo_dB)
    SNR = 10^(SNR_teo_dB(k) / 10); % SNR linear
    
    totalErrors = 0;
    totalErrors_q = 0;

    for i = 1:MC
        % Gerar bits aleatórios
        dataBits = randi([0 1], N_s, 1);
        
        % Modulação BPSK
        mSignal = 2 * dataBits - 1; % Mapear 0 -> -1, 1 -> 1
        
        % Adicionar ruído AWGN
        % Calcular a SNR com RIS desligada
        h = (1/sqrt(2) * (randn(1, N_max) + 1i * randn(1, N_max)));
        g = (1/sqrt(2) * (randn(1, N_max) + 1i * randn(1, N_max)));
        ephi_off = exp(1i * 0);
        
        % SNR com RIS desligada
        off = (abs(sum(h .* ephi_off .* g)))^2;
        gamma_off = (off * Es) / N0; % Linear
        
        % SNR com RIS ligada
        phi = -angle(h) - angle(g);
        ephi = exp(1i * phi); % Ajuste ótimo
        on = (abs(sum(h .* ephi .* g)))^2;
        gamma_on = (on * Es) / N0; % Linear
        
        % Adicionar ruído
        y = awgn((gamma_off * mSignal), SNR_teo_dB(k), 'measured');
        y_q = awgn((gamma_on * mSignal), SNR_teo_dB(k), 'measured');

        % Demodulação BPSK manual
        dBits = y > 0;
        dBits_q = y_q > 0;
        
        % Calcular erros de bit
        errors = sum(dataBits ~= dBits);
        errors_q = sum(dataBits ~= dBits_q);

        totalErrors = totalErrors + errors;
        totalErrors_q = totalErrors_q + errors_q;
    end

    % Calcular BER
    BER_sim(k) = totalErrors / (MC * N_s);
    BER_sim_q(k) = totalErrors_q / (MC * N_s);
    
    % Calcular eficiência espectral
    eta_ideal(k) = log2(1 + SNR); % Eficiência ideal
    eta_quantizado(k) = (1 - BER_sim_q(k)) * log2(1 + SNR); % Eficiência quantizada
end

% Plotar
figure;
hold on;

% Plotar caso ideal
plot(SNR_teo_dB, eta_ideal, '-o', 'LineWidth', 2, 'DisplayName', 'Caso Ideal');

% Plotar caso quantizado
plot(SNR_teo_dB, eta_quantizado, '-x', 'LineWidth', 2, 'DisplayName', 'Caso Quantizado');

hold off;
xlabel('Eb/N0 (dB)', 'Interpreter', 'latex');
ylabel('Eficiência Espectral (bps/Hz)');
title('Eficiência Espectral em Função do $E_b/N_0$', 'Interpreter', 'latex');
legend show;
grid on;
