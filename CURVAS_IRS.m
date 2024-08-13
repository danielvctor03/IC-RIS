close all;
clear all;
clc;

% Inicialização dos parâmetros
N_max = 16; % Número máximo de elementos da RIS
Es = 1; % Energia média do sinal
N0 = 1; % Densidade espectral do ruído
MC = 1e4; % Iterações de Monte Carlo
N_s = 1e3; % Quantidade de simbolos transmitidos

SNR_teo_dB = -10:2:20;
SNR_teo_lin = 10.^(SNR_teo_dB/10);

gamma_on_avg = zeros(1, N_max); % Média da SNR com RIS ligada
gamma_off_avg = zeros(1, N_max); % Média da SNR com RIS desligada

BER_sim = zeros(1, length(SNR_teo_dB));

% Loop de Monte Carlo
for i = 1:N_max
    gamma_on_mc = zeros(1, MC);
    gamma_off_mc = zeros(1, MC);
    
    for j = 1:MC
        % Coeficientes dos canais
        h = (1/sqrt(2) * (randn(1, i) + 1i * randn(1, i)));
        theta = angle(h); % Fase de h
        g = (1/sqrt(2) * (randn(1, i) + 1i * randn(1, i)));
        psi = angle(g); % Fase de g 
        %%
        
        %%
        
        % SNR com RIS desligada
        phi_off = 0;
        ephi_off = exp(1i * (phi_off));
        off = (abs(sum(h .* ephi_off .* g)))^2;
        gamma_off = (off * Es) / N0; % Linear
        gamma_off_mc(j) = 10 * log10(gamma_off); % Em dB

        % SNR com RIS ligada
        phi = -(theta + psi);
        ephi = exp(1i * (phi)); % Ajuste ótimo
        on = (abs(sum(h .* ephi .* g)))^2;
        gamma_on = (on * Es) / N0; % Linear
        gamma_on_mc(j) = 10 * log10(gamma_on); % Em dB
    end

    % Média das SNRs para cada número de elementos
    gamma_on_avg(i) = mean(gamma_on_mc);
    gamma_off_avg(i) = mean(gamma_off_mc);
end

% Símbolos BPSK, ruído AWGN, demodulação e calculo da BER
for aux = 1:length(SNR_teo_dB)
    SNR = SNR_teo_dB(aux);
    totalErrors = 0;
    totalBits = 0;

    totalErrors_2 = 0;
    totalBits_2 = 0;
    
    for i = 1:MC
        % Gerando bits aleatórios
        dataBits = randi([0 1], N_s, 1);
        
        % Modulação BPSK
        modulatedSignal = 2*dataBits - 1; % Mapear 0 -> -1, 1 -> 1
        
        % Adicionando ruído AWGN
        %receivedSignal = awgn(modulatedSignal, SNR, 'measured');

        %y = h.*ephi.*g.*receivedSignal; % Sinal após influencia da RIS
        y = awgn((h*diag(ephi)*g.')*modulatedSignal, SNR, 'measured');

        % Demodulação BPSK manual
        demodulatedBits = y > 0;
        
        
        % Calculo dos erros de bit
        errors = sum(dataBits ~= demodulatedBits);
        

        totalErrors = totalErrors + errors;
        

        totalBits = totalBits + N_s;
        
        
    end

   

    % Calculo da BER
    BER_sim(aux) = totalErrors / totalBits;
    

end


% Plotagem da curva BER 
figure;
semilogy(SNR_teo_dB, BER_sim, 'o-', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Taxa de Erro de Bit (BER)');
title('BER vs SNR para Modulação BPSK');
grid on;

% Plotagem da curva SNR vs n° elementos
figure;
plot(1:N_max, gamma_on_avg, '-o', 'LineWidth', 2, 'DisplayName', 'RIS Ligada');
hold on;
plot(1:N_max, gamma_off_avg, '-x', 'LineWidth', 2, 'DisplayName', 'RIS Desligada');
hold off;
legend show;
xlabel('Número de Elementos da RIS');
ylabel('SNR (dB)');
title('SNR em função do Número de Elementos da RIS');
grid on;



