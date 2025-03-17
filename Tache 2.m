% Paramètres du signal
nb_bits = 10000;  % Nombre de bits
Ts = 1e-6;        % Période symbole (1 µs)
Fe = 20e6;        % Fréquence d'échantillonnage (20 MHz)
NFFT = 512;      % Nombre de points pour la FFT

% Génération des bits aléatoires
bits = randi([0 1], 1, nb_bits);

% Génération du signal
t = 0:1/Fe:nb_bits*Ts-1/Fe;
s_l = generate_signal(bits, t, Ts);

% Calcul de la DSP théorique
f = linspace(-Fe/2, Fe/2, NFFT);
psd_theo = theoretical_psd(f, Ts);

% Estimation de la DSP avec Mon_Welch
psd_mon_welch= Mon_Welch(s_l, NFFT, Fe);
f_mon_welch = linspace(-Fe/2, Fe/2, NFFT);
% Estimation de la DSP avec pwelch
[psd_pwelch, f_pwelch] = pwelch(s_l, ones(1,NFFT), 0, NFFT, Fe, 'centered');

% Affichage des résultats
figure;
semilogy(f/1e6, psd_theo, 'r', 'LineWidth', 2);
hold on;
semilogy(f_mon_welch/1e6, psd_mon_welch, 'b--', 'LineWidth', 2);
semilogy(f_pwelch/1e6, psd_pwelch, 'g-.', 'LineWidth', 2);
xlabel('Fréquence (MHz)');
ylabel('DSP (W/Hz)');
legend('Théorique', 'Mon\_Welch', 'pwelch');
title('Comparaison des DSP théorique et estimées');
grid on;
xlim([-10 10]);
ylim([1e-8 1e-4]);  % Ajustez ces valeurs si nécessaire

function y= Mon_Welch(x, NFFT, Fe)
    L = length(x);
    K = floor(L/NFFT);
    segments = reshape(x(1:K*NFFT), NFFT, K);
    Y = fft(segments, NFFT);
    Pxx = mean(abs(Y).^2, 2) / (NFFT * Fe);
    y = fftshift(Pxx);
end

function s_l = generate_signal(bits, t, Ts)
    p_t = @(t) (abs(t) <= Ts/2) .* (1 - 2*abs(t)/Ts);
    s_l = 0.5 * ones(size(t));
    for k = 1:length(bits)
        A_k = 1 - 2*bits(k);
        s_l = s_l + A_k * p_t(t - (k-1)*Ts);
    end
end

function psd = theoretical_psd(f, Ts)
    % Calcul de la DSP théorique du signal ADS-B
    dc_component = 0.25 * (abs(f) < 1/2*Ts);
    envelope = 0.25 * Ts * sinc(f*Ts/2).^4;
    psd = dc_component;
    for n = -5:5
        psd = psd + envelope .* (abs(f - n/Ts) < 1/Ts);
    end
end