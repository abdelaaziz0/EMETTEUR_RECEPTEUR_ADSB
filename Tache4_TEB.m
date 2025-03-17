% Paramètres
N0 = 112;                % Nombre de bits par trame
Ts = 10^(-6);            % Durée symbole
fe = 20 * 10^6;          % Fréquence d'échantillonnage
Te = 1/fe;               % Période d'échantillonnage
Fse = Ts/Te;             % Facteur de sur-échantillonnage
Nb_trames = 1000;        % Nombre de trames à simuler par point Eb/N0
Eb_N0_dB = 0:1:10;       % Plage de Eb/N0 en dB

% Génération du préambule
Tp = 8e-6;               % Durée du préambule (8 µs)
unite = 0.5e-6;          % Durée d'une unité (0.5 µs)
s_p = zeros(1, floor(Tp/Te));
for k = 1:length(s_p)
    if (k >= 1 && k < unite/Te) || (k >= 2*unite/Te && k < 3*unite/Te) || ...
       (k >= 7*unite/Te && k < 8*unite/Te) || (k >= 9*unite/Te && k < 10*unite/Te)
        s_p(k) = 1;
    end
end

% Initialisation du vecteur TEB
TEB = zeros(size(Eb_N0_dB));

disp('Début de la simulation');

% Boucle sur les valeurs de Eb/N0
for i_EbN0 = 1:length(Eb_N0_dB)
    fprintf('Simulation pour Eb/N0 = %.2f dB\n', Eb_N0_dB(i_EbN0));
    
    Eb_N0 = 10^(Eb_N0_dB(i_EbN0)/10);
    sigma_nl = sqrt(1/(2*Eb_N0));  % Écart-type du bruit
    
    nb_erreurs = 0;
    nb_bits_total = 0;
    
    % Boucle sur les trames
    for i_trame = 1:Nb_trames
        if mod(i_trame, 100) == 0
            fprintf('  Traitement de la trame %d/%d\n', i_trame, Nb_trames);
        end
        
        % Génération des bits
        b_k = randi([0,1], 1, N0);
        
        % Modulation PPM
        s_l = zeros(N0 * Fse, 1);
        for k = 1:length(s_l)
            s_l(k) = PPM(k* Te, b_k, Ts);
        end
        s_l(s_l>1)=1;
        
        % Ajout du préambule
        s_l_with_preamble = [s_p' ; s_l];
        
        % Modélisation des distorsions
        delta_t = Te * rand() * 100;  % Délai de propagation aléatoire entre 0 et 100Te
        delta_f = 0*(rand() * 2 - 1) * 1000;  % Décalage en fréquence aléatoire entre -1kHz et 1kHz
        phi_0 = 0*2*pi*rand();  % Déphasage aléatoire entre 0 et 2π
        
        % Application des distorsions
        s_l_distorted = [zeros(round(delta_t/Te), 1); s_l_with_preamble; zeros(100, 1)];
        t = (0:length(s_l_distorted)-1)*Te;
        s_l_distorted = s_l_distorted .* exp(-1j * 2*pi * delta_f * t' + 1j * phi_0);
        
        % Ajout du bruit
        n_l = sigma_nl * (randn(size(s_l_distorted)) + 1j*randn(size(s_l_distorted))) / sqrt(2);  
        y_l = s_l_distorted + n_l;
        
        % Synchronisation
        [y_l, delta_t_est] = synchronize_signal(y_l, s_p, Te);
        
        % Suppression du préambule
        preamble_length = length(s_p);
        y_l = y_l(preamble_length+1:end);
        
        % Correction de la fréquence (simplifiée)
        t = (0:length(y_l)-1)*Te;
        y_l = y_l .* exp(1j * 2*pi * delta_f * t');
        
        % Filtrage adapté
        g_a = ones(Ts/(2*Te), 1);
        r_l_1 = conv(y_l, g_a);
        
        % Sur-échantillonnage avec Ts/2
        r_l = zeros(length(b_k) * 2, 1);
        for k = 1:length(r_l)
            r_l(k) = r_l_1(10 * (k-1) + 10);
        end
        
        % Décision
        v0 = Ts/2;
        b_k_recus = zeros(1, length(b_k));
        for k = 1:length(b_k)
            rk = [r_l(2*k-1); r_l(2*k)];
            v1 = v0 * [0; 1];
            v2 = v0 * [1; 0];
            d1 = norm(rk - v1)^2;
            d2 = norm(rk - v2)^2;
            if d1 < d2
                b_k_recus(k) = 0;
            else
                b_k_recus(k) = 1;
            end
        end
        
        % Comptage des erreurs
        nb_erreurs = nb_erreurs + sum(b_k ~= b_k_recus);
        nb_bits_total = nb_bits_total + length(b_k);
    end
    
    % Calcul du TEB
    TEB(i_EbN0) = nb_erreurs / nb_bits_total;
    fprintf('TEB pour Eb/N0 = %.2f dB : %.4e\n', Eb_N0_dB(i_EbN0), TEB(i_EbN0));
end

disp('Simulation terminée');

% Calcul du TEB théorique
TEB_theo = 0.5 * erfc(sqrt(10.^(Eb_N0_dB/10)));

% Affichage des résultats
figure;
semilogy(Eb_N0_dB, TEB, 'o-', Eb_N0_dB, TEB_theo, 'r-');
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Taux d''Erreur Binaire');
legend('TEB simulé', 'TEB théorique');
title('TEB en fonction de Eb/N0');

% Calcul de la perte en dB pour TEB = 10^-3
TEB_target = 1e-3;
Eb_N0_theo = interp1(log10(TEB_theo), Eb_N0_dB, log10(TEB_target));

% Trouver le premier TEB non nul
non_zero_TEB = TEB(TEB > 0);
non_zero_Eb_N0 = Eb_N0_dB(TEB > 0);

if ~isempty(non_zero_TEB)
    Eb_N0_simu = interp1(log10(non_zero_TEB), non_zero_Eb_N0, log10(TEB_target));
    perte_dB = Eb_N0_simu - Eb_N0_theo;
    fprintf('Perte en dB pour TEB = 10^-3 : %.2f dB\n', perte_dB);
else
    fprintf('Impossible de calculer la perte en dB : tous les TEB sont nuls.\n');
end
% Fonction de synchronisation
function [y_l_sync, delta_t_est] = synchronize_signal(y_l, s_p, Te)
    % Longueur du préambule
    preamble_length = length(s_p);
    
    % Initialisation du vecteur de corrélation
    correlation = zeros(1, length(y_l) - preamble_length + 1);
    
    % Transposer s_p pour qu'il soit un vecteur colonne
    s_p_col = s_p(:);
    
    % Calcul de la corrélation normalisée
    for i = 1:(length(y_l) - preamble_length + 1)
        segment = y_l(i:i + preamble_length - 1);
        
        if length(segment) == preamble_length
            product = segment .* conj(s_p_col);
            correlation_value = abs(sum(product)) / ...
                (sqrt(sum(abs(s_p_col).^2)) * sqrt(sum(abs(segment).^2)));
            correlation(i) = correlation_value;
        else
            correlation(i) = 0;
        end
    end
    
    % Trouver l'indice du maximum de corrélation
    [~, max_index] = max(correlation);
    
    % Calculer le délai estimé
    delta_t_est = (max_index - 1) * Te;
    
    % Au lieu de couper le signal, on le décale circulairement
    y_l_sync = circshift(y_l, -max_index + 1);
    
    % Normalisation de l'amplitude si nécessaire
    y_l_sync = y_l_sync / max(abs(y_l_sync));
end

% Fonction PPM
function sl = PPM(t, A, Ts)
    sl = zeros(size(t));
    N = length(A);    
    
    for k = 1:N
        t_start = (k-1) * Ts;
        t_mid = t_start + Ts/2;
        t_end = k * Ts;
        
        mask_0 = (t > t_mid) & (t <= t_end);
        mask_1 = (t > t_start) & (t <= t_mid);
        
        if A(k) == 0
            sl = sl + mask_0;
        else
            sl = sl + mask_1;
        end
    end
end