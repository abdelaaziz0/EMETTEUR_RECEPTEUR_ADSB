%% Définition des paramètres et des échantillons
N0 = 112;                % Nombre de bits
Ts = 10^(-6);         % Durée symbole
fe = 20 * 10^6;       % Fréquence d'échantillonnage réelle
Te = 1/fe;            % Période d'échantillonnage
Fse = Ts/Te;          % Facteur de sur-échantillonnage
b_k = randi([0,1], 1, N0);    % Bits générés
sigma_nl = 0.5;       % Écart-type du bruit

%% Génération du préambule
Tp = 8e-6; % Durée du préambule (8 µs)
unite = 0.5e-6; % Durée d'une unité (0.5 µs)
s_p = zeros(1, floor(Tp/Te));
for k = 1:length(s_p)
    if (k >= 1 && k < unite/Te) || (k >= 2*unite/Te && k < 3*unite/Te) || ...
       (k >= 7*unite/Te && k < 8*unite/Te) || (k >= 9*unite/Te && k < 10*unite/Te)
        s_p(k) = 1;
    end
end

%% La modulation PPM
s_l = zeros(N0 * Fse, 1);
for k = 1:length(s_l)
    s_l(k) = PPM(k* Te, b_k, Ts);  % Signal PPM
end
s_l(s_l>1)=1;

%% Préambule
s_l_with_preamble = [s_p' ; s_l];

%%
% Modélisation des distorsions
delta_t = Te * rand() * 100; % Délai de propagation aléatoire entre 0 et 100Te
delta_f = (rand() * 2 - 1) * 1000; % Décalage en fréquence aléatoire entre -1kHz et 1kHz
phi_0 = 2*pi*rand(); % Déphasage aléatoire entre 0 et 2π
alpha = 1; % Pas d'atténuation

% Application des distorsions
s_l_distorted = [zeros(round(delta_t/Te), 1); s_l_with_preamble; zeros(100, 1)];
t = (0:length(s_l_distorted)-1)*Te;
s_l_distorted = alpha * s_l_distorted .* exp(-1j * 2*pi * delta_f * t' + 1j * phi_0);

% Ajout du bruit n_l(t)
n_l = sigma_nl * (randn(size(s_l_distorted)) + 1j*randn(size(s_l_distorted))) / sqrt(2);  
y_l = s_l_distorted + n_l;

% Synchronisation
y_l1=y_l;
[y_l, delta_t_est] = synchronize_signal(y_l, s_p, Te);

preamble_length = length(s_p);
y_l = y_l(preamble_length+1:end);

%% Filtrage adapté
g_a = ones(Ts/(2*Te), 1);           % Réponse impulsionnelle du filtre
r_l_1 = conv(y_l, g_a);     % Filtrage avec convolution

%% Sur-échantillonnage avec Ts/2
r_l = zeros(length(b_k) * 2, 1);
for k = 1:length(r_l)
    r_l(k) = r_l_1(10 * (k-1) + 10);
end

%% Comparaison des bits originaux avec les bits estimés
v0 = Ts/2;  % Valeur de référence
b_k_recus = zeros(1, length(b_k));

for k = 1:length(b_k)
    rk = [r_l(2*k-1); r_l(2*k)];
    v1 = v0 * [0; 1];  % Vecteur pour bit 0
    v2 = v0 * [1; 0];  % Vecteur pour bit 1
    
    d1 = norm(rk - v1)^2;  % Distance au carré pour bit 0
    d2 = norm(rk - v2)^2;  % Distance au carré pour bit 1
    
    if d1 < d2
        b_k_recus(k) = 0;
    else
        b_k_recus(k) = 1;
    end
end

%% Comparaison des bits originaux avec les bits reçus
errors = sum(b_k ~= b_k_recus);
disp(['Nombre d''erreurs: ', num2str(errors)]);
% Affichage des résultats
disp(['Délai réel: ', num2str(delta_t*1e6), ' µs']);
disp(['Délai estimé: ', num2str(delta_t_est*1e6), ' µs']);

%% Visualisations originales
figure;

% Bits générés
subplot(5,1,1);
stem(b_k, 'filled');
title('Bits générés');
xlabel('Index du bit');
ylabel('Valeur du bit');
grid on;

% Signal modulé PPM
subplot(5,1,2);
t_ppm = (0:length(s_l)-1) * Te;
plot(t_ppm, s_l);
title('Signal modulé PPM');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Signal avec bruit (après synchronisation)
subplot(5,1,3);
t_sync = (0:length(y_l)-1) * Te;
plot(t_sync, real(y_l));
title('Signal synchronisé avec bruit');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Signal après filtrage adapté
subplot(5,1,4);
t_filtered = (0:length(r_l_1)-1) * Te;
plot(t_filtered, real(r_l_1));
title('Signal après filtrage adapté');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Signal suréchantillonné
subplot(5,1,5);
stem(real(r_l));
title('Signal suréchantillonné');
xlabel('Index');
ylabel('Amplitude');
grid on;

%% Nouvelles visualisations
figure;
subplot(5,1,1);
plot(real(s_l_with_preamble));
title('Signal original avec préambule');
xlabel('Échantillons'); ylabel('Amplitude');

subplot(5,1,2);
plot(real(y_l1));
title('Signal reçu (avec bruit et distorsions)');
xlabel('Échantillons'); ylabel('Amplitude');

subplot(5,1,3);
plot(real(y_l));  % Nous utilisons y_l ici car c'est le signal après synchronisation
title('Signal synchronisé en temps');
xlabel('Échantillons'); ylabel('Amplitude');

subplot(5,1,4);
plot(real(y_l));  % Nous n'avons pas de y_l_freq_corrected, donc nous utilisons y_l
title('Signal synchronisé en temps et en fréquence');
xlabel('Échantillons'); ylabel('Amplitude');

subplot(5,1,5);
stem(1:N0, b_k);
hold on;
stem(1:N0, b_k_recus, 'r');
title('Bits originaux (bleu) vs Bits décodés (rouge)');
xlabel('Index du bit'); ylabel('Valeur du bit');

% Affichage des bits reçus après décision
figure;
stem(b_k_recus, 'filled');
title('Bits reçus après décision');
xlabel('Index du bit');
ylabel('Valeur du bit');
grid on;


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
    
    % Correction du délai
    start_index = round(delta_t_est / Te);
    
    % Éviter les index hors limites
    if start_index < 1
        start_index = 1;
    elseif start_index > length(y_l)
        start_index = length(y_l);
    end
    
    y_l_sync = y_l(start_index:end);
    
    % Normalisation de l'amplitude
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