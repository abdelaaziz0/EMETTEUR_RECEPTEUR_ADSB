%% Paramètres CRC
crc_polynomial = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]; % Polynôme CRC
crc_encoder = comm.CRCGenerator('Polynomial', crc_polynomial);
crc_decoder = comm.CRCDetector('Polynomial', crc_polynomial);

%% Définition des paramètres et des échantillons
N0 = 88;                % Nombre de bits
Ts = 10^(-6);          % Durée symbole
Fs = 1/Ts;             % Fréquence d'échantillonnage
fe = 20 * 10^6;        % Fréquence d'échantillonnage réelle
Te = 1/fe;             % Période d'échantillonnage
Fse = Ts * fe;         % Facteur de sur-échantillonnage
b_k =  randi([0 1], 1, N0); % Bits générés (bits utiles sans CRC)
sigma_nl = 1;        % Carré de l'Écart-type du bruit

%% Ajout du codage CRC au message original
b_k_crc = crc_encoder(b_k');  

%% La modulation PPM
s_l = zeros(length(b_k_crc) * (Ts/Te), 1);
for k = 1:length(s_l)
    s_l(k) = PPM(k*Te, b_k_crc', Ts);  % Modulation du message avec CRC (transposition inversée)
end

%% Ajout du bruit n_l(t)
n_l = sigma_nl * randn(size(s_l));  % Génération du bruit
y_l = s_l + n_l;                    % Ajout du bruit au signal modulé

%% Filtrage adapté
g_a = ones(Ts/(2*Te), 1);           % Filtre adapté
r_l_1 = conv(y_l, g_a);             % Filtrage du signal reçu

%% Sur-échantillonnage avec Ts/2
r_l = zeros(length(b_k_crc) * 2, 1);
for k = 1:length(r_l)
    r_l(k) = r_l_1(10 * (k-1) + 10); 
end

%% Comparaison des bits originaux avec les bits estimés
vr_values = zeros(length(b_k_crc), 2); 
b_k_recus = zeros(1, length(b_k_crc)); 

for k = 1:length(b_k_crc)
    vr1 = norm([r_l(2*k - 1), (r_l(2*k) - 1)]);
    vr2 = norm([(r_l(2*k - 1) - 1), r_l(2*k)]);

    vr_values(k, 1) = vr1;
    vr_values(k, 2) = vr2;

    if vr1 < vr2
        b_k_recus(k) = 0;
    else
        b_k_recus(k) = 1;
    end
end

%% Détection d'erreurs avec le décodeur CRC
[b_k_decoded, has_errors] = crc_decoder(b_k_recus');

%% Comparaison des bits originaux avec les bits reçus après décodage
errors = sum(b_k_crc ~= b_k_recus');
disp(['Nombre d''erreurs: ', num2str(errors)]);

if has_errors
    disp('Erreur détectée lors de la réception.');
else
    disp('Aucune erreur détectée.');
end

%% Création d'une figure avec plusieurs subplots
t = (0:length(s_l)-1) * Te;  % Temps pour l'affichage du signal modulé
t1 = (0:length(r_l)-1) * Te;
figure;

% Bits générés avec CRC
subplot(5,1,1);
stem(b_k_crc, 'filled');
title('Bits générés avec CRC');
xlabel('Index du bit');
ylabel('Valeur du bit');
grid on;

% Signal modulé PPM
subplot(5,1,2);
plot(t, s_l);
title('Signal modulé PPM');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Signal avec bruit
subplot(5,1,3);
plot(t, y_l);
title('Signal avec bruit');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;


% Signal après filtrage adapté
subplot(5,1,4);
plot(t, r_l_1(1:length(t)));  
title('Signal après filtrage adapté');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Signal suréchantillonné
subplot(5,1,5);
stem(r_l);
title('Signal suréchantillonné');
xlabel('Index');
ylabel('Amplitude');
grid on;


%% Affichage des bits reçus après décision
figure;
stem(b_k_recus, 'filled');
title('Bits reçus après décision');
xlabel('Index du bit');
ylabel('Valeur du bit');
grid on;
