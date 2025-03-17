%% Définition des paramètres
N_b = 1000; % Taille du paquet d'information binaire
Ts = 10^(-6); % Durée symbole
fe = 20 * 10^6; % Fréquence d'échantillonnage réelle
Te = 1/fe; % Période d'échantillonnage
Fse = Ts * fe; % Facteur de sur-échantillonnage
Eb_N0_dB = 0:1:10; % Plage de Eb/N0 en dB
BER = zeros(size(Eb_N0_dB)); % Préallocation du TEB

%% Boucle sur chaque valeur de Eb/N0
for i = 1:length(Eb_N0_dB)
    Eb_N0 = 10^(Eb_N0_dB(i)/10); % Conversion Eb/N0 en linéaire
    errors_total = 0;
    bit_total = 0;
    
    % Continuer jusqu'à avoir au moins 100 erreurs ou atteindre la limite de bits simulés
    while errors_total < 100 
        % Génération des bits et du signal PPM pour cette itération
        b_k = randi([0 1], 1, N_b);
        s_l = zeros(N_b * Fse, 1);
        for k = 1:length(s_l)
            s_l(k) = PPM(k * Te, b_k, Ts);
        end
        
        % Ajout du bruit
        sigma_nl = sqrt(mean(s_l.^2)*Fse / (2 * Eb_N0));
        n_l = sigma_nl * randn(size(s_l));
        y_l = s_l + n_l;
        
        % Filtrage adapté
        g_a = ones(Ts/(2*Te), 1);
        r_l_1 = conv(y_l, g_a);
        
        % Sur-échantillonnage avec Ts/2
        r_l = zeros(length(b_k) * 2, 1);
        for k = 1:length(r_l)
            r_l(k) = r_l_1(10 *(k-1) +10); 
        end
        
        % Comparaison des bits originaux avec les bits estimés
        b_k_recus = zeros(1, N_b);
        for k = 1:N_b
            vr1 = norm([r_l(2*k - 1), (r_l(2*k) - 1)]);
            vr2 = norm([(r_l(2*k - 1) - 1), r_l(2*k)]);
            b_k_recus(k) = vr1 > vr2;
        end
        
        % Calcul du nombre d'erreurs pour cette séquence
        errors = sum(b_k ~= b_k_recus);
        errors_total = errors_total + errors;
        bit_total = bit_total + N_b;
    end
    
    % Calcul du TEB pour cette valeur de Eb/N0
    BER(i) = errors_total / bit_total;
    
    % Affichage des résultats intermédiaires
    disp(['Eb/N0 = ', num2str(Eb_N0_dB(i)), ' dB, Bits totaux simulés = ', num2str(bit_total), ', Erreurs = ', num2str(errors_total)]);
end

%% Affichage des résultats simulés
figure;
semilogy(Eb_N0_dB, BER, 'o-');
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Taux d''Erreur Binaire (BER)');
title('Performance de la modulation PPM');

% Calcul et tracé de la courbe théorique
Pb_theory = theoretical_pb(10.^(Eb_N0_dB/10));
hold on;
semilogy(Eb_N0_dB, Pb_theory, 'r--');
legend('Simulation', 'Théorique');

%% Fonction pour calculer la probabilité d'erreur théorique Pb pour PPM binaire
function Pb = theoretical_pb(Eb_N0)
    Pb = qfunc(sqrt(2*Eb_N0));
end
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
