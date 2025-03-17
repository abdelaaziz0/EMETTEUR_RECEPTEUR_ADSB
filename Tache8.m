% Tâche 8 - Traitement de signaux réels ADS-B complexes avec distorsions

% Chargement des données
load('buffers.mat');
disp(['Nombre de buffers : ', num2str(size(buffers, 2))]);
disp(['Taille de chaque buffer : ', num2str(size(buffers, 1))]);

% Paramètres
Fe = 4e6;  % Fréquence d'échantillonnage (4 MHz)
Ts = 1e-6; % Durée symbole (1 µs)
Fse = Ts * Fe; % Facteur de sur-échantillonnage
disp(['Facteur de sur-échantillonnage : ', num2str(Fse)]);

% Coordonnées de référence (à ajuster selon votre emplacement)
refLat = 44.8378; % Latitude de Bordeaux
refLon = -0.5792; % Longitude de Bordeaux
disp(['Coordonnées de référence : Lat ', num2str(refLat), ', Lon ', num2str(refLon)]);

% Initialisation des structures pour stocker les données des avions
adresses = {};
avions = struct('positions', {}, 'timestamps', {}, 'altitudes', {}, 'noms', {});
disp('Structures des avions initialisées.');

% Traitement de chaque buffer
for buffer_idx = 1:size(buffers, 2)
    % Chargement du buffer et traitement comme s_l_distorted
    s_l_distorted = buffers(:, buffer_idx);
    disp(['Traitement du buffer numéro : ', num2str(buffer_idx)]);
    
    % Synchronisation et correction des distorsions avec s_l_distorted
    [signal_sync, delta_t_est] = synchronize_signal(s_l_distorted, generate_preamble(), 1/Fe);
    disp(['Délai estimé : ', num2str(delta_t_est)]);
    
    % Démodulation PPM avec le signal synchronisé
    bits = demodulation_ppm(signal_sync, Fe, length(generate_preamble()));
    disp(['Nombre de bits démodulés : ', num2str(length(bits))]);

    % Vérification CRC et décodage
    if verifyCRC(bits)
        registre = bit2registre(bits, refLon, refLat);
        
        if ~isempty(registre) && ~isempty(registre.adresse)
            adresse = num2str(registre.adresse);
            disp(['Adresse de l\''avion : ', adresse]);
            
            % Recherche de l'index de l'avion dans notre liste
            idx = find(strcmp(adresses, adresse));
            
            if isempty(idx)
                % Nouvel avion, on l'ajoute à notre liste
                adresses{end+1} = adresse;
                idx = length(adresses);
                avions(idx).positions = [];
                avions(idx).timestamps = [];
                avions(idx).altitudes = [];
                avions(idx).noms = {};
                disp(['Nouvel avion ajouté : ', adresse]);
            end
            
            % Ajout des nouvelles données
            if ~isempty(registre.latitude) && ~isempty(registre.longitude)
                avions(idx).positions(end+1,:) = [registre.latitude, registre.longitude];
                avions(idx).timestamps(end+1) = buffer_idx;
                if isfield(registre, 'altitude')
                    avions(idx).altitudes(end+1) = registre.altitude;
                else
                    avions(idx).altitudes(end+1) = NaN;
                end
                disp(['Position : ', num2str(registre.latitude), ', ', num2str(registre.longitude), ...
                      ' Altitude : ', num2str(registre.altitude), ' pieds']);
            end
            
            if ~isempty(registre.nom)
                avions(idx).noms{end+1} = registre.nom;
                disp(['Identification : ', registre.nom]);
            end
        end
    else
        disp('CRC non valide pour ce buffer.');
    end
end

% Affichage des résultats
disp('Résumé des avions détectés :');
for i = 1:length(adresses)
    avion = avions(i);
    disp(['Avion ', adresses{i}]);
    if ~isempty(avion.noms)
        disp(['  Identification : ', avion.noms{end}]);
    end
    disp(['  Nombre de positions : ', num2str(size(avion.positions, 1))]);
    if ~isempty(avion.altitudes)
        disp(['  Altitude moyenne : ', num2str(mean(avion.altitudes, 'omitnan')), ' pieds']);
    end
end

% Affichage de la carte avec les trajectoires
figure;
hold on;
for i = 1:length(adresses)
    avion = avions(i);
    if ~isempty(avion.positions)
        plot(avion.positions(:,2), avion.positions(:,1), '.-');
        disp(['Trajectoire de l\''avion ', adresses{i}, ' tracée.']);
    end
end
affiche_carte(refLon, refLat);
xlabel('Longitude');
ylabel('Latitude');
title('Trajectoires des avions');
legend(adresses);
hold off;
disp('Affichage de la carte terminé.');


% Fonctions auxiliaires

function preamble = generate_preamble()
    % Génération du préambule ADS-B
    Tp = 8e-6;               % Durée du préambule (8 µs)
    Fe = 4e6;                % Fréquence d'échantillonnage
    Te = 1/Fe;               % Période d'échantillonnage
    unite = 0.5e-6;          % Durée d'une unité (0.5 µs)
    s_p = zeros(1, floor(Tp/Te));
    for k = 1:length(s_p)
        if (k >= 1 && k < unite/Te) || (k >= 2*unite/Te && k < 3*unite/Te) || ...
           (k >= 7*unite/Te && k < 8*unite/Te) || (k >= 9*unite/Te && k < 10*unite/Te)
            s_p(k) = 1;
        end
    end
    preamble = s_p;
end

function [y_l_sync, delta_t_est] = synchronize_signal(y_l, s_p, Te)
    preamble_length = length(s_p);
    max_delay = 200;
    
    correlation = zeros(1, length(y_l) - preamble_length + 1);
    s_p_col = s_p(:);
    
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
    
    [~, max_index] = max(correlation);
    delta_t_est = (max_index - 1) * Te;
    start_index = max_index;
    
    if start_index < 1
        start_index = 1;
    elseif start_index > length(y_l) - preamble_length
        start_index = length(y_l) - preamble_length;
    end
    
    y_l_sync = y_l(start_index:end);
    y_l_sync = y_l_sync / max(abs(y_l_sync));
    
    expected_length = 112 * 20;  % 112 bits de données * 20 (Fse)
    if length(y_l_sync) > expected_length
        y_l_sync = y_l_sync(1:expected_length);
    elseif length(y_l_sync) < expected_length
        y_l_sync = [y_l_sync; zeros(expected_length - length(y_l_sync), 1)];
    end
end

function bits = demodulation_ppm(signal, Fe, preamble_length)
    % signal          : Signal PPM reçu
    % Fe              : Fréquence d'échantillonnage
    % preamble_length : Longueur du préambule en échantillons à enlever

    Ts = 1e-6;  % Durée du symbole en secondes (ajuster si nécessaire)
    samples_per_symbol = round(Fe * Ts);  % Nombre d'échantillons par symbole
    
    % Suppression du préambule
    signal_without_preamble = signal(preamble_length+1:end);
    
    % On ne garde que les 112 premiers symboles
    num_symbols = 112;  
    bits = zeros(1, num_symbols);  % Initialiser les bits
    
    for i = 1:num_symbols
        symbol_start = (i-1) * samples_per_symbol + 1;
        symbol_end = i * samples_per_symbol;
        
        % Extraire le symbole correspondant
        symbol = signal_without_preamble(symbol_start:symbol_end);
        
        % Comparer l'énergie dans les deux moitiés du symbole
        first_half = sum(abs(symbol(1:samples_per_symbol/2)).^2);
        second_half = sum(abs(symbol(samples_per_symbol/2+1:end)).^2);
        
        % Si l'énergie dans la seconde moitié est plus grande, on considère un '1', sinon un '0'
        bits(i) = second_half > first_half;
    end
end

% Fonctions auxiliaires

function registre = bit2registre(bitPacketCRC, refLon, refLat)
    % Vérification de la longueur du vecteur d'entrée
    if length(bitPacketCRC) ~= 112
        error('Le vecteur d''entrée doit contenir 112 bits.');
    end
    
    % S'assurer que bitPacketCRC est un vecteur ligne
    bitPacketCRC = bitPacketCRC(:)';
    
    % Initialisation de la structure registre
    registre = struct('adresse',[],'format',[],'type',[],'nom',[], ...
        'altitude',[],'timeFlag',[],'cprFlag',[],'latitude',[],'longitude',[], ...
        'vitesse',[], 'direction',[], 'position_type',[]);
    
    % Vérification du CRC
    if ~verifyCRC(bitPacketCRC)
        disp('Erreur CRC détectée');
        registre = [];
        return;
    end
    
    % Extraction du format de la voie descendante (DF)
    registre.format = bin2dec(num2str(bitPacketCRC(1:5)));
    
    % Vérification que DF = 17 (message ADS-B)
    if registre.format ~= 17
        disp(['Format non ADS-B détecté: ', num2str(registre.format)]);
        registre = [];
        return;
    end
    
    % Extraction de l'adresse OACI
    registre.adresse = bin2dec(num2str(bitPacketCRC(9:32)));
    
    % Extraction du type de format (FTC)
    registre.type = bin2dec(num2str(bitPacketCRC(33:37)));
    
    % Traitement selon le type de message
    if ismember(registre.type, [9,10,11,12,13,14,15,16,17,18])
        % Message de position en vol
        registre.altitude = decodeAltitude(bitPacketCRC(41:52));
        registre.timeFlag = bitPacketCRC(53);
        registre.cprFlag = bitPacketCRC(54);
        latCPR = bin2dec(num2str(bitPacketCRC(55:71)));
        lonCPR = bin2dec(num2str(bitPacketCRC(72:88)));
        [registre.latitude, registre.longitude] = decodeCPR(latCPR, lonCPR, registre.cprFlag, refLat, refLon);
        registre.position_type = 'vol';
    elseif ismember(registre.type, [5,6,7,8])
        % Message de position au sol
        registre.mouvement = bin2dec(num2str(bitPacketCRC(38:44)));
        registre.statut = bitPacketCRC(45);
        registre.cprFlag = bitPacketCRC(54);
        latCPR = bin2dec(num2str(bitPacketCRC(55:71)));
        lonCPR = bin2dec(num2str(bitPacketCRC(72:88)));
        [registre.latitude, registre.longitude] = decodeCPR(latCPR, lonCPR, registre.cprFlag, refLat, refLon);
        registre.position_type = 'sol';
    elseif ismember(registre.type, [1,2,3,4])
        % Message d'identification
        registre.nom = decodeIdentification(bitPacketCRC(41:88));
    elseif registre.type == 19
        % Message de vitesse
        subtype = bin2dec(num2str(bitPacketCRC(38:40)));
        if subtype == 1 || subtype == 2
            % Vitesse au sol
            ew_dir = bitPacketCRC(46);
            ew_velocity = bin2dec(num2str(bitPacketCRC(47:56)));
            ns_dir = bitPacketCRC(57);
            ns_velocity = bin2dec(num2str(bitPacketCRC(58:67)));
            
            if ew_dir == 1
                ew_velocity = -ew_velocity;
            end
            if ns_dir == 1
                ns_velocity = -ns_velocity;
            end
            
            registre.vitesse = sqrt(ew_velocity^2 + ns_velocity^2);
            registre.direction = mod(atan2(ew_velocity, ns_velocity) * 180/pi, 360);
        elseif subtype == 3 || subtype == 4
            % Vitesse air
            registre.vitesse = bin2dec(num2str(bitPacketCRC(47:56)));
            heading_status = bitPacketCRC(57);
            if heading_status == 1
                registre.direction = bin2dec(num2str(bitPacketCRC(58:67))) * 360 / 1024;
            end
        end
    else
        disp(['Type de message non traité: ', num2str(registre.type)]);
        registre = [];
        return;
    end
end

function valid = verifyCRC(bits)
    % Polynôme générateur du CRC-24 pour ADS-B
    generator = [1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1];
    
    % Extraction des bits de données et du CRC reçu
    data = bits(1:88);
    receivedCRC = bits(89:end);
    
    % Calcul du CRC
    crc = crc24(data, generator);
    
    % Comparaison du CRC calculé avec le CRC reçu
    valid = isequal(crc, receivedCRC);
end

function crc = crc24(data, generator)
    data = [data, zeros(1, 24)];  % Ajout de 24 zéros à la fin des données
    
    for i = 1:length(data)-24
        if data(i)
            data(i:i+24) = xor(data(i:i+24), generator);
        end
    end
    
    crc = data(end-23:end);
end

function altitude = decodeAltitude(bits)
    bits = [bits(1:7) bits(9:end)];
    alt = bin2dec(num2str(bits));
    altitude = 25 * alt - 1000;
end

function nom = decodeIdentification(bits)
    nom = '';
    for i = 1:8
        charBits = bits((i-1)*6+1 : i*6);
        charCode = bin2dec(num2str(charBits));
        nom = [nom, decodeChar(charCode)];
    end
    nom = strtrim(nom);
end

function character = decodeChar(code)
    charTable = '#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######';
    if isscalar(code) && isnumeric(code) && code >= 0 && code <= 63
        character = charTable(code + 1);
    else
        character = '#';
    end
end

function [lat, lon] = decodeCPR(latCPR, lonCPR, cprFlag, lat_ref, lon_ref)
    NZ = 15;
    
    lat_cpr = latCPR / 2^17;
    lon_cpr = lonCPR / 2^17;
    
    Dlat = 360 / (4 * NZ - cprFlag);
    j = floor(lat_ref / Dlat) + floor(0.5 + mod(lat_ref, Dlat) / Dlat - lat_cpr);
    lat = Dlat * (j + lat_cpr);
    
    if lat >= 87 || lat <= -87
        Dlon = 360;
    else
        Dlon = 360 / max(1, cprNL(lat) - cprFlag);
    end
    m = floor(lon_ref / Dlon) + floor(0.5 + mod(lon_ref, Dlon) / Dlon - lon_cpr);
    lon = Dlon * (m + lon_cpr);
    
    if lat < -90 || lat > 90 || lon < -180 || lon > 180
        lat = NaN;
        lon = NaN;
    end
end

function NL = cprNL(lat)
    if abs(lat) >= 87
        NL = 1;
    elseif lat == 0
        NL = 59;
    else
        NL = floor(2*pi / acos(1 - (1 - cos(pi/(2*15))) / cos(lat*pi/180)^2));
    end
end
function [vitesse, direction] = calculVitesse(pos1, pos2, t1, t2)
    % Conversion des positions en radians
    lat1 = deg2rad(pos1(1));
    lon1 = deg2rad(pos1(2));
    lat2 = deg2rad(pos2(1));
    lon2 = deg2rad(pos2(2));
    
    % Rayon moyen de la Terre en km
    R = 6371;
    
    % Calcul de la distance orthodromique
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distance = R * c;
    
    % Calcul du temps écoulé en heures
    temps = (t2 - t1) / 3600; % Supposons que t1 et t2 sont en secondes
    
    % Calcul de la vitesse en km/h
    vitesse = distance / temps;
    
    % Conversion en nœuds (1 km/h = 0.539957 nœuds)
    vitesse = vitesse * 0.539957;
    
    % Calcul de la direction
    direction = atan2(sin(lon2-lon1)*cos(lat2), ...
                      cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1));
    direction = mod(rad2deg(direction), 360);
end
function affiche_carte(REF_LON, REF_LAT)
    % Définition des limites de la carte (ajustez selon vos besoins)
    lon_min = -1.3581;
    lon_max = 0.7128;
    lat_min = 44.4542;
    lat_max = 45.1683;
    
    
    % Affichage du point de référence
    plot(REF_LON, REF_LAT, 'r*', 'MarkerSize', 10);
    text(REF_LON + 0.05, REF_LAT, 'Référence', 'Color', 'r');
    
    % Définition des limites de l'axe
    xlim([lon_min, lon_max]);
    ylim([lat_min, lat_max]);
    
    % Ajout d'une grille
    grid on;
    
    
    % Personnalisation de l'apparence
    set(gca, 'FontSize', 10);
    box on;
end
