% Script principal pour le décodage ADS-B amélioré

% Chargement des données
load('adsb_msgs.mat');
disp(['Taille de adsb_msgs: ', num2str(size(adsb_msgs))]);

% Coordonnées de référence (à ajuster selon votre emplacement)
refLat = 44.8378; % Latitude de Bordeaux
refLon = -0.5792; % Longitude de Bordeaux

% Initialisation des structures pour stocker les données des avions
adresses = {};
avions = struct('positions', {}, 'timestamps', {}, 'vitesses', {}, 'directions', {}, 'altitudes', {}, 'noms', {}, 'trames_paires', {}, 'trames_impaires', {});

% Traitement de chaque trame
for i = 1:size(adsb_msgs, 2)
    trame = adsb_msgs(:, i);
    registre = bit2registre(trame, refLon, refLat);
    
    if ~isempty(registre) && ~isempty(registre.adresse)
        adresse = num2str(registre.adresse);
        
        % Recherche de l'index de l'avion dans notre liste
        idx = find(strcmp(adresses, adresse));
        
        if isempty(idx)
            % Nouvel avion, on l'ajoute à notre liste
            adresses{end+1} = adresse;
            idx = length(adresses);
            avions(idx).positions = [];
            avions(idx).timestamps = [];
            avions(idx).vitesses = [];
            avions(idx).directions = [];
            avions(idx).altitudes = [];
            avions(idx).noms = {};
            avions(idx).trames_paires = [];
            avions(idx).trames_impaires = [];
        end
        
        % Traitement selon le type de message
        if ~isempty(registre.latitude) && ~isempty(registre.longitude)
            % Position (en vol ou au sol)
            avions(idx).positions(end+1,:) = [registre.latitude, registre.longitude];
            avions(idx).timestamps(end+1) = i;
            if isfield(registre, 'altitude')
                avions(idx).altitudes(end+1) = registre.altitude;
            end
            if isfield(registre, 'position_type') && strcmp(registre.position_type, 'sol')
                disp(['Avion ', adresse, ' - Position au sol: ', num2str(registre.latitude), ', ', num2str(registre.longitude)]);
            else
                disp(['Avion ', adresse, ' - Position en vol: ', num2str(registre.latitude), ', ', num2str(registre.longitude), ' à ', num2str(registre.altitude), ' pieds']);
            end
        end
        
        if isfield(registre, 'vitesse') && ~isempty(registre.vitesse)
            % Vitesse
            avions(idx).vitesses(end+1) = registre.vitesse;
            if isfield(registre, 'direction') && ~isempty(registre.direction)
                avions(idx).directions(end+1) = registre.direction;
                disp(['Avion ', adresse, ' - Vitesse: ', num2str(registre.vitesse), ' noeuds, Direction: ', num2str(registre.direction), '°']);
            else
                disp(['Avion ', adresse, ' - Vitesse: ', num2str(registre.vitesse), ' noeuds']);
            end
        end
        
        if ~isempty(registre.nom)
            % Identification
            avions(idx).noms{end+1} = registre.nom;
            disp(['Avion ', adresse, ' - Identification: ', registre.nom]);
        end
        
        % Stockage des trames pour le décodage CPR amélioré
        if isfield(registre, 'cprFlag')
            if registre.cprFlag == 0
                avions(idx).trames_paires(end+1,:) = trame;
            else
                avions(idx).trames_impaires(end+1,:) = trame;
            end
        end
    end
end

% Affichage des résultats
for i = 1:length(adresses)
    avion = avions(i);
    disp(['Résumé pour l''avion ', adresses{i}]);
    if ~isempty(avion.noms)
        disp(['  Identification: ', avion.noms{end}]);
    end
    disp(['  Nombre de positions: ', num2str(size(avion.positions, 1))]);
    if ~isempty(avion.vitesses)
        disp(['  Vitesses: ', num2str(avion.vitesses)]);
        disp(['  Vitesse moyenne: ', num2str(mean(avion.vitesses)), ' noeuds']);
    else
        disp('  Aucune information de vitesse disponible');
    end
    if ~isempty(avion.altitudes)
        disp(['  Altitude moyenne: ', num2str(mean(avion.altitudes)), ' pieds']);
    end
end

% Affichage de la carte avec les trajectoires
figure;
hold on;
for i = 1:length(adresses)
    avion = avions(i);
    if ~isempty(avion.positions)
        plot(avion.positions(:,2), avion.positions(:,1), '.-');
    end
end
xlabel('Longitude');
ylabel('Latitude');
title('Trajectoires des avions');
hold off;

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

