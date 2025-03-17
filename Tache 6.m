% Script principal pour le décodage ADS-B avec affichage de carte

% Chargement des données
load('adsb_msgs.mat');
disp(['Taille de adsb_msgs: ', num2str(size(adsb_msgs))]);

% Initialisation des vecteurs pour stocker les positions
latitudes = [];
longitudes = [];
altitudes = [];
timestamps = [];

% Coordonnées de référence (à ajuster selon votre emplacement)
refLat = 44.8378; % Latitude de Bordeaux
refLon = -0.5792; % Longitude de Bordeaux

% Traitement de chaque trame
for i = 1:size(adsb_msgs, 2)
    disp(['Traitement de la trame ', num2str(i)]);
    trame = adsb_msgs(:, i);
    registre = bit2registre(trame, refLon, refLat);
    
    if ~isempty(registre)
        if ~isempty(registre.latitude) && ~isempty(registre.longitude)
            disp(['  Latitude décodée: ', num2str(registre.latitude), ', Longitude décodée: ', num2str(registre.longitude)]);
            if ~isnan(registre.latitude) && ~isnan(registre.longitude)
                latitudes = [latitudes; registre.latitude];
                longitudes = [longitudes; registre.longitude];
                if ~isempty(registre.altitude)
                    altitudes = [altitudes; registre.altitude];
                    disp(['  Altitude: ', num2str(registre.altitude)]);
                else
                    altitudes = [altitudes; NaN];
                    disp('  Altitude non disponible');
                end
                timestamps = [timestamps; i];
            else
                disp('  Position invalide (NaN)');
            end
        elseif ~isempty(registre.nom)
            disp(['Identification d''avion: ', registre.nom]);
        else
            disp('  Aucune information de position ou d''identification');
        end
    else
        disp('Trame invalide ou non traitée');
    end
end

disp(['Nombre de positions valides: ', num2str(length(latitudes))]);

% Affichage des données brutes
disp('Premières 5 positions:');
disp([timestamps(1:min(5,end)), latitudes(1:min(5,end)), longitudes(1:min(5,end)), altitudes(1:min(5,end))]);

% Affichage de la carte avec la trajectoire
if ~isempty(latitudes) && ~isempty(longitudes)
    affiche_carte(refLon, refLat);
    hold on;
    plot(longitudes, latitudes, 'b.-');
    plot(longitudes(1), latitudes(1), 'go', 'MarkerSize', 10); % Point de départ
    plot(longitudes(end), latitudes(end), 'ro', 'MarkerSize', 10); % Point d'arrivée
    hold off;
    disp('Trajectoire affichée sur la carte');
else
    warning('Pas de données de position valides pour tracer la trajectoire');
end


% Fonctions

function registre = bit2registre(bitPacketCRC, refLon, refLat)
    registre = bit2registre_(bitPacketCRC, refLon, refLat);
end

function registre = bit2registre_(bitPacketCRC, refLon, refLat)
    % Vérification de la longueur du vecteur d'entrée
    if length(bitPacketCRC) ~= 112
        error('Le vecteur d''entrée doit contenir 112 bits.');
    end
    
    % S'assurer que bitPacketCRC est un vecteur ligne
    bitPacketCRC = bitPacketCRC(:)';
    
    % Initialisation de la structure registre
    registre = struct('adresse',[],'format',[],'type',[],'nom',[], ...
        'altitude',[],'timeFlag',[],'cprFlag',[],'latitude',[],'longitude',[]);
    
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
    elseif ismember(registre.type, [1,2,3,4])
        % Message d'identification
        idBits = bitPacketCRC(41:88);
        disp(['ID bits: ', num2str(idBits)]);
        registre.nom = decodeIdentification(idBits);
    else
        disp(['Type de message non traité: ', num2str(registre.type)]);
        registre = [];
        return;
    end
end

function valid = verifyCRC(bits)
    % Implémentation de la vérification CRC
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
    % Implémentation du calcul CRC-24
    data = [data, zeros(1, 24)];  % Ajout de 24 zéros à la fin des données
    
    for i = 1:length(data)-24
        if data(i)
            data(i:i+24) = xor(data(i:i+24), generator);
        end
    end
    
    crc = data(end-23:end);
end

function altitude = decodeAltitude(bits)
    % Suppression du 8ème bit
    bits = [bits(1:7) bits(9:end)];
    
    % Conversion en décimal
    alt = bin2dec(num2str(bits));
    
    % Calcul de l'altitude en pieds
    altitude = 25 * alt - 1000;
end

function nom = decodeIdentification(bits)
    nom = '';
    if length(bits) ~= 48
        error('Le vecteur d''identification doit contenir exactement 48 bits.');
    end
    
    disp(['Bits bruts d''identification: ', num2str(bits)]);
    for i = 1:8
        charBits = bits((i-1)*6+1 : i*6);
        charCode = bin2dec(num2str(charBits));
        disp(['Caractère ', num2str(i), ' code: ', num2str(charCode)]);
        nom = [nom, decodeChar(charCode)];
    end
    nom = strtrim(nom);
end


function character = decodeChar(code)
    % Table de correspondance basée sur l'Annexe D
    charTable = '#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######';
    
    disp(['Code reçu pour decodeChar: ', num2str(code)]);
    
    % Vérification que le code est dans la plage valide
    if isscalar(code) && isnumeric(code) && code >= 0 && code <= 63
        character = charTable(code + 1);
        disp(['Caractère décodé: ', character]);
    else
        character = '#'; % Caractère non valide
        disp('Code invalide, caractère # utilisé');
    end
end

function [lat, lon] = decodeCPR(latCPR, lonCPR, cprFlag, lat_ref, lon_ref)
    NZ = 15;
    
    % Conversion des valeurs CPR
    lat_cpr = latCPR / 2^17;
    lon_cpr = lonCPR / 2^17;
    
    % Calcul de la latitude
    Dlat = 360 / (4 * NZ - cprFlag);
    j = floor(lat_ref / Dlat) + floor(0.5 + mod(lat_ref, Dlat) / Dlat - lat_cpr);
    lat = Dlat * (j + lat_cpr);
    
    % Calcul de la longitude
    if lat >= 87 || lat <= -87
        Dlon = 360;
    else
        Dlon = 360 / max(1, cprNL(lat) - cprFlag);
    end
    m = floor(lon_ref / Dlon) + floor(0.5 + mod(lon_ref, Dlon) / Dlon - lon_cpr);
    lon = Dlon * (m + lon_cpr);
    
    % Vérification de la validité des coordonnées
    if lat < -90 || lat > 90 || lon < -180 || lon > 180
        lat = NaN;
        lon = NaN;
    end
end

function NL = cprNL(lat)
    % Implémentation de la fonction NL(x) de l'Annexe C
    if abs(lat) >= 87
        NL = 1;
    elseif lat == 0
        NL = 59;
    else
        NL = floor(2*pi / acos(1 - (1 - cos(pi/(2*15))) / cos(lat*pi/180)^2));
    end
end

function [] = affiche_carte(REF_LON, REF_LAT)
% Plot trajectoire + logo avion
STYLES = {'-','--',':'};
STYLES_HEAD = {'x','o','<'};
COLORS = lines(6);
COLORS(4,:)=[];
figure(1);
hold on;
%Bdx
x = linspace(-1.3581,0.7128,1024);
y = linspace(44.4542,45.1683,1024);
[X,Y] = meshgrid(x,y(end:-1:1));
im = imread('fond.jpg');
image(x,y(end:-1:1),im);
plot(REF_LON,REF_LAT,'.r','MarkerSize',20);
text(REF_LON+0.05,REF_LAT,0,'Actual pos','color','b')
set(gca,'YDir','normal')
xlabel('Longitude en degres');
ylabel('Lattitude en degres');
zlim([0,4e4]);
% Bdx
xlim([-1.3581,0.7128]);
ylim([44.4542,45.1683]);
end