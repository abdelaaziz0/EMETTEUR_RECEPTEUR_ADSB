# Simulation d’un émetteur / récepteur ADS-B en MATLAB

🚀 **Projet TS229 – Année 2024/2025**

---

## 📋 Description
Ce projet vise à simuler sous **MATLAB** un émetteur/récepteur ADS-B (Automatic Dependent Surveillance Broadcast), puis à adapter le récepteur afin d’effectuer un **décodage en temps réel** des avions survolant la zone d’expérimentation.

L'objectif principal est de mettre en pratique les concepts étudiés en **communications numériques, codage de canal et traitement du signal** afin de concevoir un système de réception et d'affichage des données ADS-B.

---

## 🔧 Prérequis
📌 **Logiciel** : MATLAB avec les Toolboxes **Signal Processing** et **Communications**

📌 **Connaissances** : Modulation, DSP, synchronisation en télécommunications

---

## 📦 Installation
```sh
# Cloner le dépôt
git clone https://github.com/abdelaaziz0/EMETTEUR_RECEPTEUR_ADSB.git
```

Lancer MATLAB et ajouter le répertoire du projet :
```matlab
addpath(genpath('ADSB_Simulation'))
```

---

## 🚀 Utilisation
### 1️⃣ Exécution de la simulation
Pour exécuter une simulation complète de la couche physique ADS-B :
```matlab
simulate_ADSB()
```

### 2️⃣ Décodage de données réelles
Pour utiliser le récepteur SDR en temps réel :
```matlab
decode_ADSB_real_time()
```

Le programme affiche alors les trajectoires des avions détectés sur une carte.

---

## 📝 Organisation du projet
### 📌 Modules principaux
Le projet est organisé en plusieurs tâches distinctes correspondant aux différentes couches du système ADS-B :

- **Tâche 1 : Couche physique ADS-B** (modulation PPM, simulation de la chaîne de communication)
- **Tâche 2 : Analyse spectrale** (densité spectrale de puissance via l’algorithme de Welch)
- **Tâche 3 : Codage et décodage de canal** (ajout et vérification du CRC)
- **Tâche 4 : Synchronisation temporelle** (alignement des signaux en réception)
- **Tâche 5 : Synchronisation fréquentielle** (correction des erreurs Doppler)
- **Tâche 6-7 : Implémentation de la couche MAC** (décodage des messages ADS-B)
- **Tâche 8-12 : Applications avancées** (traitement en temps réel, affichage des trajectoires, vue 3D)

---

## 🔍 Fonctionnement détaillé
### 1️⃣ Modulation PPM
L’ADS-B utilise une **modulation par position d’impulsion (PPM)** où chaque bit est encodé par un signal spécifique :
```matlab
[p0, p1] = generate_PPM_signals();
```

### 2️⃣ Traitement du signal et synchronisation
L’algorithme de synchronisation repose sur **l’intercorrélation avec le préambule** ADS-B :
```matlab
delta_t = estimate_sync_time(signal, preamble);
```

### 3️⃣ Décodage des trames ADS-B
Les trames ADS-B sont extraites et analysées pour récupérer les informations d’avion :
```matlab
registre = bit2registre(trame_ADSB);
```

### 4️⃣ Affichage en temps réel des trajectoires
L’interface graphique permet d’afficher en **temps réel** la position des avions détectés :
```matlab
plot_ADSB_trajectories(trajectoire_data);
```

---

## 📊 Tests et Validation
🛠 **Vérification de la chaîne de communication**
```matlab
test_modulation_PPM()
```

🛠 **Comparaison de la DSP théorique et expérimentale**
```matlab
compare_DSP_Welch()
```

🛠 **Vérification du décodage des trames ADS-B**
```matlab
test_decodage_ADSB()
```

---

## 📄 Licence
Ce projet est sous licence MIT. Voir le fichier `LICENSE` pour plus de détails.
