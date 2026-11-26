# config.py — paramètres du projet
# C'est le SEUL fichier à modifier pour changer de ville ou de période.
# Tous les scripts lisent leurs paramètres ici.

from pathlib import Path

# Dossier racine du dépôt (le dossier qui contient ce fichier)
RACINE = Path(__file__).resolve().parent

# --- Ville étudiée ---
VILLE = "bruxelles"

# Emprise de recherche en WGS84 : (lon_min, lat_min, lon_max, lat_max)
# Volontairement plus large que la Région (marge de quelques km) :
# la limite exacte viendra de la couche régionale (audit, étape 1.7).
EMPRISE_WGS84 = (4.20, 50.74, 4.52, 50.94)
CRS_CALCUL = "EPSG:31370"   # Lambert 72 : CRS de calcul des vecteurs (phase 2, étape 2.1)

# --- Période ---
ANNEES = list(range(2020, 2026))        # analyse : 2020 à 2025 inclus
ANNEES_AUDIT = list(range(2013, 2026))  # audit Landsat : 2013 à 2025 (repli possible)
MOIS_ETE = [6, 7, 8]                    # juin, juillet, août

# --- Dossiers (chemins relatifs à la racine du dépôt) ---
DATA_BRUT = RACINE / "data" / "brut" / VILLE        # téléchargements, jamais modifiés
DATA_TRAVAIL = RACINE / "data" / "travail" / VILLE  # extractions légères, fichiers intermédiaires
AUDIT = RACINE / "data" / "travail" / VILLE / "audit"
RESULTATS = RACINE / "05_resultats"                 # CSV qui font foi
FIGURES = RACINE / "06_figures"

# --- Espaces verts accessibles : filtre principal (phase 2, étape 2.2) ---
# Définition « zones de carence » de Bruxelles Environnement : accessible de droit,
# plus de 50 % de couverture végétale, aménagé pour le public.
EV_STATUT = ["Public", "Public/privé"]  # « espace à caractère public » (Plan Nature 2016), pas la propriété du sol
EV_PCT_VEG_MIN = 50                     # couverture strictement supérieure (« plus de 50 % »)
EV_EAU_COMPTEE = True                   # couverture = pct_veg + pct_water (IBSA, Focus n° 63)
EV_SOUS_TYPES = ["Bois et forêt", "Espace vert aménagé", "Cimetière", "Abords de voirie"]

# --- Réseau piéton (phase 2, étape 2.3) ---
SEUIL_ENTREE_DIFFUSE_M = 10   # entrée diffuse : voie à moins de 10 m d'un espace vert sans entrée ponctuelle
ENTREE_PENETRATION_M = 5      # entrée ponctuelle : la voie atteint l'intérieur de l'espace rétréci de 5 m
EV_FERMES_PGS_ID = ["11006"]   # bois du CERIA (Anderlecht) : fermé, vérifié sur Street View le 24/11/2026

# --- Population (phase 2, étape 2.4) ---
SEUIL_POP_SECTEUR = 50        # secteurs de moins de 50 hab. exclus des statistiques (gardés dans les totaux)