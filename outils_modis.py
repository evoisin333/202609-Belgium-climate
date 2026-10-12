# -*- coding: utf-8 -*-
"""
outils_modis.py
Fonctions communes a 03_modis_ndvi.py et 20_phenologie_ndvi.py : appariement
des fichiers AppEEARS, masque qualite et date reelle de chaque pixel.

Source : Didan, K., Barreto Munoz, A. (2019). MODIS Vegetation Index User's
Guide (MOD13 Series), Version 3.10, Collection 6.1. NASA LP DAAC.
  - pixel_reliability : 0 bon, 1 marginal ("useful, but look at other QA
    information"), 2 neige/glace, 3 nuageux, -1 non calcule ;
  - VI_Quality (16 bits, commune au NDVI et a l'EVI) : bits 2-5 indice
    d'utilite (0000 meilleur ... 1100 plus bas niveau utile, 1101 "too low
    to be useful", 1110-1111 inutilisable), bit 8 nuage adjacent, bit 10
    nuages melanges, bit 14 neige/glace possible, bit 15 ombre possible ;
  - composite_day_of_the_year : jour d'acquisition du pixel (1-366, -1 = vide).
"""

import glob
import os
import re

import numpy as np
import pandas as pd
import rasterio

COUCHES_QA = ("pixel_reliability", "VI_Quality", "composite_day_of_the_year")
VALEUR_VIDE = -3000   # valeur de remplissage du NDVI et de l'EVI
ECHELLE = 0.0001      # facteur d'echelle du NDVI et de l'EVI
UTILITE_MAX = 0b1100  # au-dela : "too low to be useful" selon le guide


def debut_du_fichier(chemin):
    """(annee, jour julien) du premier jour du composite, lus dans le nom."""
    m = re.search(r"doy(\d{4})(\d{3})", os.path.basename(chemin))
    return int(m.group(1)), int(m.group(2))


def fichiers_apparies(dossier_modis, indice, annees):
    """Composites presents dans les 4 dossiers (indice + 3 couches QA),
    dont l'annee de debut est dans annees (tuple debut, fin) ;
    renvoie une liste triee de (annee, jour, {couche: chemin})."""
    couches = (indice,) + COUCHES_QA
    par_couche = {}
    for c in couches:
        fichiers = glob.glob(os.path.join(dossier_modis, c, "*.tif"))
        if not fichiers:
            raise SystemExit("Aucun fichier dans %s" % os.path.join(dossier_modis, c))
        par_couche[c] = {debut_du_fichier(f): f for f in fichiers}
    communs = set.intersection(*(set(d) for d in par_couche.values()))
    for c in couches:
        manquants = len(par_couche[c]) - len(communs)
        if manquants:
            print("  !! %s : %d composite(s) sans equivalent dans les autres couches, ignores"
                  % (c, manquants))
    return [(a, j, {c: par_couche[c][(a, j)] for c in couches})
            for (a, j) in sorted(communs) if annees[0] <= a <= annees[1]]


def lire(chemin):
    with rasterio.open(chemin) as src:
        return src.read(1)


def masque_qualite(valeur, fiabilite, vi_quality, details=None):
    """Pixels gardes : fiabilite 0, ou fiabilite 1 si VI_Quality ne signale
    ni nuage (adjacent ou melange), ni ombre, ni neige, et si l'indice
    d'utilite n'est pas 'trop bas'. details (dict) : compteurs optionnels
    des pixels marginaux retires par chaque critere (non exclusifs)."""
    q = vi_quality.astype(np.uint16)          # lecture des bits sans signe
    rejets = {
        "utilite trop basse": ((q >> 2) & 0b1111) > UTILITE_MAX,
        "nuage adjacent": ((q >> 8) & 1) == 1,
        "nuages melanges": ((q >> 10) & 1) == 1,
        "neige possible": ((q >> 14) & 1) == 1,
        "ombre possible": ((q >> 15) & 1) == 1,
    }
    marginal = (fiabilite == 1) & (valeur != VALEUR_VIDE)
    rejete = np.zeros(valeur.shape, dtype=bool)
    for nom, r in rejets.items():
        rejete |= r
        if details is not None:
            details[nom] = details.get(nom, 0) + int(np.sum(marginal & r))
    if details is not None:
        details["marginaux"] = details.get("marginaux", 0) + int(marginal.sum())
        details["marginaux gardes"] = details.get("marginaux gardes", 0) + int(np.sum(marginal & ~rejete))
    return (valeur != VALEUR_VIDE) & ((fiabilite == 0) | ((fiabilite == 1) & ~rejete))


_MOIS = {}


def date_des_pixels(annee_debut, jour_debut, jour_pixel):
    """Annee et mois d'acquisition de chaque pixel. Un composite qui
    commence en fin d'annee peut contenir des pixels de janvier suivant :
    leur jour julien est alors plus petit que le jour de debut."""
    annee = np.where(jour_pixel < jour_debut, annee_debut + 1, annee_debut)
    mois = np.zeros(jour_pixel.shape, dtype=np.int8)     # 0 = pas de date
    for a in np.unique(annee):
        if a not in _MOIS:
            _MOIS[a] = pd.date_range("%d-01-01" % a, "%d-12-31" % a).month.values
        table = _MOIS[a]
        sel = (annee == a) & (jour_pixel >= 1) & (jour_pixel <= len(table))
        mois[sel] = table[jour_pixel[sel] - 1]
    return annee, mois
