# -*- coding: utf-8 -*-
"""
20_phenologie_ndvi.py
A lancer depuis la racine du projet, environnement geo :
    python 20_phenologie_ndvi.py          (NDVI)
    python 20_phenologie_ndvi.py EVI      (EVI, meme analyse)

La baisse du NDVI moyen d'ete sur les cultures vient-elle d'une vegetation
moins vigoureuse, ou d'un calendrier qui avance (senescence et recolte plus
precoces) ? On regarde la tendance mois par mois :
  - si le printemps et le maximum ne baissent pas mais que juillet-aout
    baissent : plutot un decalage de calendrier ;
  - si tous les mois et le maximum baissent : plutot une perte de vigueur.

Pour chaque composite MODIS de 16 jours (2001-2024) :
  - meme masque qualite que 03_modis_ndvi.py (outils_modis.py : fiabilite 0,
    ou 1 si VI_Quality ne signale ni nuage, ni ombre, ni neige ; valeur
    -3000 retiree, facteur 0,0001) ; NDVI et EVI partagent ces couches ;
  - chaque pixel est date par son jour d'acquisition reel
    (composite_day_of_the_year) : un composite a cheval sur deux mois
    alimente les deux ;
  - composite ignore pour une classe si moins de 50 % de ses pixels sont valides.
Indicateurs par classe et par annee :
  - NDVI moyen de chaque mois, mars a novembre (9 indicateurs), moyenne de
    tous les pixels valides de la classe acquis ce mois-la ;
  - NDVI moyen de chaque saison : MAM, JJA, SON (3 indicateurs), seulement
    si au moins 3 composites valides y contribuent (comme 03) ;
  - NDVI maximum de mars a novembre (1 indicateur) : plus forte moyenne de
    composite, chaque composite etant rattache au mois de la plupart de ses
    pixels valides.
Tendance de chaque serie : pente de Sen, Mann-Kendall 3PW, puis correction
FDR 10 % sur l'ensemble des tests (outils_tendance.py). Au moins 18 annees.

Le JJA calcule ici (moyenne spatiale de chaque composite, puis moyenne de la
saison) peut differer legerement de celui de 16_analyse.py (moyenne de la
saison par pixel, puis moyenne spatiale) : l'ordre des moyennes change.

L'EVI sert de controle pour les forets : le NDVI sature quand le couvert est
dense (Huete et al. 2002), l'EVI beaucoup moins.

Entrees :
  01_raw/MODIS/NDVI/*.tif (ou EVI/), et les couches qualite
  01_raw/MODIS/pixel_reliability, VI_Quality, composite_day_of_the_year
  01_raw/corine_8classes_grille_modis.tif   (15_ndvi_corine.py)
  04_stats/corine_8classes_libelles.csv
Sorties :
  05_resultats/<indice>_indicateurs_par_classe.csv   series annuelles
  05_resultats/phenologie_<indice>_par_classe.csv    tendances
  06_figures/fig5_tendance_ndvi_mensuelle (fig5b_tendance_evi_mensuelle pour l'EVI)
  en PNG et PDF, et la version anglaise dans 06_figures_en
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray as rxr

from outils_modis import ECHELLE, date_des_pixels, fichiers_apparies, lire, masque_qualite
from outils_tendance import fdr, p_3pw, tendance_serie

RACINE = os.getcwd()   # racine du projet
INDICE = sys.argv[1].upper() if len(sys.argv) > 1 else "NDVI"
if INDICE not in ("NDVI", "EVI"):
    raise SystemExit("Indice inconnu : %s (NDVI ou EVI)" % INDICE)
MODIS = os.path.join(RACINE, "01_raw", "MODIS")
SUFFIXE = INDICE.lower()
NOM_FIG = "fig5_tendance_ndvi_mensuelle" if INDICE == "NDVI" else "fig5b_tendance_evi_mensuelle"
print("indice : %s" % INDICE)
CLASSES_TIF = os.path.join(RACINE, "01_raw", "corine_8classes_grille_modis.tif")
LIBELLES_CSV = os.path.join(RACINE, "04_stats", "corine_8classes_libelles.csv")
RESULTATS = os.path.join(RACINE, "05_resultats")
os.makedirs(RESULTATS, exist_ok=True)

ANNEES = (2001, 2024)
MOIS = [3, 4, 5, 6, 7, 8, 9, 10, 11]
SAISONS = {"MAM": [3, 4, 5], "JJA": [6, 7, 8], "SON": [9, 10, 11]}
PART_VALIDE_MIN = 0.5   # part minimale de pixels valides d'une classe par composite
ANNEES_MIN = 18         # annees minimum pour une tendance
COMPOSITES_MIN = 3      # composites valides minimum pour une moyenne de saison
SEUIL_FDR = 0.10


# ===========================================================================
# 1. Classes CORINE sur la grille MODIS
# ===========================================================================
if not os.path.exists(CLASSES_TIF):
    raise SystemExit("Introuvable : %s\n-> lancer d'abord 15_ndvi_corine.py" % CLASSES_TIF)
libelles = pd.read_csv(LIBELLES_CSV).set_index("code")["libelle"].to_dict()

composites = fichiers_apparies(MODIS, INDICE, ANNEES)
print("%d composites apparies (%s + 3 couches qualite)" % (len(composites), INDICE))

with rxr.open_rasterio(composites[0][2][INDICE]) as g:
    gabarit = g.squeeze(drop=True).load()
with rxr.open_rasterio(CLASSES_TIF) as g:
    cl = g.squeeze(drop=True).load()
# meme grille : appariement des centres de pixels, a un demi-pixel pres
demi = abs(float(gabarit.x[1] - gabarit.x[0])) / 2
try:
    cl = cl.sel(x=gabarit.x, y=gabarit.y, method="nearest", tolerance=demi)
except KeyError:
    raise SystemExit("La grille des classes CORINE ne recouvre pas la grille MODIS")
classes = np.nan_to_num(cl.values, nan=0).astype("int64")
classes[(classes < 1) | (classes > 8)] = 0
n_classe = np.bincount(classes.ravel(), minlength=9)
print("pixels par classe :", {libelles[k]: int(n_classe[k]) for k in range(1, 9)})

# ===========================================================================
# 2. NDVI moyen de chaque classe, composite par composite
# ===========================================================================
print("\n=== Lecture des composites ===")
lignes = []   # une ligne par composite, classe et mois d'acquisition
for i, (an, jour, f) in enumerate(composites):
    n = lire(f[INDICE])
    ok = masque_qualite(n, lire(f["pixel_reliability"]), lire(f["VI_Quality"]))
    annee, mois = date_des_pixels(an, jour, lire(f["composite_day_of_the_year"]))
    ok &= (mois > 0) & (classes > 0)
    nb_tot = np.bincount(classes[ok], minlength=9)
    garder = (n_classe > 0) & (nb_tot >= PART_VALIDE_MIN * n_classe)
    cle = annee * 100 + mois          # annee et mois d'acquisition, ex. 200307
    for am in np.unique(cle[ok]):
        a, m = divmod(int(am), 100)
        sel = ok & (cle == am)
        c = classes[sel]
        somme = np.bincount(c, weights=n[sel] * ECHELLE, minlength=9)
        nb = np.bincount(c, minlength=9)
        for k in range(1, 9):
            if garder[k] and nb[k] > 0:
                lignes.append({"COMPOSITE": i, "CLASSE": libelles[k], "ANNEE": a,
                               "MOIS": m, "SOMME": somme[k], "NB": nb[k]})
    if (i + 1) % 100 == 0:
        print("  %d/%d" % (i + 1, len(composites)))
comp = pd.DataFrame(lignes)

# ===========================================================================
# 3. Indicateurs annuels
# ===========================================================================
def moyennes(d):
    """Moyenne de tous les pixels valides par classe et annee, et nombre de
    composites qui y contribuent."""
    g = d.groupby(["CLASSE", "ANNEE"]).agg(SOMME=("SOMME", "sum"), NB=("NB", "sum"),
                                           N_COMP=("COMPOSITE", "nunique"))
    g["VALEUR"] = g["SOMME"] / g["NB"]
    return g


series = []
for m in MOIS:
    g = moyennes(comp[comp["MOIS"] == m])
    series.append(g["VALEUR"].reset_index().assign(INDICATEUR="M%02d" % m))
for s, mois_s in SAISONS.items():
    g = moyennes(comp[comp["MOIS"].isin(mois_s)])
    g = g[g["N_COMP"] >= COMPOSITES_MIN]
    series.append(g["VALEUR"].reset_index().assign(INDICATEUR=s))
# maximum : moyenne de chaque composite, rattachee au mois de la plupart de ses pixels
par_comp = comp.groupby(["COMPOSITE", "CLASSE"]).agg(SOMME=("SOMME", "sum"), NB=("NB", "sum"))
par_comp["VALEUR"] = par_comp["SOMME"] / par_comp["NB"]
principal = comp.loc[comp.groupby(["COMPOSITE", "CLASSE"])["NB"].idxmax(),
                     ["COMPOSITE", "CLASSE", "ANNEE", "MOIS"]].set_index(["COMPOSITE", "CLASSE"])
par_comp = par_comp.join(principal)
g = par_comp[par_comp["MOIS"].between(3, 11)].groupby(["CLASSE", "ANNEE"])["VALEUR"].max()
series.append(g.rename("VALEUR").reset_index().assign(INDICATEUR="MAX"))
series = pd.concat(series, ignore_index=True)
series.round(4).to_csv(os.path.join(RESULTATS, "%s_indicateurs_par_classe.csv" % SUFFIXE),
                       index=False, sep=";", decimal=",")

# ===========================================================================
# 4. Tendances
# ===========================================================================
print("\n=== Tendances (pente %s / decennie, 3PW, FDR %.0f %%) ===" % (INDICE, SEUIL_FDR * 100))
res = []
for (classe, ind), g in series.groupby(["CLASSE", "INDICATEUR"]):
    g = g.sort_values("ANNEE")
    if len(g) < ANNEES_MIN:
        continue
    pente, p_mk, p_pw, p_tf = tendance_serie(g["ANNEE"].values.astype(float),
                                             g["VALEUR"].values.astype(float))
    res.append({"CLASSE": classe, "INDICATEUR": ind, "N": len(g), "PENTE_DEC": pente * 10,
                "P_MK": p_mk, "P_3PW": float(p_3pw(p_pw, p_tf))})
res = pd.DataFrame(res)
res["P_FDR"] = fdr(res["P_3PW"].values)
res["SIGNIFICATIF"] = np.where(res["P_FDR"] < SEUIL_FDR, "oui", "non")
res = res.round({"PENTE_DEC": 4, "P_MK": 4, "P_3PW": 4, "P_FDR": 4})
res.to_csv(os.path.join(RESULTATS, "phenologie_%s_par_classe.csv" % SUFFIXE),
           index=False, sep=";", decimal=",")

ORDRE = ["M%02d" % m for m in MOIS] + list(SAISONS) + ["MAX"]
tab = res.pivot(index="CLASSE", columns="INDICATEUR", values="PENTE_DEC").reindex(columns=ORDRE)
sig = res.pivot(index="CLASSE", columns="INDICATEUR", values="SIGNIFICATIF").reindex(columns=ORDRE)
pval = res.pivot(index="CLASSE", columns="INDICATEUR", values="P_3PW").reindex(columns=ORDRE)
print("  %-22s" % "" + "".join("%8s" % c for c in ORDRE))
for classe in tab.index:
    print("  %-22s" % classe + "".join(
        "%8s" % ("" if np.isnan(v) else "%+.3f%s" % (v, "*" if s == "oui" else ""))
        for v, s in zip(tab.loc[classe], sig.loc[classe])))
print("  (* significatif apres correction ; case vide : moins de %d annees)" % ANNEES_MIN)

# ===========================================================================
# 5. Figure : pente par mois, puis par saison et maximum
# ===========================================================================
# vert-bleu = verdissement, brun = baisse (palette BrBG de ColorBrewer, lisible
# par les daltoniens)
VERT, BRUN, ENCRE, ENCRE_2, GRILLE, SURFACE = ("#018571", "#a6611a", "#0b0b0b",
                                               "#52514e", "#e5e4e0", "#fcfcfb")
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": GRILLE, "xtick.color": ENCRE_2, "ytick.color": ENCRE_2})
CLASSES_FIG = ["Cultures", "Prairies", "Foret feuillus"]
TEXTES = {
    "fr": dict(dossier="06_figures", noms={"Cultures": "Cultures", "Prairies": "Prairies",
                                           "Foret feuillus": "Forêt de feuillus"},
               mois=["mars", "avr.", "mai", "juin", "juil.", "août", "sept.", "oct.", "nov."],
               fin=["MAM", "JJA", "SON", "max"], y="Pente du %s (/ décennie)" % INDICE,
               titre="Tendance du %s par mois et par saison, 2001-2024" % INDICE,
               sous="Pente de Sen   |   * significative après correction (Mann-Kendall 3PW, FDR 10 %%)   |   "
                    "○ p < 0,05 pour ce test seul (exploratoire)   |   max : %s maximum de mars à novembre" % INDICE),
    "en": dict(dossier="06_figures_en", noms={"Cultures": "Cropland", "Prairies": "Grassland",
                                              "Foret feuillus": "Broadleaf forest"},
               mois=["Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"],
               fin=["MAM", "JJA", "SON", "max"], y="%s slope (/ decade)" % INDICE,
               titre="%s trend by month and season, 2001–2024" % INDICE,
               sous="Sen's slope   |   * significant after correction (Mann-Kendall 3PW, FDR 10 %%)   |   "
                    "○ p < 0.05 for this test alone (exploratory)   |   max: maximum %s from March to November" % INDICE),
}
CLASSES_FIG = [c for c in CLASSES_FIG if c in tab.index]
if not CLASSES_FIG:
    raise SystemExit("Aucune des classes de la figure n'a de tendance calculee")
lim = np.nanmax(np.abs(tab.loc[CLASSES_FIG].values)) * 1.25
for langue, t in TEXTES.items():
    os.makedirs(os.path.join(RACINE, t["dossier"]), exist_ok=True)
    fig, axes = plt.subplots(1, len(CLASSES_FIG), figsize=(14, 4), sharey=True)
    for ax, classe in zip(axes, CLASSES_FIG):
        if classe not in tab.index:
            continue
        x = np.r_[np.arange(9), np.arange(4) + 9.8]   # un espace avant saisons et max
        v = tab.loc[classe].values
        s = sig.loc[classe].values
        pv = pval.loc[classe].values
        ax.bar(x, np.nan_to_num(v), width=0.75,
               color=[VERT if vi > 0 else BRUN for vi in np.nan_to_num(v)])
        for xi, vi, si, pi in zip(x, v, s, pv):
            if si == "oui":
                marque, taille = "*", 13
            elif pi < 0.05:
                marque, taille = "○", 9
            else:
                continue
            ax.annotate(marque, (xi, vi), xytext=(0, 2 if vi > 0 else -2),
                        textcoords="offset points", ha="center",
                        va="bottom" if vi > 0 else "top",
                        fontsize=taille, fontweight="bold", color=ENCRE)
        ax.axhline(0, color=ENCRE_2, linewidth=0.8)
        ax.axvline(9.4, color=GRILLE, linewidth=1)
        ax.set_xticks(x)
        ax.set_xticklabels(t["mois"] + t["fin"], rotation=45, ha="right")
        ax.set_ylim(-lim, lim)
        ax.set_title(t["noms"][classe], fontweight="bold")
        ax.grid(axis="y", color=GRILLE, linewidth=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(t["y"])
    fig.suptitle(t["titre"], fontsize=13, fontweight="bold")
    fig.text(0.5, 0.9, t["sous"], ha="center", color=ENCRE_2, fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(RACINE, t["dossier"], "%s.%s" % (NOM_FIG, ext)),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

print("\n-> 05_resultats/%s_indicateurs_par_classe.csv" % SUFFIXE)
print("-> 05_resultats/phenologie_%s_par_classe.csv" % SUFFIXE)
print("-> 06_figures/ et 06_figures_en/%s.png / .pdf" % NOM_FIG)
