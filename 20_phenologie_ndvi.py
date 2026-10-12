# -*- coding: utf-8 -*-
"""
20_phenologie_ndvi.py
A lancer depuis la racine du projet, environnement geo :
    python 20_phenologie_ndvi.py

La baisse du NDVI moyen d'ete sur les cultures vient-elle d'une vegetation
moins vigoureuse, ou d'un calendrier qui avance (senescence et recolte plus
precoces) ? On regarde la tendance mois par mois :
  - si le printemps et le maximum ne baissent pas mais que juillet-aout
    baissent : plutot un decalage de calendrier ;
  - si tous les mois et le maximum baissent : plutot une perte de vigueur.

Pour chaque composite MODIS de 16 jours (2001-2024) :
  - meme masque qualite que 03_modis_ndvi.py (pixel_reliability 0 ou 1,
    valeur -3000 retiree, facteur 0,0001), meme date (debut + 8 jours) ;
  - NDVI moyen de chaque classe CORINE (8 classes, grille MODIS) ;
  - composite ignore pour une classe si moins de 50 % de ses pixels sont valides.
Indicateurs par classe et par annee :
  - NDVI moyen de chaque mois, mars a novembre (9 indicateurs) ;
  - NDVI moyen de chaque saison : MAM, JJA, SON (3 indicateurs) ;
  - NDVI maximum de mars a novembre (1 indicateur).
Tendance de chaque serie : pente de Sen, Mann-Kendall 3PW, puis correction
FDR 10 % sur l'ensemble des tests (outils_tendance.py). Au moins 18 annees.

Le JJA calcule ici (moyenne spatiale de chaque composite, puis moyenne de la
saison) peut differer legerement de celui de 16_analyse.py (moyenne de la
saison par pixel, puis moyenne spatiale) : l'ordre des moyennes change.

Entrees :
  01_raw/MODIS/NDVI/*.tif, 01_raw/MODIS/pixel_reliability/*.tif
  01_raw/corine_8classes_grille_modis.tif   (15_ndvi_corine.py)
  04_stats/corine_8classes_libelles.csv
Sorties :
  05_resultats/ndvi_indicateurs_par_classe.csv   series annuelles
  05_resultats/phenologie_ndvi_par_classe.csv    tendances
  06_figures/fig5_tendance_ndvi_mensuelle.png/.pdf (et 06_figures_en)
"""

import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rioxarray as rxr

from outils_tendance import fdr, p_3pw, tendance_serie

RACINE = os.getcwd()   # racine du projet
RAW_NDVI = os.path.join(RACINE, "01_raw", "MODIS", "NDVI")
RAW_REL = os.path.join(RACINE, "01_raw", "MODIS", "pixel_reliability")
CLASSES_TIF = os.path.join(RACINE, "01_raw", "corine_8classes_grille_modis.tif")
LIBELLES_CSV = os.path.join(RACINE, "04_stats", "corine_8classes_libelles.csv")
RESULTATS = os.path.join(RACINE, "05_resultats")
os.makedirs(RESULTATS, exist_ok=True)

ANNEES = (2001, 2024)
MOIS = [3, 4, 5, 6, 7, 8, 9, 10, 11]
SAISONS = {"MAM": [3, 4, 5], "JJA": [6, 7, 8], "SON": [9, 10, 11]}
PART_VALIDE_MIN = 0.5   # part minimale de pixels valides d'une classe par composite
ANNEES_MIN = 18         # annees minimum pour une tendance
SEUIL_FDR = 0.10


def date_du_fichier(chemin):
    """Date mediane du composite de 16 jours (comme 03_modis_ndvi.py)."""
    m = re.search(r"doy(\d{4})(\d{3})", os.path.basename(chemin))
    an, jour = int(m.group(1)), int(m.group(2))
    return pd.Timestamp(f"{an}-01-01") + pd.Timedelta(days=jour - 1 + 8)


# ===========================================================================
# 1. Classes CORINE sur la grille MODIS
# ===========================================================================
if not os.path.exists(CLASSES_TIF):
    raise SystemExit("Introuvable : %s\n-> exporter corine_grille_ndvi depuis ArcGIS "
                     "(voir 15_ndvi_corine.py)" % CLASSES_TIF)
libelles = pd.read_csv(LIBELLES_CSV).set_index("code")["libelle"].to_dict()

ndvi_f = {date_du_fichier(f): f for f in glob.glob(os.path.join(RAW_NDVI, "*.tif"))}
rel_f = {date_du_fichier(f): f for f in glob.glob(os.path.join(RAW_REL, "*.tif"))}
dates = sorted(d for d in set(ndvi_f) & set(rel_f) if ANNEES[0] <= d.year <= ANNEES[1])
print("%d composites apparies, %s -> %s" % (len(dates), dates[0].date(), dates[-1].date()))

gabarit = rxr.open_rasterio(ndvi_f[dates[0]]).squeeze(drop=True)
cl = rxr.open_rasterio(CLASSES_TIF).squeeze(drop=True)
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
lignes = []
for i, d in enumerate(dates):
    n = rxr.open_rasterio(ndvi_f[d], masked=False).squeeze(drop=True).values
    r = rxr.open_rasterio(rel_f[d], masked=False).squeeze(drop=True).values
    valide = ((r == 0) | (r == 1)) & (n != -3000) & (classes > 0)
    c = classes[valide]
    somme = np.bincount(c, weights=n[valide] * 0.0001, minlength=9)
    nb = np.bincount(c, minlength=9)
    for k in range(1, 9):
        if n_classe[k] > 0 and nb[k] >= PART_VALIDE_MIN * n_classe[k]:
            lignes.append({"DATE": d, "CLASSE": libelles[k], "NDVI": somme[k] / nb[k]})
    if (i + 1) % 100 == 0:
        print("  %d/%d" % (i + 1, len(dates)))
comp = pd.DataFrame(lignes)
comp["ANNEE"] = comp["DATE"].dt.year
comp["MOIS"] = comp["DATE"].dt.month

# ===========================================================================
# 3. Indicateurs annuels
# ===========================================================================
series = []
for m in MOIS:
    g = comp[comp["MOIS"] == m].groupby(["CLASSE", "ANNEE"])["NDVI"].mean()
    series.append(g.rename("VALEUR").reset_index().assign(INDICATEUR="M%02d" % m))
for s, mois in SAISONS.items():
    g = comp[comp["MOIS"].isin(mois)].groupby(["CLASSE", "ANNEE"])["NDVI"].mean()
    series.append(g.rename("VALEUR").reset_index().assign(INDICATEUR=s))
g = comp[comp["MOIS"].between(3, 11)].groupby(["CLASSE", "ANNEE"])["NDVI"].max()
series.append(g.rename("VALEUR").reset_index().assign(INDICATEUR="MAX"))
series = pd.concat(series, ignore_index=True)
series.round(4).to_csv(os.path.join(RESULTATS, "ndvi_indicateurs_par_classe.csv"),
                       index=False, sep=";", decimal=",")

# ===========================================================================
# 4. Tendances
# ===========================================================================
print("\n=== Tendances (pente NDVI / decennie, 3PW, FDR %.0f %%) ===" % (SEUIL_FDR * 100))
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
res.to_csv(os.path.join(RESULTATS, "phenologie_ndvi_par_classe.csv"),
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
               fin=["MAM", "JJA", "SON", "max"], y="Pente du NDVI (/ décennie)",
               titre="Tendance du NDVI par mois et par saison, 2001-2024",
               sous="Pente de Sen   |   * significative après correction (Mann-Kendall 3PW, FDR 10 %)   |   "
                    "○ p < 0,05 pour ce test seul (exploratoire)   |   max : NDVI maximum de mars à novembre"),
    "en": dict(dossier="06_figures_en", noms={"Cultures": "Cropland", "Prairies": "Grassland",
                                              "Foret feuillus": "Broadleaf forest"},
               mois=["Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"],
               fin=["MAM", "JJA", "SON", "max"], y="NDVI slope (/ decade)",
               titre="NDVI trend by month and season, 2001–2024",
               sous="Sen's slope   |   * significant after correction (Mann-Kendall 3PW, FDR 10 %)   |   "
                    "○ p < 0.05 for this test alone (exploratory)   |   max: maximum NDVI from March to November"),
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
        fig.savefig(os.path.join(RACINE, t["dossier"], "fig5_tendance_ndvi_mensuelle.%s" % ext),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

print("\n-> 05_resultats/ndvi_indicateurs_par_classe.csv")
print("-> 05_resultats/phenologie_ndvi_par_classe.csv")
print("-> 06_figures/ et 06_figures_en/fig5_tendance_ndvi_mensuelle.png / .pdf")
