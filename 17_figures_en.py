# -*- coding: utf-8 -*-
"""
17_figures_en.py
A lancer depuis la racine du projet, environnement geo :
    python 17_figures_en.py
Necessite matplotlib :  conda install -c conda-forge matplotlib

Version anglaise de 17_figures.py. Produit dans 06_figures_en\\ quatre
figures, en PNG 300 dpi et en PDF vectoriel :

  fig1_seasonal_temperature : series 1961-2024 par saison, avec tendance
                              et les deux normales OMM
  fig2_normals_change       : ecart 1991-2020 moins 1961-1990, par unite
  fig3_ndvi_by_landcover    : NDVI estival 2001-2024 par classe CORINE
  fig4_correlations         : correlations NDVI / climat

Choix graphiques : une seule echelle par panneau (jamais deux axes y),
palette validee pour le daltonisme, valeurs ecrites en toutes lettres la ou
elles portent le message plutot que sur chaque point.

Les pentes affichees (figures 1 et 3) sont des moindres carres sur la moyenne
nationale : elles different des pentes de Sen par pixel des cartes de tendance.
"""

import os

import matplotlib
matplotlib.use("Agg")                      # pas d'affichage : ecriture directe
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

RACINE = os.getcwd()   # racine du projet
RESULTATS = os.path.join(RACINE, "05_resultats")
SORTIE = os.path.join(RACINE, "06_figures_en")
os.makedirs(SORTIE, exist_ok=True)

# --- Palette ----------------------------------------------------------------
BLEU = "#2a78d6"
ORANGE = "#eb6834"
ROUGE = "#e34948"
SURFACE = "#fcfcfb"
ENCRE = "#0b0b0b"
ENCRE_2 = "#52514e"
GRILLE = "#e5e4e0"
GRIS_FOND = "#d8d7d2"
DIVERGENTE = LinearSegmentedColormap.from_list("bleu_rouge", [BLEU, "#f0efec", ROUGE])

SAISONS = ["DJF", "MAM", "JJA", "SON"]
NOM_SAISON = {"DJF": "Winter (DJF)", "MAM": "Spring (MAM)",
              "JJA": "Summer (JJA)", "SON": "Autumn (SON)"}

# Les libelles stockes dans les CSV sont sans accents (contrainte des champs
# texte ecrits sous ArcGIS) : ils sont traduits ici, a l'affichage seulement
NOM_CLASSE = {
    "Artificialise": "Artificial surfaces",
    "Cultures": "Cropland",
    "Prairies": "Grassland",
    "Foret feuillus": "Broadleaf forest",
    "Foret resineux": "Coniferous forest",
    "Foret mixte": "Mixed forest",
    "Milieux semi-naturels": "Semi-natural areas",
    "Zones humides et eau": "Wetlands and water",
}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": ENCRE, "axes.labelcolor": ENCRE_2,
    "xtick.color": ENCRE_2, "ytick.color": ENCRE_2,
    "axes.edgecolor": GRILLE, "grid.color": GRILLE,
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False,
})


def lire(nom):
    return pd.read_csv(os.path.join(RESULTATS, nom + ".csv"), sep=";", decimal=",")


def enregistrer(fig, nom):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(SORTIE, "%s.%s" % (nom, ext)),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  -> %s.png / .pdf" % nom)


def tendance(annees, valeurs):
    """Pente par decennie, par moindres carres."""
    a, b = np.polyfit(annees, valeurs, 1)
    return a * 10, a * np.asarray(annees) + b


def axes_propres(ax):
    ax.grid(axis="y", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


series = lire("series_nationales")
print("=== Figures ===")

# ===========================================================================
# Figure 1 : temperature saisonniere, un panneau par saison
# ===========================================================================
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax, s in zip(axes.ravel(), SAISONS):
    d = series[(series["INDICATEUR"] == "t2m") & (series["SAISON"] == s)]
    d = d.sort_values("ANNEE")
    an, va = d["ANNEE"].values, d["VALEUR"].values

    ax.plot(an, va, color=BLEU, linewidth=1.6, label="Seasonal mean")
    pente, droite = tendance(an, va)
    ax.plot(an, droite, color=ENCRE_2, linewidth=1.6, linestyle="--",
            label="Linear trend")

    # Les deux normales OMM, en segments horizontaux
    for a1, a2, style in [(1961, 1990, ":"), (1991, 2020, "-")]:
        m = va[(an >= a1) & (an <= a2)].mean()
        ax.plot([a1, a2], [m, m], color=ORANGE, linewidth=2.4, linestyle=style,
                solid_capstyle="butt",
                label="Normal %d–%d" % (a1, a2) if s == "DJF" else None)

    n1 = va[(an >= 1961) & (an <= 1990)].mean()
    n2 = va[(an >= 1991) & (an <= 2020)].mean()
    ax.set_title("%s   %+.2f °C between normals" % (NOM_SAISON[s], n2 - n1))
    ax.annotate("%+.2f °C / decade" % pente, xy=(0.03, 0.92),
                xycoords="axes fraction", color=ENCRE_2, fontsize=8.5)
    ax.set_ylabel("Temperature (°C)")
    axes_propres(ax)

handles, labels = axes.ravel()[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, -0.03))
fig.suptitle("Mean seasonal temperature in Belgium, 1961–2024",
             fontsize=13, fontweight="bold", y=0.98)
fig.text(0.5, 0.925, "ERA5-Land, area-weighted national mean",
         ha="center", color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0.02, 1, 0.92])
enregistrer(fig, "fig1_seasonal_temperature")

# ===========================================================================
# Figure 2 : ecart entre normales, un panneau par unite
# ===========================================================================
UNITES_FIG = [
    ("degC", "Temperature (°C)",
     [("t2m", s) for s in SAISONS] + [("txx", "JJA")]),
    ("jours", "Number of days",
     [("su", "JJA"), ("su30", "JJA"), ("tr", "JJA"), ("fd", "DJF"), ("fd", "MAM")]),
    ("mm", "Precipitation (mm)",
     [("prcptot", s) for s in SAISONS]),
    ("m3/m3", "Soil moisture (m³/m³)",
     [("swvl1", s) for s in SAISONS]),
]
NOMS = {"t2m": "Mean temperature", "txx": "Hottest day (TXx)",
        "su": "Days > 25 °C", "su30": "Days ≥ 30 °C",
        "tr": "Tropical nights", "fd": "Frost days",
        "prcptot": "Wet-day total (E-OBS)", "swvl1": "Soil moisture 0–7 cm"}

fig, axes = plt.subplots(1, 4, figsize=(14, 4.6))
for ax, (unite, titre, liste) in zip(axes, UNITES_FIG):
    etiquettes, ecarts = [], []
    for ind, sai in liste:
        d = series[(series["INDICATEUR"] == ind) & (series["SAISON"] == sai)]
        v, a = d["VALEUR"].values, d["ANNEE"].values
        n1 = v[(a >= 1961) & (a <= 1990)].mean()
        n2 = v[(a >= 1991) & (a <= 2020)].mean()
        etiquettes.append("%s\n%s" % (NOMS.get(ind, ind), sai))
        ecarts.append(n2 - n1)

    y = np.arange(len(etiquettes))
    couleurs = [ROUGE if e > 0 else BLEU for e in ecarts]
    ax.barh(y, ecarts, color=couleurs, edgecolor=SURFACE, linewidth=1.2, height=0.62)
    ax.axvline(0, color=ENCRE_2, linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(etiquettes, fontsize=8)
    ax.invert_yaxis()
    ax.set_title(titre)
    ax.set_xlabel("Change 1991–2020 vs 1961–1990")
    ax.grid(axis="x", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    # valeur ecrite au bout de chaque barre : la lecture ne depend pas de la couleur
    marge = max(abs(min(ecarts)), abs(max(ecarts))) * 0.04
    for yi, e in zip(y, ecarts):
        fmt = "%+.3f" if unite == "m3/m3" else "%+.1f"
        ax.annotate(fmt % e, xy=(e + (marge if e >= 0 else -marge), yi),
                    va="center", ha="left" if e >= 0 else "right",
                    fontsize=8, color=ENCRE)
    lim = max(abs(min(ecarts)), abs(max(ecarts))) * 1.45
    ax.set_xlim(-lim, lim)

fig.suptitle("Change between the two climate normals (1961–1990 → 1991–2020)",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.905, "Red: increase   |   Blue: decrease", ha="center",
         color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.88])
enregistrer(fig, "fig2_normals_change")

# ===========================================================================
# Figure 3 : NDVI estival par classe d'occupation du sol
# ===========================================================================
ndvi = lire("ndvi_par_classe")
ete = ndvi[ndvi["SAISON"] == "JJA"].sort_values(["CLASSE", "ANNEE"])
classes = sorted(ete["CLASSE"].unique())

fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharex=True, sharey=True)
for ax, cl in zip(axes.ravel(), classes):
    # toutes les classes en fond gris : chaque panneau se lit dans le contexte
    for autre in classes:
        g = ete[ete["CLASSE"] == autre]
        ax.plot(g["ANNEE"], g["VALEUR"], color=GRIS_FOND, linewidth=0.9, zorder=1)

    d = ete[ete["CLASSE"] == cl]
    an, va = d["ANNEE"].values, d["VALEUR"].values
    ax.plot(an, va, color=BLEU, linewidth=1.9, zorder=3)
    pente, droite = tendance(an, va)
    ax.plot(an, droite, color=ENCRE_2, linewidth=1.4, linestyle="--", zorder=4)

    ax.set_title(NOM_CLASSE.get(cl, cl), fontsize=9.5)
    ax.annotate("%+.3f / decade" % pente, xy=(0.04, 0.06),
                xycoords="axes fraction", fontsize=8, color=ENCRE_2)
    axes_propres(ax)

for ax in axes[:, 0]:
    ax.set_ylabel("NDVI")
fig.suptitle("Summer NDVI by land-cover class, 2001–2024",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.925,
         "MODIS MOD13Q1, CORINE Land Cover 2018 classes — other classes in grey",
         ha="center", color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.91])
enregistrer(fig, "fig3_ndvi_by_landcover")

# ===========================================================================
# Figure 4 : correlations NDVI / climat
# ===========================================================================
correl = lire("correlation_ndvi_climat")
NOM_CLIM = {"t2m": "Temperature", "tp": "Precipitation", "swvl1": "Soil moisture"}

fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=True)
for ax, s in zip(axes, SAISONS):
    d = correl[correl["SAISON"] == s]
    tab = d.pivot(index="CLASSE", columns="CLIMAT", values="R")
    tab = tab.reindex(columns=[c for c in ("t2m", "tp", "swvl1") if c in tab.columns])
    sig = d.pivot(index="CLASSE", columns="CLIMAT", values="SIGNIFICATIF")
    sig = sig.reindex(columns=tab.columns)

    ax.imshow(tab.values, cmap=DIVERGENTE, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(tab.columns)))
    ax.set_xticklabels([NOM_CLIM.get(c, c) for c in tab.columns],
                       rotation=30, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(tab.index)))
    ax.set_yticklabels([NOM_CLASSE.get(c, c) for c in tab.index], fontsize=8.5)
    ax.set_title(NOM_SAISON[s])
    ax.set_xticks(np.arange(-0.5, len(tab.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(tab.index), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(axis="y", length=0)

    # la valeur est ecrite dans chaque case : la couleur appuie, elle ne porte pas seule
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            r = tab.values[i, j]
            if np.isnan(r):
                continue
            etoile = "*" if sig.values[i, j] == "oui" else ""
            ax.text(j, i, "%+.2f%s" % (r, etoile), ha="center", va="center",
                    fontsize=8, fontweight="bold" if etoile else "normal",
                    color="#ffffff" if abs(r) > 0.55 else ENCRE)

fig.suptitle("Correlation between NDVI and climate by land-cover class, 2001–2024",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.9,
         "Pearson coefficient, n = 24 years   |   * significant at 5 % (|r| > 0.404)   |"
         "   Red: positive, Blue: negative",
         ha="center", color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.86])
enregistrer(fig, "fig4_correlations")

print("\nTermine. Figures dans : %s" % SORTIE)