# -*- coding: utf-8 -*-
"""
17_figures.py
A lancer depuis la racine du projet, environnement geo :
    python scripts/17_figures.py
Necessite matplotlib :  conda install -c conda-forge matplotlib

Produit dans 06_figures\\ quatre figures, en PNG 300 dpi (rapport Word,
StoryMap) et en PDF vectoriel (impression) :

  fig1_temperature_saisonniere : series 1961-2024 par saison, avec tendance
                                 et les deux normales OMM
  fig2_ecarts_normales         : ecart 1991-2020 moins 1961-1990, par unite
  fig3_ndvi_par_classe         : NDVI estival 2001-2024 par classe CORINE
  fig4_correlations            : correlations NDVI / climat

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
SORTIE = os.path.join(RACINE, "06_figures")
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
NOM_SAISON = {"DJF": "Hiver (DJF)", "MAM": "Printemps (MAM)",
              "JJA": "Été (JJA)", "SON": "Automne (SON)"}

# Les libelles stockes dans les CSV sont sans accents (contrainte des champs
# texte ecrits sous ArcGIS) : ils sont reaccentues ici, a l'affichage seulement
NOM_CLASSE = {
    "Artificialise": "Artificialisé",
    "Cultures": "Cultures",
    "Prairies": "Prairies",
    "Foret feuillus": "Forêt de feuillus",
    "Foret resineux": "Forêt de résineux",
    "Foret mixte": "Forêt mixte",
    "Milieux semi-naturels": "Milieux semi-naturels",
    "Zones humides et eau": "Zones humides et eau",
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

    ax.plot(an, va, color=BLEU, linewidth=1.6, label="Moyenne saisonnière")
    pente, droite = tendance(an, va)
    ax.plot(an, droite, color=ENCRE_2, linewidth=1.6, linestyle="--",
            label="Tendance linéaire")

    # Les deux normales OMM, en segments horizontaux
    for a1, a2, style in [(1961, 1990, ":"), (1991, 2020, "-")]:
        m = va[(an >= a1) & (an <= a2)].mean()
        ax.plot([a1, a2], [m, m], color=ORANGE, linewidth=2.4, linestyle=style,
                solid_capstyle="butt",
                label="Normale %d-%d" % (a1, a2) if s == "DJF" else None)

    n1 = va[(an >= 1961) & (an <= 1990)].mean()
    n2 = va[(an >= 1991) & (an <= 2020)].mean()
    ax.set_title("%s   %+.2f °C entre normales" % (NOM_SAISON[s], n2 - n1))
    ax.annotate("%+.2f °C / décennie" % pente, xy=(0.03, 0.92),
                xycoords="axes fraction", color=ENCRE_2, fontsize=8.5)
    ax.set_ylabel("Température (°C)")
    axes_propres(ax)

handles, labels = axes.ravel()[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, -0.03))
fig.suptitle("Température moyenne saisonnière en Belgique, 1961-2024",
             fontsize=13, fontweight="bold", y=0.98)
fig.text(0.5, 0.925, "ERA5-Land, moyenne nationale pondérée par la surface",
         ha="center", color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0.02, 1, 0.92])
enregistrer(fig, "fig1_temperature_saisonniere")

# ===========================================================================
# Figure 2 : ecart entre normales, un panneau par unite
# ===========================================================================
UNITES_FIG = [
    ("degC", "Température (°C)",
     [("t2m", s) for s in SAISONS] + [("txx", "JJA")]),
    ("jours", "Nombre de jours",
     [("su", "JJA"), ("su30", "JJA"), ("tr", "JJA"), ("fd", "DJF"), ("fd", "MAM")]),
    ("mm", "Précipitations (mm)",
     [("prcptot", s) for s in SAISONS]),
    ("m3/m3", "Humidité du sol (m³/m³)",
     [("swvl1", s) for s in SAISONS]),
]
NOMS = {"t2m": "Temp. moyenne", "txx": "Jour le plus chaud (TXx)",
        "su": "Jours > 25 °C", "su30": "Jours ≥ 30 °C",
        "tr": "Nuits tropicales", "fd": "Jours de gel",
        "prcptot": "Cumul (E-OBS)", "swvl1": "Humidité 0-7 cm"}

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
    ax.set_xlabel("Écart 1991-2020 / 1961-1990")
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

fig.suptitle("Évolution entre les deux normales climatiques (1961-1990 → 1991-2020)",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.905, "Rouge : hausse   |   Bleu : baisse", ha="center",
         color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.88])
enregistrer(fig, "fig2_ecarts_normales")

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
    ax.annotate("%+.3f / décennie" % pente, xy=(0.04, 0.06),
                xycoords="axes fraction", fontsize=8, color=ENCRE_2)
    axes_propres(ax)

for ax in axes[:, 0]:
    ax.set_ylabel("NDVI")
fig.suptitle("NDVI estival par classe d'occupation du sol, 2001-2024",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.925,
         "MODIS MOD13Q1, classes CORINE Land Cover 2018 — en gris, les autres classes",
         ha="center", color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.91])
enregistrer(fig, "fig3_ndvi_par_classe")

# ===========================================================================
# Figure 4 : correlations NDVI / climat
# ===========================================================================
correl = lire("correlation_ndvi_climat")
NOM_CLIM = {"t2m": "Température", "tp": "Précipitations", "swvl1": "Humidité du sol"}

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

fig.suptitle("Corrélation entre NDVI et climat par classe d'occupation du sol, 2001-2024",
             fontsize=13, fontweight="bold")
fig.text(0.5, 0.9,
         "Coefficient de Pearson, n = 24 ans   |   * significatif à 5 % (|r| > 0,404)   |"
         "   Rouge : positif, Bleu : négatif",
         ha="center", color=ENCRE_2, fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.86])
enregistrer(fig, "fig4_correlations")

print("\nTermine. Figures dans : %s" % SORTIE)