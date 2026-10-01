# -*- coding: utf-8 -*-
"""
16_analyse.py
A lancer depuis la racine du projet, environnement geo :
    python 16_analyse.py

Exploite les CSV de statistiques zonales produits sous ArcGIS (14 et 15) pour construire
les chiffres du rapport. Ne depend plus d'ArcGIS ni de la licence.

Produit dans 05_resultats\\ :
  - normales_par_province.csv  : 1961-1990 vs 1991-2020 vs 2015-2024, par province
  - series_nationales.csv      : moyenne belge annuelle, par indicateur et saison
  - ndvi_par_classe.csv        : NDVI national par classe d'occupation du sol
  - correlation_ndvi_climat.csv: correlation NDVI / temperature, precipitations, humidite du sol

Details de format geres ici :
  - CSV francais : separateur ';', virgule decimale
  - dates en jj-mm-aa sur deux chiffres : '61' = 1961, pas 2061
  - moyennes nationales ponderees par COUNT (nombre de pixels), donc par surface
"""

import os
import sys

import pandas as pd

RACINE = os.getcwd()   # racine du projet
STATS = os.path.join(RACINE, "04_stats")
SORTIE = os.path.join(RACINE, "05_resultats")
os.makedirs(SORTIE, exist_ok=True)

# Periodes de reference OMM + decennie illustrative
PERIODES = {
    "1961-1990": (1961, 1990),
    "1991-2020": (1991, 2020),
    "2015-2024": (2015, 2024),
}

# CSV de statistiques zonales : nom de fichier -> (source, saison)
CLIMAT = {
    "zs_era5_DJF": ("ERA5", "DJF"), "zs_era5_MAM": ("ERA5", "MAM"),
    "zs_era5_JJA": ("ERA5", "JJA"), "zs_era5_SON": ("ERA5", "SON"),
    "zs_txx_JJA": ("txx", "JJA"), "zs_su_JJA": ("su", "JJA"),
    "zs_tr_JJA": ("tr", "JJA"), "zs_su30_JJA": ("su30", "JJA"),
    "zs_fd_DJF": ("fd", "DJF"), "zs_fd_MAM": ("fd", "MAM"),
    "zs_prcptot_DJF": ("prcptot", "DJF"), "zs_prcptot_MAM": ("prcptot", "MAM"),
    "zs_prcptot_JJA": ("prcptot", "JJA"), "zs_prcptot_SON": ("prcptot", "SON"),
}

UNITES = {
    "t2m": "degC", "tp": "mm", "swvl1": "m3/m3", "swvl2": "m3/m3",
    "txx": "degC", "su": "jours", "tr": "jours", "su30": "jours",
    "fd": "jours", "prcptot": "mm", "ndvi": "-",
}


def annee_depuis_date(texte):
    """'15-07-61 00:00:00' -> 1961 ; '1961-07-15' -> 1961.

    Les annees sur deux chiffres sont ambigues : ici les donnees couvrent
    1961-2024, donc un '61' est 1961 et un '01' est 2001. Seuil a 30.
    """
    date = str(texte).strip().split(" ")[0]
    morceaux = date.replace("/", "-").split("-")
    if len(morceaux[0]) == 4:              # aaaa-mm-jj
        return int(morceaux[0])
    aa = int(morceaux[-1])                 # jj-mm-aa
    if len(morceaux[-1]) == 4:
        return aa
    return 1900 + aa if aa > 30 else 2000 + aa


def charger(nom):
    """Lit un CSV de statistiques zonales et le met en forme."""
    chemin = os.path.join(STATS, nom + ".csv")
    if not os.path.exists(chemin):
        print("  MANQUANT : %s" % os.path.basename(chemin))
        return None
    d = pd.read_csv(chemin, sep=";", decimal=",")
    d.columns = [c.strip() for c in d.columns]
    d["ANNEE"] = d["StdTime"].map(annee_depuis_date)
    for c in ("NUTS_ID", "CLASSE", "Variable"):
        if c in d.columns and d[c].dtype == object:
            d[c] = d[c].astype(str).str.strip()
    return d


def moyenne_nationale(d, cles):
    """Moyenne ponderee par le nombre de pixels : une moyenne de surface,
    pas une moyenne des provinces (qui donnerait le meme poids a Bruxelles
    qu'au Hainaut)."""
    d = d.copy()
    d["_num"] = d["MEAN"] * d["COUNT"]
    g = d.groupby(cles).agg(_num=("_num", "sum"), _den=("COUNT", "sum"))
    g["VALEUR"] = g["_num"] / g["_den"]
    return g["VALEUR"].reset_index()


# ===========================================================================
# 1. Chargement
# ===========================================================================
print("=== Chargement des CSV ===")
climat = []
for fichier, (source, saison) in CLIMAT.items():
    d = charger(fichier)
    if d is None:
        continue
    d["SAISON"] = saison
    # ERA5 contient 4 variables, les indices E-OBS une seule
    d["INDICATEUR"] = d["Variable"].str.lower()
    climat.append(d[["NUTS_ID", "INDICATEUR", "SAISON", "ANNEE", "MEAN", "COUNT"]])
    print("  %-16s %5d lignes  %s" % (fichier, len(d), sorted(d["INDICATEUR"].unique())))

if not climat:
    sys.exit("Aucun CSV lu : verifier le dossier 04_stats")
climat = pd.concat(climat, ignore_index=True)

ndvi = []
for s in ["DJF", "MAM", "JJA", "SON"]:
    d = charger("zs_ndvi_%s" % s)
    if d is None:
        continue
    d["SAISON"] = s
    ndvi.append(d[["NUTS_ID", "CLASSE", "CLC_CODE", "SAISON", "ANNEE", "MEAN", "COUNT"]])
    print("  zs_ndvi_%s      %5d lignes" % (s, len(d)))
ndvi = pd.concat(ndvi, ignore_index=True) if ndvi else None

print("\n  annees climat : %d-%d   annees NDVI : %s" % (
    climat["ANNEE"].min(), climat["ANNEE"].max(),
    "%d-%d" % (ndvi["ANNEE"].min(), ndvi["ANNEE"].max()) if ndvi is not None else "-"))

# ===========================================================================
# 2. Normales par province
# ===========================================================================
print("\n=== Normales par province ===")
lignes = []
for (prov, ind, sai), g in climat.groupby(["NUTS_ID", "INDICATEUR", "SAISON"]):
    ligne = {"NUTS_ID": prov, "INDICATEUR": ind, "SAISON": sai,
             "UNITE": UNITES.get(ind, "")}
    for nom, (a1, a2) in PERIODES.items():
        sel = g[(g["ANNEE"] >= a1) & (g["ANNEE"] <= a2)]
        # une normale demande la periode complete
        ligne[nom] = sel["MEAN"].mean() if len(sel) == (a2 - a1 + 1) else None
    if ligne["1961-1990"] is not None and ligne["1991-2020"] is not None:
        ligne["ECART"] = ligne["1991-2020"] - ligne["1961-1990"]
    lignes.append(ligne)

normales = pd.DataFrame(lignes).round(3)
normales.to_csv(os.path.join(SORTIE, "normales_par_province.csv"),
                index=False, sep=";", decimal=",")
print("  -> normales_par_province.csv  (%d lignes)" % len(normales))

# Synthese nationale, pour le texte du rapport
print("\n  Ecart 1991-2020 moins 1961-1990, moyenne simple des provinces")
print("  (non ponderee, pour controle ; series ponderees : series_nationales.csv) :")
synth = normales.dropna(subset=["ECART"]).groupby(["INDICATEUR", "SAISON", "UNITE"])["ECART"].mean()
for (ind, sai, u), v in synth.items():
    print("    %-9s %-4s %+7.2f %s" % (ind, sai, v, u))

# ===========================================================================
# 3. Series nationales annuelles
# ===========================================================================
print("\n=== Series nationales ===")
series = moyenne_nationale(climat, ["INDICATEUR", "SAISON", "ANNEE"])
series["UNITE"] = series["INDICATEUR"].map(UNITES)
series = series.sort_values(["INDICATEUR", "SAISON", "ANNEE"]).round(3)
series.to_csv(os.path.join(SORTIE, "series_nationales.csv"),
              index=False, sep=";", decimal=",")
print("  -> series_nationales.csv  (%d lignes)" % len(series))

t = series[(series["INDICATEUR"] == "t2m") & (series["SAISON"] == "JJA")]
print("  controle t2m JJA : 1961 = %.2f   2003 = %.2f   2024 = %.2f degC" % (
    t[t["ANNEE"] == 1961]["VALEUR"].iloc[0],
    t[t["ANNEE"] == 2003]["VALEUR"].iloc[0],
    t[t["ANNEE"] == 2024]["VALEUR"].iloc[0]))

# ===========================================================================
# 4. NDVI par classe d'occupation du sol
# ===========================================================================
if ndvi is not None:
    print("\n=== NDVI par classe ===")
    ndvi_nat = moyenne_nationale(ndvi, ["CLASSE", "SAISON", "ANNEE"])
    ndvi_nat = ndvi_nat.sort_values(["CLASSE", "SAISON", "ANNEE"]).round(4)
    ndvi_nat.to_csv(os.path.join(SORTIE, "ndvi_par_classe.csv"),
                    index=False, sep=";", decimal=",")
    print("  -> ndvi_par_classe.csv  (%d lignes)" % len(ndvi_nat))

    # =======================================================================
    # 5. Correlation NDVI / climat, 2001-2024
    # =======================================================================
    print("\n=== Correlation NDVI / climat (2001-2024, n=24) ===")
    print("  |r| > 0.404 = significatif a 5 %%\n")

    clim_nat = moyenne_nationale(climat, ["INDICATEUR", "SAISON", "ANNEE"])
    resultats = []
    for (classe, saison), g in ndvi_nat.groupby(["CLASSE", "SAISON"]):
        g = g.set_index("ANNEE")["VALEUR"]
        for ind in ("t2m", "tp", "swvl1"):
            c = clim_nat[(clim_nat["INDICATEUR"] == ind) & (clim_nat["SAISON"] == saison)]
            c = c.set_index("ANNEE")["VALEUR"]
            commun = g.index.intersection(c.index)
            if len(commun) < 10:
                continue
            r = g.loc[commun].corr(c.loc[commun])
            resultats.append({"CLASSE": classe, "SAISON": saison, "CLIMAT": ind,
                              "N": len(commun), "R": round(r, 3),
                              "SIGNIFICATIF": "oui" if abs(r) > 0.404 else "non"})

    correl = pd.DataFrame(resultats)
    correl.to_csv(os.path.join(SORTIE, "correlation_ndvi_climat.csv"),
                  index=False, sep=";", decimal=",")
    print("  -> correlation_ndvi_climat.csv  (%d lignes)" % len(correl))

    print("\n  Correlations significatives en ete (JJA) :")
    ete = correl[(correl["SAISON"] == "JJA") & (correl["SIGNIFICATIF"] == "oui")]
    if len(ete) == 0:
        print("    aucune")
    for _, r in ete.sort_values("R").iterrows():
        print("    %-22s %-6s r = %+.3f" % (r["CLASSE"], r["CLIMAT"], r["R"]))

print("\nTermine. Resultats dans : %s" % SORTIE)