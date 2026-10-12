"""Beschrijvende statistiek en groepsvergelijkingen ("Tabel 1")."""

import numpy as np
import pandas as pd
from scipy import stats

from . import config


def _is_normaal(waarden):
    """Shapiro-Wilk; bij n < 3 of n > 5000 wordt normaliteit niet getoetst."""
    waarden = waarden.dropna()
    if len(waarden) < 3 or len(waarden) > 5000:
        return len(waarden) > 5000  # bij grote n: centrale limietstelling
    return stats.shapiro(waarden).pvalue >= config.ALPHA_NORMALITEIT


def _fmt_p(p):
    if p is None or np.isnan(p):
        return ""
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def _beschrijf_numeriek(waarden, normaal):
    waarden = waarden.dropna()
    if waarden.empty:
        return "–"
    if normaal:
        return f"{waarden.mean():.1f} ({waarden.std():.1f})"
    q1, med, q3 = waarden.quantile([0.25, 0.5, 0.75])
    return f"{med:.1f} [{q1:.1f}–{q3:.1f}]"


def _toets_numeriek(groepen, normaal):
    groepen = [g.dropna() for g in groepen]
    if len(groepen) < 2 or any(len(g) < config.MIN_N_PER_GROEP for g in groepen):
        return np.nan, ""
    if len(groepen) == 2:
        if normaal:
            return stats.ttest_ind(*groepen, equal_var=False).pvalue, "Welch t-toets"
        return stats.mannwhitneyu(*groepen).pvalue, "Mann-Whitney U"
    if normaal:
        return stats.f_oneway(*groepen).pvalue, "one-way ANOVA"
    return stats.kruskal(*groepen).pvalue, "Kruskal-Wallis"


def _toets_categorisch(kruistabel):
    kruistabel = kruistabel.loc[kruistabel.sum(axis=1) > 0, kruistabel.sum(axis=0) > 0]
    if kruistabel.shape[0] < 2 or kruistabel.shape[1] < 2:
        return np.nan, ""
    chi2, p, _, verwacht = stats.chi2_contingency(kruistabel)
    if kruistabel.shape == (2, 2) and (verwacht < 5).any():
        return stats.fisher_exact(kruistabel).pvalue, "Fisher exact"
    return p, "chi-kwadraat"


def tabel_1(df, groep_kolom=None):
    """
    Maak een 'Tabel 1' met kenmerken van de patiënten.

    Numeriek: gemiddelde (SD) bij normale verdeling, anders mediaan [IQR].
    Categorisch: n (%). Met een groepskolom worden de groepen vergeleken
    en wordt per variabele een p-waarde en de gebruikte toets vermeld.
    """
    groep_kolom = groep_kolom if groep_kolom in df.columns else None
    if groep_kolom:
        niveaus = list(df[groep_kolom].dropna().unique())
        niveaus.sort(key=str)
        groepen = {niv: df[df[groep_kolom] == niv] for niv in niveaus}
    else:
        groepen = {}

    rijen = []

    def rij(variabele, totaal, per_groep=None, p=np.nan, toets="", ontbr=""):
        r = {"variabele": variabele, f"totaal (n={len(df)})": totaal}
        for niv, d in groepen.items():
            r[f"{niv} (n={len(d)})"] = (per_groep or {}).get(niv, "")
        if groepen:
            r["p-waarde"] = _fmt_p(p)
            r["toets"] = toets
        r["ontbrekend"] = ontbr
        rijen.append(r)

    for kolom in config.NUMERIEK:
        if kolom not in df.columns or kolom == groep_kolom:
            continue
        normaal = _is_normaal(df[kolom])
        maat = "gem. (SD)" if normaal else "mediaan [IQR]"
        per_groep = {niv: _beschrijf_numeriek(d[kolom], normaal) for niv, d in groepen.items()}
        p, toets = (_toets_numeriek([d[kolom] for d in groepen.values()], normaal)
                    if groepen else (np.nan, ""))
        rij(f"{kolom}, {maat}", _beschrijf_numeriek(df[kolom], normaal),
            per_groep, p, toets, ontbr=int(df[kolom].isna().sum()))

    for kolom in config.CATEGORISCH:
        if kolom not in df.columns or kolom == groep_kolom:
            continue
        p, toets = (_toets_categorisch(pd.crosstab(df[kolom], df[groep_kolom]))
                    if groepen else (np.nan, ""))
        rij(f"{kolom}, n (%)", "", None, p, toets, ontbr=int(df[kolom].isna().sum()))
        for waarde in sorted(df[kolom].dropna().unique(), key=str):
            def n_pct(d):
                geldig = d[kolom].notna().sum()
                n = (d[kolom] == waarde).sum()
                return f"{n} ({n / geldig * 100:.1f}%)" if geldig else "–"
            rij(f"    {waarde}", n_pct(df),
                {niv: n_pct(d) for niv, d in groepen.items()})

    return pd.DataFrame(rijen)


def correlaties(df):
    """Spearman-correlatiematrix van de numerieke variabelen."""
    kolommen = [k for k in config.NUMERIEK if k in df.columns]
    return df[kolommen].corr(method="spearman").round(2)
