"""Visualisaties van de patiëntendataset (opgeslagen als PNG)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # geen scherm nodig; figuren worden alleen opgeslagen
import matplotlib.pyplot as plt  # noqa: E402

from . import config  # noqa: E402

KLEUR = "#2E4F82"


def _opslaan(fig, map_, naam):
    pad = Path(map_) / f"{naam}.png"
    fig.tight_layout()
    fig.savefig(pad, dpi=150)
    plt.close(fig)
    return pad


def histogrammen(df, map_):
    kolommen = [k for k in config.NUMERIEK if k in df.columns]
    if not kolommen:
        return None
    n_kol = min(3, len(kolommen))
    n_rij = -(-len(kolommen) // n_kol)
    fig, assen = plt.subplots(n_rij, n_kol, figsize=(4 * n_kol, 3.2 * n_rij), squeeze=False)
    for as_, kolom in zip(assen.flat, kolommen):
        as_.hist(df[kolom].dropna(), bins=20, color=KLEUR, edgecolor="white")
        as_.set_title(kolom)
        as_.set_ylabel("aantal patiënten")
    for as_ in list(assen.flat)[len(kolommen):]:
        as_.set_visible(False)
    return _opslaan(fig, map_, "histogrammen")


def boxplots_per_groep(df, groep_kolom, map_):
    if not groep_kolom or groep_kolom not in df.columns:
        return None
    kolommen = [k for k in config.NUMERIEK if k in df.columns]
    niveaus = sorted(df[groep_kolom].dropna().unique(), key=str)
    n_kol = min(3, len(kolommen))
    n_rij = -(-len(kolommen) // n_kol)
    fig, assen = plt.subplots(n_rij, n_kol, figsize=(4 * n_kol, 3.5 * n_rij), squeeze=False)
    for as_, kolom in zip(assen.flat, kolommen):
        data = [df.loc[df[groep_kolom] == niv, kolom].dropna() for niv in niveaus]
        as_.boxplot(data, tick_labels=[str(n) for n in niveaus], patch_artist=True,
                    boxprops={"facecolor": "#D6E0F0"}, medianprops={"color": KLEUR})
        as_.set_title(kolom)
    for as_ in list(assen.flat)[len(kolommen):]:
        as_.set_visible(False)
    fig.suptitle(f"Numerieke variabelen per {groep_kolom}")
    return _opslaan(fig, map_, f"boxplots_per_{groep_kolom}")


def staafdiagrammen(df, map_):
    kolommen = [k for k in config.CATEGORISCH if k in df.columns]
    if not kolommen:
        return None
    n_kol = min(3, len(kolommen))
    n_rij = -(-len(kolommen) // n_kol)
    fig, assen = plt.subplots(n_rij, n_kol, figsize=(4 * n_kol, 3.2 * n_rij), squeeze=False)
    for as_, kolom in zip(assen.flat, kolommen):
        telling = df[kolom].astype(object).fillna("ontbrekend").astype(str).value_counts()
        as_.bar(telling.index, telling.values, color=KLEUR)
        as_.set_title(kolom)
        as_.tick_params(axis="x", rotation=30)
    for as_ in list(assen.flat)[len(kolommen):]:
        as_.set_visible(False)
    return _opslaan(fig, map_, "staafdiagrammen")


def correlatie_heatmap(corr, map_):
    if corr.shape[0] < 2:
        return None
    fig, as_ = plt.subplots(figsize=(1 + 0.9 * len(corr), 0.8 + 0.8 * len(corr)))
    beeld = as_.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
    as_.set_xticks(range(len(corr)), corr.columns, rotation=45, ha="right")
    as_.set_yticks(range(len(corr)), corr.index)
    for i in range(len(corr)):
        for j in range(len(corr)):
            waarde = corr.iat[i, j]
            as_.text(j, i, f"{waarde:.2f}", ha="center", va="center",
                     color="white" if abs(waarde) > 0.6 else "black", fontsize=8)
    fig.colorbar(beeld, ax=as_, label="Spearman ρ")
    as_.set_title("Correlaties tussen numerieke variabelen")
    return _opslaan(fig, map_, "correlaties")


def ontbrekende_data(df, map_):
    pct = df.isna().mean().sort_values() * 100
    fig, as_ = plt.subplots(figsize=(6, 0.35 * len(pct) + 1))
    as_.barh(pct.index, pct.values, color=KLEUR)
    as_.set_xlabel("% ontbrekend")
    as_.set_xlim(0, max(10, pct.max() * 1.1))
    as_.set_title("Ontbrekende waarden per variabele")
    return _opslaan(fig, map_, "ontbrekende_waarden")


def maak_alle_figuren(df, corr, groep_kolom, map_):
    paden = [
        histogrammen(df, map_),
        boxplots_per_groep(df, groep_kolom, map_),
        staafdiagrammen(df, map_),
        correlatie_heatmap(corr, map_),
        ontbrekende_data(df, map_),
    ]
    return [p for p in paden if p]
