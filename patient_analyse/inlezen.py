"""Inlezen en opschonen van de patiëntendataset."""

from pathlib import Path

import numpy as np
import pandas as pd

from . import config


def lees_bestand(pad):
    """Lees een CSV- of Excelbestand in als DataFrame (alles als tekst)."""
    pad = Path(pad)
    if not pad.exists():
        raise FileNotFoundError(f"Bestand niet gevonden: {pad}")

    if pad.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(pad, sheet_name=config.EXCEL_WERKBLAD, dtype=str)
    else:
        # sep=None laat pandas zelf ; of , als scheidingsteken herkennen
        # (Nederlandse Excel-exports gebruiken vaak ;).
        df = pd.read_csv(pad, sep=None, engine="python", dtype=str,
                         encoding="utf-8-sig")

    df.columns = [str(k).strip().lower().replace(" ", "_") for k in df.columns]
    return df


def _naar_getal(reeks):
    """Zet tekst om naar getallen; accepteert ook een komma als decimaalteken."""
    tekst = reeks.astype(str).str.strip().str.replace(",", ".", regex=False)
    return pd.to_numeric(tekst.replace({"": np.nan, "nan": np.nan}), errors="coerce")


def schoon_op(df):
    """
    Schoon de dataset op volgens de instellingen in config.py.

    Geeft (opgeschoonde DataFrame, logboek als DataFrame) terug. Het logboek
    bevat iedere wijziging, zodat het opschonen reproduceerbaar en
    controleerbaar blijft.
    """
    df = df.copy()
    log = []

    def noteer(stap, kolom, rij, oud, nieuw):
        # rij + 2: regel 1 is de kopregel en pandas telt vanaf 0
        log.append({"stap": stap, "kolom": kolom, "rij_in_bestand": rij + 2,
                    "oude_waarde": oud, "nieuwe_waarde": nieuw})

    # Lege rijen verwijderen
    leeg = df.isna().all(axis=1) | (df.fillna("").apply(lambda r: r.str.strip()) == "").all(axis=1)
    for i in df.index[leeg]:
        noteer("lege rij verwijderd", None, i, None, None)
    df = df[~leeg]

    # Ontbrekende kolommen melden (niet stilzwijgend negeren)
    verwacht = set(config.NUMERIEK) | set(config.CATEGORISCH) | set(config.DATUMS)
    if config.ID_KOLOM:
        verwacht.add(config.ID_KOLOM)
    ontbrekend = sorted(verwacht - set(df.columns))
    if ontbrekend:
        print(f"LET OP: kolommen uit config.py niet gevonden in data: {ontbrekend}")

    # Dubbele patiënten verwijderen (eerste voorkomen blijft staan)
    if config.ID_KOLOM in df.columns:
        df[config.ID_KOLOM] = df[config.ID_KOLOM].str.strip()
        dubbel = df.duplicated(subset=config.ID_KOLOM, keep="first")
        for i in df.index[dubbel]:
            noteer("dubbele patiënt verwijderd", config.ID_KOLOM, i,
                   df.at[i, config.ID_KOLOM], None)
        df = df[~dubbel]
    else:
        dubbel = df.duplicated(keep="first")
        for i in df.index[dubbel]:
            noteer("dubbele rij verwijderd", None, i, None, None)
        df = df[~dubbel]

    # Numerieke kolommen
    for kolom, (minimum, maximum) in config.NUMERIEK.items():
        if kolom not in df.columns:
            continue
        ruw = df[kolom]
        getal = _naar_getal(ruw)

        niet_te_lezen = getal.isna() & ruw.notna() & (ruw.astype(str).str.strip() != "")
        for i in df.index[niet_te_lezen]:
            noteer("geen geldig getal", kolom, i, ruw[i], np.nan)

        buiten = pd.Series(False, index=df.index)
        if minimum is not None:
            buiten |= getal < minimum
        if maximum is not None:
            buiten |= getal > maximum
        for i in df.index[buiten]:
            noteer(f"buiten bereik [{minimum}, {maximum}]", kolom, i, getal[i], np.nan)

        df[kolom] = getal.mask(buiten)

    # Categorische kolommen
    for kolom, hercodering in config.CATEGORISCH.items():
        if kolom not in df.columns:
            continue
        ruw = df[kolom]
        schoon = ruw.astype("string").str.strip()
        sleutel = schoon.str.lower()
        schoon = schoon.mask(sleutel.isin(["", "nan", "onbekend", "?"]))
        if hercodering:
            # Vergelijk zonder hoofdletters; niet-gehercodeerde waarden blijven zoals ze zijn
            schoon = sleutel.map(hercodering).astype("string").fillna(schoon)
        veranderd = schoon.fillna("<NA>") != ruw.astype("string").fillna("<NA>")
        for i in df.index[veranderd]:
            noteer("hercodering", kolom, i, ruw[i], schoon[i])
        df[kolom] = schoon.astype("category")

    # Datumkolommen
    for kolom in config.DATUMS:
        if kolom not in df.columns:
            continue
        ruw = df[kolom]
        datum = pd.to_datetime(ruw, dayfirst=True, errors="coerce", format="mixed")
        mislukt = datum.isna() & ruw.notna() & (ruw.astype(str).str.strip() != "")
        for i in df.index[mislukt]:
            noteer("geen geldige datum", kolom, i, ruw[i], pd.NaT)
        toekomst = datum > pd.Timestamp.today()
        for i in df.index[toekomst]:
            noteer("datum in de toekomst", kolom, i, ruw[i], pd.NaT)
        df[kolom] = datum.mask(toekomst)

    df = df.reset_index(drop=True)
    logboek = pd.DataFrame(log, columns=["stap", "kolom", "rij_in_bestand", "oude_waarde", "nieuwe_waarde"])
    return df, logboek


def ontbrekende_waarden(df):
    """Overzicht van het aantal en percentage ontbrekende waarden per kolom."""
    n = df.isna().sum()
    return pd.DataFrame({
        "ontbrekend_n": n,
        "ontbrekend_%": (n / len(df) * 100).round(1),
    }).sort_values("ontbrekend_n", ascending=False)
