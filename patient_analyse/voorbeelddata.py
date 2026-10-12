"""
Genereer een FICTIEVE voorbeelddataset om de code mee te testen.

De data zijn willekeurig gegenereerd en bevatten bewust 'vuile' waarden
(komma's als decimaalteken, tikfouten, uitschieters, dubbele patiënten,
ontbrekende waarden) zodat je kunt zien wat het opschonen doet.

Gebruik:  python -m patient_analyse.voorbeelddata
"""

from pathlib import Path

import numpy as np
import pandas as pd


def maak_voorbeelddata(n=300, seed=42):
    rng = np.random.default_rng(seed)
    behandeling = rng.choice(["standaard", "interventie"], n)
    effect = np.where(behandeling == "interventie", -6, 0)

    df = pd.DataFrame({
        "patient_id": [f"P{i:04d}" for i in range(1, n + 1)],
        "leeftijd": rng.normal(64, 12, n).round().clip(18, 99),
        "geslacht": rng.choice(["M", "V"], n),
        "bmi": rng.normal(27.5, 4.5, n).round(1),
        "systolische_bd": (rng.normal(142, 18, n) + effect).round(),
        "diastolische_bd": rng.normal(84, 10, n).round(),
        "hba1c": rng.lognormal(np.log(48), 0.2, n).round(),
        "roker": rng.choice(["ja", "nee"], n, p=[0.25, 0.75]),
        "diagnose": rng.choice(["hartfalen", "COPD", "diabetes", "pneumonie"], n),
        "behandeling": behandeling,
        "opnameduur_dagen": rng.gamma(2, 2.5, n).round(),
        "heropname_30d": rng.choice(["ja", "nee"], n, p=[0.15, 0.85]),
        "opnamedatum": (pd.Timestamp("2024-01-01")
                        + pd.to_timedelta(rng.integers(0, 600, n), unit="D")).strftime("%d-%m-%Y"),
    }).astype(object)

    # Bewust 'vuile' data toevoegen
    def kies(k):
        return rng.choice(n, k, replace=False)

    for i in kies(25):
        df.at[i, "bmi"] = str(df.at[i, "bmi"]).replace(".", ",")   # komma als decimaal
    for i in kies(8):
        df.at[i, "geslacht"] = rng.choice(["man", "vrouw", " m", "F", "onbekend"])
    df.loc[kies(3), "leeftijd"] = [3, 154, 230]                       # onmogelijke leeftijden
    df.loc[kies(2), "systolische_bd"] = [1400, 15]                     # invoerfouten
    df.loc[kies(2), "bmi"] = ["n.v.t.", "abc"]                        # tekst in getalkolom
    df.loc[kies(2), "opnamedatum"] = ["31-02-2024", "01-01-2099"]      # ongeldige datum
    for kolom, k in [("bmi", 15), ("hba1c", 30), ("roker", 12), ("diastolische_bd", 6)]:
        df.loc[kies(k), kolom] = np.nan                               # ontbrekende waarden

    dubbel = df.sample(4, random_state=seed)                          # dubbele patiënten
    return pd.concat([df, dubbel], ignore_index=True)


if __name__ == "__main__":
    pad = Path("data/voorbeeld_patienten.csv")
    pad.parent.mkdir(exist_ok=True)
    maak_voorbeelddata().to_csv(pad, sep=";", index=False)
    print(f"Fictieve voorbeelddataset opgeslagen in {pad}")
