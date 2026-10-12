"""
Volledige pipeline: inlezen → opschonen → analyse → figuren.

Gebruik (vanuit de hoofdmap van het project):
    python -m patient_analyse.main
    python -m patient_analyse.main --bestand pad/naar/data.xlsx --groep behandeling
"""

import argparse
from pathlib import Path

import pandas as pd

from . import analyse, config, figuren, inlezen


def main():
    parser = argparse.ArgumentParser(description="Analyse van een patiëntendataset")
    parser.add_argument("--bestand", default=config.DATA_BESTAND, help="CSV- of Excelbestand")
    parser.add_argument("--groep", default=config.GROEP_KOLOM,
                        help="kolom voor groepsvergelijking ('geen' = uitschakelen)")
    parser.add_argument("--uitvoer", default=config.UITVOER_MAP, help="map voor resultaten")
    args = parser.parse_args()
    groep = None if str(args.groep).lower() in ("geen", "none", "") else args.groep

    uitvoer = Path(args.uitvoer)
    (uitvoer / "figuren").mkdir(parents=True, exist_ok=True)

    print(f"1. Inlezen: {args.bestand}")
    ruw = inlezen.lees_bestand(args.bestand)
    print(f"   {len(ruw)} rijen, {len(ruw.columns)} kolommen")

    print("2. Opschonen")
    df, logboek = inlezen.schoon_op(ruw)
    print(f"   {len(df)} patiënten over, {len(logboek)} wijzigingen in het logboek")
    if not logboek.empty:
        for stap, n in logboek["stap"].value_counts().items():
            print(f"     - {stap}: {n}")

    ontbrekend = inlezen.ontbrekende_waarden(df)

    print("3. Analyse")
    tabel1 = analyse.tabel_1(df, groep)
    corr = analyse.correlaties(df)

    excel = uitvoer / "resultaten.xlsx"
    with pd.ExcelWriter(excel) as schrijver:
        tabel1.to_excel(schrijver, sheet_name="Tabel 1", index=False)
        ontbrekend.to_excel(schrijver, sheet_name="Ontbrekende waarden")
        corr.to_excel(schrijver, sheet_name="Correlaties")
        logboek.to_excel(schrijver, sheet_name="Opschoonlogboek", index=False)
    df.to_csv(uitvoer / "opgeschoonde_data.csv", index=False)

    print("4. Figuren")
    paden = figuren.maak_alle_figuren(df, corr, groep, uitvoer / "figuren")

    print("\nTabel 1:")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(tabel1.to_string(index=False))
    print(f"\nKlaar. Resultaten in '{uitvoer}/':")
    print(f"  - {excel.name} (Tabel 1, ontbrekende waarden, correlaties, opschoonlogboek)")
    print("  - opgeschoonde_data.csv")
    for p in paden:
        print(f"  - figuren/{p.name}")


if __name__ == "__main__":
    main()
