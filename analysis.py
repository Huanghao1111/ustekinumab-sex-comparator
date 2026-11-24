from pathlib import Path
import argparse
import math
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[0]
INDEX = "ustekinumab"
IL23 = {"ustekinumab", "guselkumab", "risankizumab", "tildrakizumab"}
BIOLOGICS = IL23 | {"secukinumab", "ixekizumab", "bimekizumab", "brodalumab",
                    "adalimumab", "certolizumab", "etanercept", "infliximab"}
COMPARATORS = {
    "C4": {"guselkumab", "risankizumab", "tildrakizumab"},
    "C5": {"secukinumab", "ixekizumab"},
    "C6": {"secukinumab", "ixekizumab", "bimekizumab"},
}
OUTCOMES = {
    "depression_pt": {"depression"},
    "anxiety_pt": {"anxiety"},
    "suicidal_ideation_pt": {"suicidal ideation"},
}


def norm(value):
    return re.sub(r"\s+", " ", str(value or "").strip().lower())


def read_frame(data_dir, database):
    base = data_dir / database
    if database == "faers":
        reports = pd.read_csv(base / "reports.csv.gz", usecols=["report_id", "received_date", "sex_standard", "age_y"], dtype={"report_id": str})
        indications = pd.read_csv(base / "indications.csv.gz", usecols=["report_id", "indi_pt"], dtype={"report_id": str}).rename(columns={"indi_pt": "indication_pt"})
    else:
        reports = pd.read_csv(base / "reports.csv.gz", usecols=["report_id", "received_date", "sex", "age_y"], dtype={"report_id": str})
        reports["sex_standard"] = reports["sex"].map({"female": "F", "male": "M"}).fillna("")
        indications = pd.read_csv(base / "indications.csv.gz", usecols=["report_id", "indication_pt"], dtype={"report_id": str})
    drugs = pd.read_csv(base / "drug.csv.gz", usecols=["report_id", "drug_group"], dtype={"report_id": str}).drop_duplicates()
    reactions = pd.read_csv(base / "reactions.csv.gz", usecols=["report_id", "pt"], dtype={"report_id": str}).drop_duplicates()

    reports["received_date"] = pd.to_datetime(reports["received_date"], errors="coerce")
    reports = reports[reports["received_date"].between("2012-10-01", "2026-06-30") & reports["sex_standard"].isin(["F", "M"])].copy()
    reactions["pt"] = reactions["pt"].fillna("").map(norm)
    indications["indication_pt"] = indications["indication_pt"].fillna("").map(norm)

    non_ae = pd.read_csv(ROOT / "dictionaries/non_adverse_event_pt_dictionary.csv")["exact_or_prefix_pt"].dropna().map(norm)
    non_ae_rx = re.compile(r"^(?:" + "|".join(re.escape(x) for x in non_ae) + ")", re.I)
    valid = reactions["pt"].ne("") & ~reactions["pt"].str.contains(non_ae_rx, na=False)
    reports = reports[reports["report_id"].isin(set(reactions.loc[valid, "report_id"]))]

    derm_terms = pd.read_csv(ROOT / "dictionaries/dermatologic_indication_dictionary.csv")["term"].dropna().map(norm)
    derm_rx = re.compile(r"^(?:" + "|".join(re.escape(x) for x in derm_terms) + ")", re.I)
    derm_ids = set(indications.loc[indications["indication_pt"].str.contains(derm_rx, na=False), "report_id"])
    psoriasis_ids = set(indications.loc[indications["indication_pt"].str.contains("psoriasis", case=False, na=False), "report_id"])
    reports["dermatologic_indication"] = reports["report_id"].isin(derm_ids)
    reports["psoriasis_indication"] = reports["report_id"].isin(psoriasis_ids)

    keep = set(reports["report_id"])
    drugs = drugs[drugs["report_id"].isin(keep)]
    reactions = reactions[reactions["report_id"].isin(keep)]
    drug_sets = drugs.groupby("report_id")["drug_group"].agg(frozenset).rename("drugs")
    frame = reports.merge(drug_sets, left_on="report_id", right_index=True, how="inner")
    for outcome, terms in OUTCOMES.items():
        ids = set(reactions.loc[reactions["pt"].isin(terms), "report_id"])
        frame[outcome] = frame["report_id"].isin(ids).astype(np.int8)
    return frame


def study_arms(frame, comparator):
    drugs = frame["drugs"]
    has_index = drugs.map(lambda x: INDEX in x)
    indication = frame["dermatologic_indication"] if comparator == "C1" else frame["psoriasis_indication"]
    target = BIOLOGICS if comparator == "C3" else COMPARATORS.get(comparator)
    if target is None:
        overlap = drugs.map(lambda x: INDEX in x and bool(x - {INDEX}))
        comparator_drugs = drugs.map(lambda x: bool(x - {INDEX}))
    else:
        overlap = drugs.map(lambda x: INDEX in x and bool((x - {INDEX}) & target))
        comparator_drugs = drugs.map(lambda x: bool((x - {INDEX}) & target))
    return has_index & indication & ~overlap, ~has_index & indication & comparator_drugs


def estimate(frame, database, comparator, outcome, adult_only):
    work = frame[frame["age_y"].ge(18)] if adult_only else frame
    index, comparator_mask = study_arms(work, comparator)
    row = {"database": database, "comparator": comparator, "outcome": outcome,
           "scenario": "known_adults" if adult_only else "all_ages"}
    sex_estimates = {}
    for sex in ["F", "M"]:
        i = index & work["sex_standard"].eq(sex)
        c = comparator_mask & work["sex_standard"].eq(sex)
        a, ni = int(work.loc[i, outcome].sum()), int(i.sum())
        cc, nc = int(work.loc[c, outcome].sum()), int(c.sum())
        b, d = ni - a, nc - cc
        row.update({f"a_{sex}": a, f"b_{sex}": b, f"c_{sex}": cc, f"d_{sex}": d,
                    f"n_index_{sex}": ni, f"n_comparator_{sex}": nc})
        if min(a, b, cc, d) > 0:
            log_ror = math.log((a * d) / (b * cc))
            se = math.sqrt(1 / a + 1 / b + 1 / cc + 1 / d)
            sex_estimates[sex] = (log_ror, se)
            row[f"ROR_{sex}"] = math.exp(log_ror)
            row[f"ROR_{sex}_lci"] = math.exp(log_ror - 1.96 * se)
            row[f"ROR_{sex}_uci"] = math.exp(log_ror + 1.96 * se)
        else:
            sex_estimates[sex] = None
            row.update({f"ROR_{sex}": np.nan, f"ROR_{sex}_lci": np.nan, f"ROR_{sex}_uci": np.nan})
    if all(sex_estimates.values()):
        log_rorr = sex_estimates["F"][0] - sex_estimates["M"][0]
        se = math.sqrt(sex_estimates["F"][1] ** 2 + sex_estimates["M"][1] ** 2)
        row.update({"estimable": True, "log_RORR": log_rorr, "SE_log_RORR": se,
                    "RORR": math.exp(log_rorr), "RORR_lci": math.exp(log_rorr - 1.96 * se),
                    "RORR_uci": math.exp(log_rorr + 1.96 * se),
                    "p_value": math.erfc(abs(log_rorr / se) / math.sqrt(2))})
        row["D_index"] = math.exp(math.log((row["a_F"] / row["b_F"]) / (row["a_M"] / row["b_M"])))
        row["D_comparator"] = math.exp(math.log((row["c_F"] / row["d_F"]) / (row["c_M"] / row["d_M"])))
    else:
        row.update({"estimable": False, "log_RORR": np.nan, "SE_log_RORR": np.nan,
                    "RORR": np.nan, "RORR_lci": np.nan, "RORR_uci": np.nan, "p_value": np.nan,
                    "D_index": np.nan, "D_comparator": np.nan})
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=ROOT / "data/processed")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for database in ["faers", "canada"]:
        frame = read_frame(args.data, database)
        for adult in [False, True]:
            for comparator in ["C1", "C2", "C3", "C4", "C5", "C6"]:
                for outcome in OUTCOMES:
                    rows.append(estimate(frame, database, comparator, outcome, adult))
    result = pd.DataFrame(rows)
    result.to_csv(args.out / "core_estimates.csv", index=False)
    print(f"Wrote {len(result)} estimates to {args.out / 'core_estimates.csv'}")


if __name__ == "__main__":
    main()
