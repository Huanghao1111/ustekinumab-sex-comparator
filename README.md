# Ustekinumab sex-comparator analysis

Case/non-case analysis of ustekinumab depression, anxiety, and suicidal ideation reporting in FAERS and Canada Vigilance.

## Data

Download `processed_cohorts.zip` from the accompanying Mendeley Data dataset and extract its `faers/` and `canada/` folders to `data/processed/`. The source data are de-identified public adverse-event reports from:

- [FAERS quarterly extracts](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html), 2012 Q4–2026 Q2
- [Canada Vigilance](https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database.html)

## Run

```bash
pip install -r requirements.txt
python analysis.py
python bh_adjust.py
python figure.py
```

The analysis window is 2012-10-01 through 2026-06-30. FAERS uses primary/secondary suspect reports; Canada Vigilance uses suspect reports. Analyses include all ages and a known-adult (≥18 years) sensitivity analysis. No automatic 0.5 continuity correction is used.

Comparators: C1, all dermatologic drugs; C2, psoriasis-indication drugs; C3, psoriasis biologics; C4, guselkumab/risankizumab/tildrakizumab; C5, secukinumab/ixekizumab; C6, C5 plus bimekizumab. C2–C6 use psoriasis-indication reports and exclude index/comparator overlap; brodalumab is excluded from the active class comparators.

Outputs are written to `results/`: `core_estimates.csv`, `core_estimates_bh.csv`, and Figure 1 in PDF/PNG format.
