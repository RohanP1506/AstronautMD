# 🚀 AstronautMD — Individualized Astronaut Health Risk Dashboard

**Team:** Jeremy John, Rohan Pandit, Kai Mayberry
**Hackathon:** Torchlight Biosovereignty Hackathon 2026 — Track 2: Individualized Risk Profile

---

## Problem

The Inspiration4 mission generated hundreds of biological measurements per astronaut multiple timepoints, before and after the flight. With only 4 crew members, we only showcased the ability to focus on what happened to this specific astronaut, in this specific biological domain, relative to their own body before flight.

---

## Approach

We integrated five NASA OSDR datasets spanning the multi-omics data of spaceflight health:

| Dataset | Source | What it measures |
|---|---|---|
| Complete Blood Count | OSD-569 | Immune cells, red blood cells, platelets |
| Whole blood RNA-seq | GLDS-561 | Gene expression across biological domains |
| Direct m6A sequencing | GLDS-561 | Epitranscriptomic modification at R+1 |
| Serum cytokines (EVE panel) | OSD-575 | 23 circulating immune signaling proteins |
| Urine cytokines (Alamar panel) | OSD-656 | Urinary inflammation and oxidative stress markers |

Each astronaut's postflight measurements are compared to their own personal preflight mean across L-92, L-44, and L-3 — making every astronaut their own control. Findings are mapped onto five canonical biological domains: immune regulation, inflammation, oxidative stress, DNA damage response, and mitochondrial function. Within each domain, markers are weighted by their importance to each domain, and the R+1 timepoint is downweighted to avoid construed data through physiological stress of landing. Final domain scores are put into a 1–5 risk tier per astronaut and surfaced through an interactive dashboard.

---

## Scoring

Each data source contributes domain scores through the same pipeline:

1. **Compute log2FC** relative to each astronaut's personal preflight mean (L-92, L-44, L-3). CBC uses clinical reference ranges instead: each analyte is flagged as abnormal if it falls outside the `RANGE_MIN`/`RANGE_MAX` bounds in the dataset, since CBC markers are on incomparable physical scales across analytes.
2. **Flag significant markers**: |log2FC| ≥ 1 (2x fold-change) for cytokines, RNA-seq, and urine; absolute shift ≥ 10 % for m6A, which is in probability space rather than ratio space.
3. **Compute a weighted fraction** of significant markers per domain per timepoint. Markers are weighted 1.0 (primary driver of the domain) or 0.5 (secondary contributor).
4. **Weight timepoints** — R+1 is weighted 0.5x, all later timepoints 1.0x, to reduce the influence of acute reentry physiology on the final score.
5. **Average across data sources** — all sources contributing to a domain are averaged at the fraction level before binning, so the final score reflects convergent signal across modalities rather than any single dataset.
6. **Place into 1–5 risk tiers:**

| Fraction of markers abnormal | Risk Score | Label |
|---|---|---|
| < 10% | 1 | Low |
| 10–25% | 2 | Mild |
| 25–50% | 3 | Moderate |
| 50–75% | 4 | Significant |
| > 75% | 5 | High |

---

## Results

### The Strongest Finding: Space Anemia Signal Across Three Independent Data Layers

The most compelling result in the dataset is the convergence of signal on a single gene — **HBM (Hemoglobin M)** — across three completely independent measurement modalities at R+1:

- **RNA-seq:** HBM is significantly downregulated (log2FC −1.47) at R+1 vs preflight baseline
- **m6A sequencing:** HBM shows simultaneous epigenetic modification at R+1
- **CBC:** Hemoglobin and RBC counts trend toward or below clinical reference bounds in multiple crew members post-landing

Transcriptional suppression and epitranscriptomic modification of the same erythroid gene, corroborated by clinical blood counts, is the strongest multi-layer biological signal in the dataset and directly supports the space anemia narrative.

### Crew-Wide Domain Risk Scores

DNA damage response and mitochondrial function domains show the highest scores overall, driven by RNA-seq and m6A data at R+1, with notable variability between individuals — C002 and C003 score higher than C001 and C004, suggesting that the cellular stress response to spaceflight is not uniform across crew members. Serum cytokine inflammation markers remain elevated as late as R+82 and R+194, pointing to prolonged rather than acute post-landing inflammatory activity. Immune regulation scores are generally mild across the crew by R+194, suggesting partial recovery of the adaptive immune arm within six months of landing.

> 📊 Figures, radar charts, and the full interactive dashboard are in the `outputs/` directory.

---

## How to Run It

Requires Node.js.

```bash
cd Astronaut-Dashboard
npm install
npm run dev
```

To regenerate scores from raw NASA data:

```bash
pip install pandas numpy matplotlib seaborn openpyxl
python download_data.py
python run_all_scores.py
```