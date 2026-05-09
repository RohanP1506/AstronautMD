# 🚀 AstronautMD — Individualized Astronaut Health Risk Dashboard

**Team:** Jeremy John, Rohan Pandit, Kai Mayberry
**Hackathon:** Torchlight Biosovereignty Hackathon 2026 — Track 2: Individualized Risk Profile

---

## Problem

The Inspiration4 mission generated hundreds of biological measurements per astronaut across 8 timepoints — from 92 days before launch through 194 days after landing. This data is publicly available through NASA's Open Science Data Repository, but it lives across dozens of raw spreadsheets with no integrated interpretation layer.

Existing spaceflight biology analyses report group-level findings. But with only 4 crew members, traditional statistics are underpowered and inappropriate. The real clinical question isn't *what happened to the crew on average* — it's **what happened to this specific astronaut, in this specific biological domain, relative to their own body before flight?**

AstronautMD was built to answer that question.

---

## Approach

We integrated five NASA OSDR datasets spanning the full multi-omics picture of spaceflight health:

| Dataset | Source | What it measures |
|---|---|---|
| Complete Blood Count | OSD-569 | Immune cells, red blood cells, platelets |
| Whole blood RNA-seq | GLDS-561 | Gene expression across biological domains |
| Direct m6A sequencing | GLDS-561 | Epitranscriptomic modification at R+1 |
| Serum cytokines (EVE panel) | OSD-575 | 23 circulating immune signaling proteins |
| Urine cytokines (Alamar panel) | OSD-656 | Urinary inflammation and oxidative stress markers |

Each astronaut's postflight measurements are compared to their own personal preflight mean across L-92, L-44, and L-3 — making every astronaut their own control. Findings are mapped onto five canonical biological domains: **immune regulation, inflammation, oxidative stress, DNA damage response, and mitochondrial function**. Within each domain, markers are weighted by their biological centrality, and the R+1 timepoint is downweighted to avoid conflating acute reentry physiology with true spaceflight-induced pathology. Final domain scores are binned into a 1–5 risk tier per astronaut and surfaced through an interactive dashboard.

---

## Results

### The Strongest Finding: Space Anemia Confirmed Across Three Independent Data Layers

The most compelling result in the dataset is the convergence of signal on a single gene — **HBM (Hemoglobin M)** — across three completely independent measurement modalities at R+1:

- **CBC:** Hemoglobin and RBC counts trend toward or below clinical reference bounds in multiple crew members post-landing
- **RNA-seq:** HBM is significantly downregulated (log2FC −1.47) at R+1 vs preflight baseline
- **m6A sequencing:** HBM shows simultaneous epigenetic modification at R+1

Transcriptional suppression and epitranscriptomic modification of the same erythroid gene, corroborated by clinical blood counts, is the strongest multi-layer biological signal in the dataset and directly supports the space anemia narrative.

### Crew-Wide Domain Risk Scores

DNA damage response and mitochondrial function domains show the highest scores overall, driven by RNA-seq and m6A data at R+1, with notable inter-individual variability — C002 and C003 score higher than C001 and C004, suggesting that the cellular stress response to spaceflight is not uniform across crew members. Serum cytokine inflammation markers remain elevated as late as R+82 and R+194, pointing to prolonged rather than acute post-landing inflammatory activity. Immune regulation scores are generally mild across the crew by R+194, suggesting partial recovery of the adaptive immune arm within six months of landing.

> 📊 Figures, radar charts, and the full interactive dashboard are in the `outputs/` directory.

---

## Impact

Space medicine is entering an era where commercial spaceflight will send people with no astronaut training — and no crew physician — into orbit. The tools to monitor and interpret their health in real time don't yet exist in an accessible form. AstronautMD is a proof of concept that publicly available multi-omics data can be integrated into an individualized, interpretable health risk profile without requiring bioinformatics expertise to read.

The within-person scoring framework scales naturally: every additional dataset or biomarker panel adds another signal to the same five-domain structure. And because scores are built from the fraction of markers showing abnormal values — not black-box model outputs — every score is fully traceable back to the underlying data.

---

## How to Run It

Requires Node.js.

```bash
cd Astronaut-Dashboard
npm install
npm run dev
```