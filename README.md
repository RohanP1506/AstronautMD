# 🚀 AstronautMD
**An individualized multi-omics health dashboard for astronaut recovery monitoring**

**Team:** Jeremy John · Rohan Pandit · Kai Mayberry

## Project Summary

AstronautMD is an individualized bioinformatics dashboard that translates raw multi-omics spaceflight data into transparent, actionable health risk profiles for each crew member. Using data from NASA's Inspiration4 mission (OSDR OSD-569, OSD-575, OSD-656), we score four astronauts across five biological domains — immune regulation, inflammation, anemia risk, mitochondrial stress, and proteostasis — and present results as a personalized, time-resolved scorecard that a non-specialist can interpret.

The core design principle: **every astronaut serves as their own control.** All scores are derived from within-person fold-changes against a personal preflight baseline, never group statistics.

---

## Motivation & Problem Statement

The Inspiration4 mission generated one of the richest multi-omics datasets ever collected in a civilian spaceflight context, spanning cytokine arrays, urine proteomics, plasma metabolomics, whole-blood RNA-seq, complete blood counts, and direct RNA m6A methylation — hundreds of values per person across eight timepoints.

Despite this richness, **the data are not currently interpretable at the individual level.** Existing analyses report crew-wide averages or differential expression statistics that collapse across all four astronauts. An individual crew member or their flight surgeon cannot look at a fold-change table and answer: *"How did my immune system respond? Is it back to normal? What should I watch?"*

AstronautMD addresses this gap by:
- Computing personalized risk scores per biological domain per timepoint
- Presenting scores on a per-astronaut dashboard, not a group summary
- Contextualizing individual findings against crew-level differential expression as a confidence indicator
- Explicitly communicating statistical limitations at the n=4 scale alongside every score

---

## Data Sources

| Dataset | OSDR ID | Modality | Used For |
|---|---|---|---|
| Serum multiplex cytokines (EVE Panel) | OSD-575 | Protein (multiplex) | Immune regulation · Inflammation |
| Urine immune/inflammatory proteomics (Alamar Panel) | OSD-656 | Protein (proximity extension) | Inflammation · Oxidative stress |
| Complete blood count | OSD-569 | Clinical labs | Anemia corroboration · NLR |
| Whole-blood RNA-seq | OSD-569 / GLDS-561 | Transcriptomics | Anemia risk · Immune dysregulation · Mitochondrial stress · Proteostasis · Transcriptional stress |
| Direct RNA m6A methylation | OSD-569 / GLDS-561 | Epitranscriptomics | Post-transcriptional regulation (in progress) |
| Plasma metabolomics DEGs | Pre-computed | Metabolomics | Crew-level confidence indicator |
| Plasma / EVP proteomics DEGs | Pre-computed | Proteomics | Crew-level confidence indicator |

---

## Scoring Framework

### Design Principles

AstronautMD uses a scoring framework based on individual-level evidence Domain composite scores are the primary output. They are computed from per-astronaut log₂ fold-changes against each person's own preflight mean. No group statistics are used or implied.

---

### Part A: Cytokine & Urine Scoring (OSD-575 / OSD-656)

Used for **immune regulation** and **inflammation** domains. Produces a time-resolved trajectory across all postflight timepoints.

#### Step 1 — Personal Baseline
Each astronaut's preflight mean is computed across all pre-launch timepoints (L-92, L-44, L-3) per marker.

#### Step 2 — Log₂ Fold-Change
```
log₂FC = log₂(postflight concentration / personal preflight mean)
```
Significance threshold: **|log₂FC| ≥ 1.0** (≥ 2x fold-change).

#### Step 3 — Domain Score (Weighted Fraction)
```
weighted_fraction = Σ(weights of significant markers) / Σ(weights of all domain markers)
```
Marker weights reflect biological centrality: **1.0** for primary domain drivers, **0.5** for secondary contributors.

> **Note on IFN-γ:** IFN-γ is included in the immune regulation domain (serum cytokines) but excluded from the RNA-seq immune dysregulation domain to avoid double-counting its contribution across modalities, which would artificially inflate cross-domain correlation.

#### Step 4 — Timepoint Weighting
R+1 is down-weighted to **0.5×** to reduce confounding from acute reentry physiological stress (G-force, fluid shifts, sleep disruption). All other postflight timepoints are weighted **1.0×**.

#### Step 5 — Risk Tier

| Tier | Label | Weighted fraction of markers affected |
|---|---|---|
| 1 | Minimal | < 10% |
| 2 | Mild | 10–25% |
| 3 | Moderate | 25–50% |
| 4 | Significant | 50–75% |
| 5 | Severe | > 75% |

---

### Part A: RNA-seq Scoring (OSD-569 / GLDS-561)

Used for **anemia risk, immune dysregulation, mitochondrial stress, proteostasis, and transcriptional stress** domains. Only the R+1 timepoint is available, so this produces a **snapshot score**, not a trajectory.

#### Why a Different Scoring Method?
Cytokine scoring uses the fraction of markers significant (bounded 0–1), which works because each domain contains ~10–13 carefully selected markers. RNA-seq domains contain hundreds of candidate genes — the fraction significant would approach zero regardless of effect size. Instead, RNA-seq domains are scored using the **weighted mean |log₂FC| of significant DEGs**, which captures both the breadth and magnitude of dysregulation.

#### Step 1 — Personal Baseline (same principle)
Per-astronaut preflight mean computed from all L- timepoints available in the RNA-seq data (L-92, L-44, L-3 for each astronaut).

#### Step 2 — Significance Filter (dual gate)
A gene is significant if it passes **both**:
- Personal |log₂FC| ≥ 1.0 (within-person 2x fold-change)
- DESeq2 adjusted p-value < 0.05 (crew-level statistical reliability)

This dual gate uses crew-level statistics only to filter noise, not to score — the score itself is always within-person.

#### Step 3 — Domain Assignment
Significant DEGs are matched against curated gene sets for each domain by gene symbol. Gene sets are derived from g:Profiler enrichment results run on the RNA-seq DEG list.

| RNA-seq Domain | Biological Focus | Gene Set Examples | Domain Weight |
|---|---|---|---|
| **Anemia Risk** | RBC biology, erythropoiesis, heme/iron metabolism | HBB, HBA1/2, ALAS2, KLF1, TFRC, GYPA/B | 1.0 |
| **Immune Dysregulation** | Adaptive immunity, cytotoxic T/NK cells, B cells, Tregs | CD3D/E/G, GZMB, PRF1, CD19, FOXP3, KLRD1 | 1.0 |
| **Mitochondrial Stress** | OXPHOS complexes, electron transport chain, membrane transport | MT-CO1/2/3, MT-ND1/2, UQCRC1/2, COX4I1 | 0.8 |
| **Proteostasis** | Ubiquitin-proteasome system, protein quality control, chaperones | UBB, PSMA1-2, PSMB1-3, HSPA1A, HSP90AA1 | 0.8 |
| **Transcriptional Stress** | Stress-responsive TFs, immediate early genes, chromatin regulators | EGR1-3, FOS/FOSB, JUN/JUNB, ATF3/4, NR4A1-3 | 0.6 |

#### Step 4 — Risk Tier (RNA-seq)

| Tier | Label | Mean |log₂FC| of significant domain genes |
|---|---|---|
| 1 | Minimal | 0 significant genes |
| 2 | Mild | > 0 genes, mean \|log₂FC\| < 1.5 |
| 3 | Moderate | mean \|log₂FC\| 1.5–2.0 |
| 4 | Significant | mean \|log₂FC\| 2.0–3.0 |
| 5 | Severe | mean \|log₂FC\| > 3.0 |

---

### CBC: Corroborating Evidence, Not a Scored Domain

Complete blood count data is presented as **raw clinical values with reference range bands**, not as log₂FC. This follows standard medical advice: a hemoglobin of 11.2 g/dL is clinically meaningful; a log₂FC of −0.3 is not. CBC is used in two ways:

1. **Anemia corroboration:** Hemoglobin and RBC count below the lower clinical reference bound at postflight timepoints provide direct clinical validation for the RNA-seq Anemia Risk domain score.
2. **NLR (Neutrophil-to-Lymphocyte Ratio):** A validated clinical stress index plotted over the full pre- and postflight trajectory. NLR > 3 = elevated; NLR > 5 = markedly elevated. Used as a corroborating inflammatory marker alongside cytokine domain scores.

CBC does **not** contribute a primary domain score to the AstronautMD dashboard — it provides supporting clinical context for RNA-seq and cytokine findings.

---

## Preliminary Results

 Serum cytokine domains (immune regulation, inflammation) are fully scored. Urine Alamar composite scores are complete. RNA-seq scoring framework is implemented and pending data integration. CBC visualizations and anemia corroboration table are complete.

---

### Serum Cytokines — Immune Regulation Domain

Postflight immune regulatory profiles are **relatively the same** across the four crew members. No uniform group-level response is observed.

**C001** shows early elevation in IL-10 and IL-13 at R+1, significant suppression at R+45, and a rebound in IL-2 by R+82. Profile largely normalizes by R+194.

**C002** exhibits the most pronounced and persistent dysregulation. IL-21 and IL-4 remain strongly suppressed at R+82 and R+194, suggesting prolonged disruption to T follicular helper and Th2 signaling. Composite marker score remains below baseline through the final timepoint.

**C003** shows an early IL-2 spike at R+1 that resolves rapidly — the regulatory cytokine profile is largely restored by R+45, the fastest apparent recovery in this domain.

**C004** shows late-onset elevation in IL-10, IL-13, and IL-4 at R+194, consistent with a delayed Th2 shift emerging well into recovery.

---

### Serum Cytokines — Inflammation Domain

**C001** shows suppression of MCP-1, IL-8, and IP-10 at R+45 and R+82, suggesting a period of dampened innate immune recruitment during mid-recovery.

**C002** displays a paradoxical inflammatory spike at R+45 (IL-8, IP-10 elevated), followed by a prolonged suppressive trend through R+194. The log₂FC trajectory crosses baseline twice — inconsistent with simple resolution.

**C003** maintains relatively modest inflammatory shifts throughout, with gradual upregulation persisting at R+82 and R+194.

**C004** is the most stable in this domain, remaining close to preflight levels across all postflight timepoints.

---

### Urine Biomarkers — OSD-656 Alamar Panel

Urine inflammation composite scores remain low (near baseline) for most astronauts postflight. Oxidative stress markers show substantially more inter-individual variation, with C002 and C004 showing marked fluctuation even in the preflight window. Persistent marker analysis identifies C002 as having the broadest sustained response (PTX3, MMP8, MMP9, MPO, S100A12, HGF across all three postflight timepoints), consistent with ongoing neutrophil activation and extracellular matrix remodeling.

---

### Key Observations

1. **Immune recovery is individual-specific.** No two astronauts follow the same trajectory in direction, magnitude, or timing.
2. **C002 shows the most sustained immune disruption** across both serum and urine modalities.
3. **Suppression, not elevation, dominates late recovery** for inflammatory markers in C001 and C002 — potentially compensatory immunosuppression, or genuine resolution; the RNA-seq data will help distinguish these.
4. **Urine oxidative stress markers capture a different biological layer** than serum cytokines, with tissue-level signals (AGER, HGF, MMPs) not reflected in the multiplex serum panel.
5. **The R+1 timepoint is an unreliable signal** for several markers due to acute reentry physiology — consistent with our decision to down-weight it in scoring.

---

## n=4 Limitations & Mitigation

| Limitation | Mitigation |
|---|---|
| No statistical power for group comparisons | Within-person log₂FC only; no group p-values reported |
| Single-timepoint RNA-seq (R+1 only) | Scored as an explicit snapshot; clearly labeled as such on dashboard; CBC and cytokine trajectories provide longitudinal context |
| Individual variability indistinguishable from noise | Significance requires \|log₂FC\| ≥ 1.0 AND (for RNA-seq) DESeq2 padj < 0.05 |
| Confounds at R+1 (reentry stress) | R+1 down-weighted to 0.5× in cytokine/urine domain scores; CBC NLR presented in full temporal context |
| No independent replication cohort | Crew-level DEG data used as a directional confidence indicator only, not as primary evidence |
| n=4 limits ability to distinguish responders from non-responders | All four astronauts presented individually; no grouping or averaging across crew |

Every score displayed on the AstronautMD dashboard includes an explicit n=4 caveat and a plain-language description of what the score measures and does not claim.

---

## Codebase Structure

```
AstronautMD/
│
├── scoring/
│   ├── cytokine_dataframe.py     # Core: data loading, personal baseline, log₂FC, domain scoring
│   ├── immune_eve.py             # OSD-575 EVE Panel — domain definitions, marker weights, scoring
│   │                             #   Note: IFN-γ included here, excluded from RNA-seq immune domain
│   ├── composite_score.py        # Domain scoring engine for urine Alamar panel
│   ├── urine_dataframes.py       # Urine data loading and preprocessing
│   ├── urine.py                  # OSD-656 Alamar Panel — analysis entry point
│   └── rnaseq_scoring.py         # RNA-seq domain scoring (5 domains, R+1 snapshot)
│                                 #   Domains: Anemia Risk, Immune Dysregulation,
│                                 #   Mitochondrial Stress, Proteostasis, Transcriptional Stress
│
├── visualization/
│   ├── plots.py                  # Serum cytokine heatmaps and time-course plots
│   └── urine_viz.py              # Urine figures: composite timeline, heatmap, persistence bars
│
├── clinical/
│   └── CBC.py                    # Raw CBC values with clinical reference ranges
│                                 #   NLR trajectory · Anemia corroboration table
│
├── genomics/
│   ├── RNAseq.ipynb              # DEG analysis and g:Profiler enrichment
│   └── m6A.ipynb                 # Direct RNA m6A methylation (OSD-569 / GLDS-561)
│
└── README.md
```

---

## Figures

| Filename | Description | Status |
|---|---|---|
| `fig_log2fc_heatmap_immune_regulation.png` | Serum log₂FC heatmap — immune regulation | ✅ Complete |
| `fig_log2fc_heatmap_inflammation.png` | Serum log₂FC heatmap — inflammation | ✅ Complete |
| `immune_regulation_time_plot.png` | Composite score + log₂FC over time — immune regulation | ✅ Complete |
| `inflammation_time_plot.png` | Composite score + log₂FC over time — inflammation | ✅ Complete |
| `fig1_composite_timeline.png` | Urine domain composite scores over time | ✅ Complete |
| `fig2_log2fc_heatmap.png` | Urine Alamar log₂FC heatmap (all markers) | ✅ Complete |
| `fig3_persistence_bars.png` | Persistent postflight immune shifts by astronaut | ✅ Complete |
| `fig_cbc_raw_values.png` | CBC raw values with clinical reference range bands | ✅ Complete |
| `fig_cbc_nlr.png` | NLR trajectory pre- and post-flight | ✅ Complete |
| `fig_rnaseq_domain_scores.png` | RNA-seq domain risk score heatmap (5 domains × 4 astronauts) | 🔄 Pending data integration |
| `fig_rnaseq_deg_counts.png` | Significant DEG counts per domain at R+1 | 🔄 Pending data integration |
| Dashboard radar chart | Per-domain risk scores per astronaut | 🔄 In progress |

---

## Changelog

**v2 (Checkpoint #2, May 8)**
- `immune_eve.py` — removed IFN-γ from RNA-seq immune domain cross-reference to eliminate double-counting across modalities
- `CBC.py` — replaced log₂FC plots with raw clinical values and reference range bands; added anemia corroboration table and NLR trajectory with clinical threshold annotations
- `rnaseq_scoring.py` — new file: five RNA-seq domains scored at R+1 using DESeq2 padj + personal log₂FC dual gate; separate scoring methodology from cytokines (weighted mean |log₂FC| of significant DEGs rather than fraction of markers); includes heatmap and DEG count visualizations

---

## Roadmap

**Checkpoint #2 (May 8) — Complete**
- [x] Serum cytokine scoring — immune regulation + inflammation
- [x] Urine Alamar composite scores — inflammation + oxidative stress
- [x] CBC raw value plots with reference ranges, NLR, anemia table
- [x] RNA-seq scoring framework implemented (5 domains)
- [x] Preliminary figures — heatmaps, timelines, persistence bars
- [x] Draft README / BioBrief structure

**Next Steps**
- [ ] Run `rnaseq_scoring.py` against GLDS-561 data to produce domain scores
- [ ] Plasma metabolomics integration for oxidative stress domain
- [ ] Crew-level confidence indicator algorithm (Part B)
- [ ] m6A epitranscriptomics cross-reference with RNA-seq DEGs
- [ ] Dashboard frontend — per-astronaut scorecard, radar chart, domain timeline plots
- [ ] Threshold documentation and n=4 limitation callouts embedded per domain on dashboard

---

*Data sourced from the NASA Open Science Data Repository (OSDR). This analysis is preliminary and descriptive. No clinical claims are made.*
