# AstronautMD

**An individualized dashboard for measuring astronaut health**

**Team:** Jeremy John, Rohan Pandit, Kai Mayberry

---

## 🚀 Problem Statement
*   Astronauts experience changes across all bodily systems during spaceflight, but raw multi-omics data is too complex for astronauts to interpret and make informed decisions on.
*   The Inspiration4 mission generated hundreds of values per astronaut of multi-omics data across 8 timepoints, but no integrated risk-scoring framework currently presents this data in a form a non-geneticist astronaut or space medicine specialist can read and act on.
*   Existing analyses report changes across groups of people — they do not give individual crew members a personalized snapshot of their own physiological response across key markers of physical health.
*   **Goal:** Produce a transparent, defensible, individualized health risk profile for each Inspiration4 crew member without overstating findings given $n=4$ limitations.

---

## 🛠 Approach

### Two-Part Scoring Framework
*   **Part A — Per-astronaut scores:** Score each crew member individually across 5 biological domains using datasets given (cytokine arrays, urine inflammation panel, metabolic panel, CBC, RNA-seq raw counts).
*   **Part B — Crew-level evidence:** Use pre-computed differential expression results (plasma metabolomics, plasma proteomics, EVP proteomics) as plausible indicators for Part A.

### Five Biological Domains
*   **Immune regulation:** Primary dataset: cytokine arrays + PBMC.
*   **Inflammation:** Primary dataset: urine inflammation panel + cytokines.
*   **Oxidative stress:** Primary dataset: plasma metabolomics.
*   **DNA damage response:** Primary dataset: whole blood RNA-seq.
*   **Mitochondrial function:** Primary dataset: metabolomics + RNA-seq mitochondrial gene set.

---

## 📊 Scoring Methodology
*   Score by how far each marker deviates from the population normal ranges.
*   Personal baseline $log_2$ fold-change scoring for marker-level data (each astronaut's pre-flight average as baseline).
*   Pathway-level evidence indicators from crew-wide differential expression.
*   **Consistency weighting:** Sustained differences across multiple timepoints contribute more than sudden spikes.
*   **Risk tiers:** Per domain based on score: **Low, Medium, High**.

---

## 🛡️ n=4 Mitigation Strategy
*   **Within-person comparisons only:** Each astronaut serves as their own baseline.
*   **No group-level p-value claims:** Only report magnitude.
*   **Individual Presentation:** All 4 astronauts presented individually on dashboard to show any variability.
*   **Confidence indicator:** Crew-level differential expression used as confidence indicator, not the primary score.

---

## 💻 Expected Output
Web-based interactive dashboard with one risk summary page per astronaut (features can change):
*   Scorecard for each astronaut per domain.
*   Radar chart with preflight vs inflight vs postflight for each domain.
*   Per-domain timeline plots showing change over time in metrics for all astronauts.
*   Confidence indicators derived from crew-level evidence. (maybe make an algorithm for this)
*   Documentation explaining each scoring threshold and why it was chosen.
*   Acknowledgment of $n=4$ limitations on each score.
