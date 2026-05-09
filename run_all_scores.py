"""
run_all_scores.py — AstronautMD Master Scoring Pipeline
Torchlight Biosovereignty Hackathon — Track 2: Individualized Risk Profile

Runs all domain scoring modules and produces:
  outputs/scores_all_sources.csv   — raw scores per astronaut/domain/source
  outputs/scores_summary.csv       — aggregated final scores per astronaut/domain
  outputs/scores.json              — JSON for dashboard
  outputs/radar_C001-C004.png      — per-astronaut radar charts
  outputs/domain_heatmap.png       — crew-wide domain heatmap
  outputs/dashboard.html           — standalone interactive dashboard

Usage:
    python run_all_scores.py
"""

import json
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')
os.makedirs('outputs', exist_ok=True)

# ── Import scoring modules ──────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'observations', 'urine'))

from cbc_score        import score_cbc
from rnaseq_score     import score_rnaseq
from m6a_score        import score_m6a
from cytokine_score   import score_cytokines
from observations.urine.composite_score  import calculate_composite_scores, DOMAIN_MAP, DOMAIN_WEIGHTS
from observations.urine.urine_dataframes import organize_urine_df, process_urine_df

# ── Domain configuration ────────────────────────────────────────────────────────
CANONICAL_DOMAINS = ['immune_regulation', 'inflammation', 'oxidative_stress',
                     'dna_damage', 'mitochondrial']

DOMAIN_LABELS = {
    'immune_regulation': 'Immune\nRegulation',
    'inflammation':      'Inflammation',
    'oxidative_stress':  'Oxidative\nStress',
    'dna_damage':        'DNA Damage\nResponse',
    'mitochondrial':     'Mitochondrial\nFunction',
}

DOMAIN_LABELS_SHORT = {
    'immune_regulation': 'Immune Regulation',
    'inflammation':      'Inflammation',
    'oxidative_stress':  'Oxidative Stress',
    'dna_damage':        'DNA Damage Response',
    'mitochondrial':     'Mitochondrial Function',
}

RISK_COLORS = {1: '#27ae60', 2: '#82c341', 3: '#f39c12', 4: '#e67e22', 5: '#e74c3c'}
RISK_LABELS = {1: 'Low', 2: 'Mild', 3: 'Moderate', 4: 'Significant', 5: 'High'}

SUBJECTS = ['C001', 'C002', 'C003', 'C004']

# ── 1. CBC pipeline ─────────────────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: CBC")
print("=" * 60)
cbc_raw = pd.read_csv('data/cbc.csv')
cbc_scores = score_cbc(cbc_raw, subject_col='SUBJECT_ID', timepoint_col='TEST_DATE')

# ── 2. RNA-seq scoring ──────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: RNA-seq")
print("=" * 60)

degs_path = os.path.join('observations', 'whole_blood', 'degs_FP1.csv')
if not os.path.exists(degs_path):
    print(f"  WARNING: {degs_path} not found. Skipping RNA-seq scoring.")
    print("  Run RNAseq.ipynb and save with: degs.to_csv('degs_FP1.csv')")
    rnaseq_scores = pd.DataFrame()
else:
    rnaseq_scores = score_rnaseq(degs_path)

# ── 3. m6A scoring ──────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: m6A")
print("=" * 60)

m6a = pd.read_excel(
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?source=datamanager'
    '&file=GLDS-561_directm6Aseq_Direct_RNA_seq_m6A_Processed_Data.xlsx',
    sheet_name='I4-FP1',
    skiprows=[0, 1, 2, 3, 4, 5, 6, 7],
    header=[0, 1],
    index_col=0
)
m6a.columns = [f"{str(a)}_{str(b)}".strip("_") for a, b in m6a.columns]

m6a_filtered = m6a[
    (m6a['I4-FP1_methylKit q-value'] < 0.05) &
    (m6a['I4-FP1_methylKit methylation difference, %'].abs() > 10)
].copy()

print(f"  m6A significant positions: {len(m6a_filtered)}")
m6a_scores = score_m6a(m6a_filtered)

# ── 4. Cytokine scoring ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Cytokine (OSD-575 EVE Panel)")
print("=" * 60)

try:
    cytokine_scores = score_cytokines()
except Exception as e:
    print(f"  WARNING: Cytokine scoring failed: {e}")
    cytokine_scores = pd.DataFrame()

# ── 5. Urine scoring ────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 5: Urine")
print("=" * 60)

urine_path = os.path.join('observations', 'urine', 'urine_data.csv')
if not os.path.exists(urine_path):
    print(f"  WARNING: {urine_path} not found. Skipping urine scoring.")
    urine_scores = pd.DataFrame()
else:
    urine_raw  = pd.read_csv(urine_path, index_col=0)
    urine_long = organize_urine_df(urine_raw.transpose())
    urine_proc = process_urine_df(urine_long, 'urine_alamar_656')
    _, urine_overall = calculate_composite_scores(urine_proc, DOMAIN_MAP, DOMAIN_WEIGHTS)
    urine_overall['source'] = 'Urine'
    urine_scores = urine_overall[['astronaut', 'domain', 'overall_fraction',
                                   'risk_score', 'source']].copy()
    print(f"  Urine scores: {len(urine_scores)} rows")

# ── 6. Combine all scores ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 6: Combining scores")
print("=" * 60)

all_parts = [df for df in [cbc_scores, rnaseq_scores, m6a_scores,
                             cytokine_scores, urine_scores] if len(df) > 0]
all_scores = pd.concat(all_parts, ignore_index=True)
all_scores.to_csv('outputs/scores_all_sources.csv', index=False)
print(f"  Total score rows: {len(all_scores)}")

summary = (
    all_scores[all_scores['domain'].isin(CANONICAL_DOMAINS)]
    .groupby(['astronaut', 'domain'])
    .agg(
        mean_fraction=('overall_fraction', 'mean'),
        n_sources=('source', 'count'),
        sources=('source', lambda x: ', '.join(sorted(x.unique()))),
    )
    .reset_index()
)

def frac_to_risk(f):
    if f < 0.10:   return 1
    elif f < 0.25: return 2
    elif f < 0.50: return 3
    elif f < 0.75: return 4
    else:           return 5

summary['final_risk'] = summary['mean_fraction'].apply(frac_to_risk)
summary.to_csv('outputs/scores_summary.csv', index=False)

print("\nFinal domain risk scores:")
pivot = summary.pivot(index='astronaut', columns='domain', values='final_risk')
print(pivot.to_string())

# ── 7. Radar charts ─────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 7: Radar charts")
print("=" * 60)

def radar_chart(astronaut, summary_df, output_path):
    N      = len(CANONICAL_DOMAINS)
    labels = [DOMAIN_LABELS[d] for d in CANONICAL_DOMAINS]
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    scores = []
    for d in CANONICAL_DOMAINS:
        row = summary_df[(summary_df['astronaut'] == astronaut) & (summary_df['domain'] == d)]
        scores.append(float(row['final_risk'].values[0]) if len(row) > 0 else 0.0)
    scores_plot = scores + scores[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    ax.fill(angles, [5] * (N + 1), alpha=0.06, color='#e74c3c')
    ax.fill(angles, [3] * (N + 1), alpha=0.06, color='#f39c12')
    ax.fill(angles, [2] * (N + 1), alpha=0.06, color='#27ae60')

    for r in [1, 2, 3, 4, 5]:
        ax.plot(angles, [r] * (N + 1), color='#cccccc', linewidth=0.5)

    ax.plot(angles, scores_plot, 'o-', linewidth=2.0, color='#2980b9', markersize=6)
    ax.fill(angles, scores_plot, alpha=0.25, color='#2980b9')

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=9, fontweight='bold')
    ax.set_ylim(0, 5)
    ax.set_yticks([1, 2, 3, 4, 5])
    ax.set_yticklabels(['1\nLow', '2', '3\nMod', '4', '5\nHigh'],
                       fontsize=6.5, color='gray')
    ax.set_title(f'Crew Member {astronaut}\nBiological Domain Risk Profile',
                 fontsize=12, fontweight='bold', pad=22)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")

for subject in SUBJECTS:
    radar_chart(subject, summary, f'outputs/radar_{subject}.png')

# ── 8. Domain heatmap ───────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 8: Domain heatmap")
print("=" * 60)

heatmap_rows = []
for subject in SUBJECTS:
    for domain in CANONICAL_DOMAINS:
        row     = summary[(summary['astronaut'] == subject) & (summary['domain'] == domain)]
        risk    = int(row['final_risk'].values[0]) if len(row) > 0 else 0
        sources = row['sources'].values[0] if len(row) > 0 else 'N/A'
        heatmap_rows.append({'astronaut': subject, 'domain': domain,
                              'risk': risk, 'sources': sources})

heatmap_df = pd.DataFrame(heatmap_rows)
pivot_data  = heatmap_df.pivot(index='astronaut', columns='domain', values='risk')
pivot_data  = pivot_data.reindex(columns=CANONICAL_DOMAINS)

fig, ax = plt.subplots(figsize=(12, 5))
sns.heatmap(
    pivot_data,
    ax=ax,
    cmap='RdYlGn_r',
    vmin=1, vmax=5,
    annot=True, fmt='.0f', annot_kws={'size': 14, 'weight': 'bold'},
    linewidths=1.5, linecolor='white',
    cbar_kws={'label': 'Risk Score (1 = Low, 5 = High)', 'shrink': 0.8},
    xticklabels=[DOMAIN_LABELS_SHORT[d] for d in CANONICAL_DOMAINS],
)
ax.set_title(
    'Individualized Biological Domain Risk Scores — Inspiration4 Crew\n'
    'Integrated from CBC, Serum Cytokines, Urine Cytokines, Whole Blood RNA-seq, and m6A',
    fontsize=12, fontweight='bold', pad=12
)
ax.set_xlabel('')
ax.set_ylabel('Crew Member', fontsize=11)
ax.tick_params(axis='x', rotation=20, labelsize=10)
ax.tick_params(axis='y', rotation=0)

plt.tight_layout()
plt.savefig('outputs/domain_heatmap.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: outputs/domain_heatmap.png")

# ── 9. Export scores.json ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 9: Exporting scores.json")
print("=" * 60)

scores_json = {}
for subject in SUBJECTS:
    astro_df = summary[summary['astronaut'] == subject]
    domain_scores = {}
    for domain in CANONICAL_DOMAINS:
        row = astro_df[astro_df['domain'] == domain]
        if len(row) > 0:
            domain_scores[domain] = {
                'risk_score': int(row['final_risk'].values[0]),
                'fraction':   round(float(row['mean_fraction'].values[0]), 3),
                'n_sources':  int(row['n_sources'].values[0]),
                'sources':    row['sources'].values[0],
            }
        else:
            domain_scores[domain] = {
                'risk_score': 0,
                'fraction':   0.0,
                'n_sources':  0,
                'sources':    'N/A',
            }
    scores_json[subject] = domain_scores

with open('outputs/scores.json', 'w') as f:
    json.dump(scores_json, f, indent=2)
print("  Saved: outputs/scores.json")

# ── 10. Generate dashboard.html ──────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 10: Generating dashboard")
print("=" * 60)

def risk_badge(score, domain):
    color = RISK_COLORS.get(score, '#95a5a6')
    label = RISK_LABELS.get(score, 'N/A')
    return (
        f'<div class="domain-card" style="border-left: 5px solid {color};">'
        f'  <div class="domain-name">{DOMAIN_LABELS_SHORT.get(domain, domain)}</div>'
        f'  <div class="risk-score" style="color:{color};">{score}/5</div>'
        f'  <div class="risk-label" style="color:{color};">{label}</div>'
        f'</div>'
    )

astronaut_panels = ''
for i, subject in enumerate(SUBJECTS):
    scores      = scores_json[subject]
    badges      = ''.join(risk_badge(scores[d]['risk_score'], d) for d in CANONICAL_DOMAINS)
    radar_path  = f'radar_{subject}.png'
    display     = 'block' if i == 0 else 'none'
    sources_list = '<br>'.join(
        f"<b>{DOMAIN_LABELS_SHORT.get(d, d)}:</b> {scores[d]['sources']}"
        for d in CANONICAL_DOMAINS if scores[d]['sources'] != 'N/A'
    )
    astronaut_panels += f'''
    <div id="panel-{subject}" class="astronaut-panel" style="display:{display};">
      <div class="panel-header">
        <h2>Crew Member {subject}</h2>
        <p class="subtitle">Individualized health risk profile — Inspiration4 mission</p>
      </div>
      <div class="section-title">Domain Risk Scores</div>
      <div class="domain-grid">{badges}</div>
      <div class="two-col">
        <div>
          <div class="section-title">Radar Profile</div>
          <img src="{radar_path}" class="figure-img" alt="Radar chart for {subject}">
        </div>
        <div>
          <div class="section-title">Data Sources per Domain</div>
          <div class="sources-box">{sources_list}</div>
          <div class="caveat-box">
            ⚠ <b>n=4 limitation:</b> Scores are within-person comparisons only.
            Each astronaut serves as their own baseline. No group-level p-values
            are reported — only deviation magnitudes. RNA-seq per-astronaut scores
            are magnitude estimates from single-sample comparisons.
          </div>
        </div>
      </div>
    </div>
    '''

scores_js = json.dumps(scores_json)

tab_buttons = ''
for i, s in enumerate(SUBJECTS):
    active = ' active' if i == 0 else ''
    tab_buttons += f'<button class="tab-btn{active}" onclick="showAstronaut(\'{s}\')">{s}</button>\n    '

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>AstronautMD — Inspiration4 Health Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
  :root {{
    --bg: #0f1117;
    --card: #1a1d27;
    --border: #2a2d3a;
    --text: #e8eaf0;
    --muted: #8a8fa8;
    --accent: #4a9eff;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: 'Segoe UI', system-ui, sans-serif; min-height: 100vh; }}
  header {{ background: linear-gradient(135deg, #1a1d27 0%, #0f1117 100%);
            border-bottom: 1px solid var(--border); padding: 20px 40px; }}
  header h1 {{ font-size: 1.8rem; color: var(--accent); }}
  header p  {{ color: var(--muted); font-size: 0.9rem; margin-top: 4px; }}
  .container {{ max-width: 1200px; margin: 0 auto; padding: 30px 20px; }}
  .tab-bar {{ display: flex; gap: 8px; margin-bottom: 24px; flex-wrap: wrap; }}
  .tab-btn {{
    padding: 10px 24px; border-radius: 8px; border: 1px solid var(--border);
    background: var(--card); color: var(--muted); cursor: pointer;
    font-size: 1rem; font-weight: 600; transition: all 0.2s;
  }}
  .tab-btn:hover  {{ border-color: var(--accent); color: var(--accent); }}
  .tab-btn.active {{ background: var(--accent); color: white; border-color: var(--accent); }}
  .astronaut-panel {{ animation: fadeIn 0.3s ease; }}
  @keyframes fadeIn {{ from {{ opacity:0; transform:translateY(8px); }} to {{ opacity:1; transform:none; }} }}
  .panel-header {{ margin-bottom: 20px; }}
  .panel-header h2 {{ font-size: 1.5rem; }}
  .subtitle {{ color: var(--muted); font-size: 0.9rem; margin-top: 4px; }}
  .section-title {{
    font-size: 0.75rem; font-weight: 700; letter-spacing: 0.1em;
    color: var(--muted); text-transform: uppercase; margin-bottom: 12px; margin-top: 24px;
  }}
  .domain-grid {{ display: grid; grid-template-columns: repeat(5, 1fr); gap: 12px; }}
  .domain-card {{
    background: var(--card); border-radius: 10px; padding: 16px 12px;
    border-left: 5px solid #555; text-align: center;
  }}
  .domain-name  {{ font-size: 0.78rem; color: var(--muted); margin-bottom: 8px; min-height: 36px; }}
  .risk-score   {{ font-size: 2rem; font-weight: 800; }}
  .risk-label   {{ font-size: 0.75rem; font-weight: 600; margin-top: 2px; }}
  .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px; margin-top: 8px; }}
  .figure-img {{ width: 100%; border-radius: 10px; border: 1px solid var(--border); }}
  .sources-box {{
    background: var(--card); border-radius: 10px; padding: 16px;
    font-size: 0.82rem; line-height: 1.8; border: 1px solid var(--border); color: var(--muted);
  }}
  .caveat-box {{
    background: #1e1a0f; border: 1px solid #5a4a10; border-radius: 10px;
    padding: 14px; font-size: 0.8rem; color: #c8a84b; margin-top: 12px; line-height: 1.6;
  }}
  .heatmap-section {{ margin-top: 40px; }}
  .heatmap-section img {{ width: 100%; border-radius: 10px; border: 1px solid var(--border); }}
  footer {{
    border-top: 1px solid var(--border); padding: 20px 40px; color: var(--muted);
    font-size: 0.8rem; text-align: center; margin-top: 40px;
  }}
  @media (max-width: 768px) {{
    .domain-grid {{ grid-template-columns: repeat(2, 1fr); }}
    .two-col {{ grid-template-columns: 1fr; }}
  }}
</style>
</head>
<body>
<header>
  <h1>🚀 AstronautMD</h1>
  <p>Individualized health risk dashboard — Inspiration4 mission &nbsp;|&nbsp; Torchlight Biosovereignty Hackathon</p>
</header>
<div class="container">
  <div class="tab-bar">
    {tab_buttons}
  </div>
  {astronaut_panels}
  <div class="heatmap-section">
    <div class="section-title">Crew-Wide Comparison</div>
    <img src="domain_heatmap.png" alt="Domain heatmap all astronauts">
  </div>
</div>
<footer>
  AstronautMD &mdash; Jeremy John, Rohan Pandit, Kai Mayberry &nbsp;|&nbsp;
  Data: NASA OSDR OSD-569, OSD-575, GLDS-561 &nbsp;|&nbsp;
  Torchlight Biosovereignty Hackathon 2026
</footer>
<script>
const SCORES = {scores_js};
const DOMAINS = {json.dumps(CANONICAL_DOMAINS)};
const DOMAIN_LABELS = {json.dumps(DOMAIN_LABELS_SHORT)};

function showAstronaut(subject) {{
  document.querySelectorAll('.astronaut-panel').forEach(p => p.style.display = 'none');
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById('panel-' + subject).style.display = 'block';
  document.querySelectorAll('.tab-btn').forEach(b => {{
    if (b.textContent.trim() === subject) b.classList.add('active');
  }});
}}
</script>
</body>
</html>"""

with open('outputs/dashboard.html', 'w', encoding='utf-8') as f:
    f.write(html)
print("  Saved: outputs/dashboard.html")

print("\n" + "=" * 60)
print("✓ All outputs complete!")
print("  outputs/scores_all_sources.csv")
print("  outputs/scores_summary.csv")
print("  outputs/scores.json")
print("  outputs/radar_C001-C004.png")
print("  outputs/domain_heatmap.png")
print("  outputs/dashboard.html")
print("\nTo view the dashboard, open outputs/ in a terminal and run:")
print("  python -m http.server 8000")
print("  Then open: http://localhost:8000/dashboard.html")
print("=" * 60)
