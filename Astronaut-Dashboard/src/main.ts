import { Chart, RadarController, RadialLinearScale, PointElement,
         LineElement, Filler, Tooltip, Legend } from 'chart.js'
import type { ScoresJSON, AstronautScores } from './types'
import { DOMAINS, DOMAIN_LABELS, RISK_LABELS, RISK_COLORS, SUBJECTS } from './types'

Chart.register(RadarController, RadialLinearScale, PointElement,
               LineElement, Filler, Tooltip, Legend)

let scores: ScoresJSON = {}
let activeSubject = SUBJECTS[0]
let radarChart: Chart | null = null
let anemiaVisible = true

const ASTRONAUT_SUMMARIES: Record<string, { summary: string; followup: string[] }> = {
  C001: {
    summary: `Your overall health profile looks good. Most areas are showing only minor changes from your pre-flight baseline. The areas with the most activity are DNA damage response and inflammation — both mild. Your body appears to be handling the return to Earth well. One thing to keep an eye on: your red blood cell markers showed a small dip, which is common after spaceflight and typically resolves on its own.`,
    followup: [
      'Schedule a routine blood panel at your next check-up to confirm red blood cell levels have returned to normal',
      'Continue standard post-flight monitoring — your profile does not indicate any urgent concerns',
      'No additional action required at this time',
    ]
  },
  C002: {
    summary: `Your profile shows two areas of elevated activity worth monitoring — DNA damage response and mitochondrial function, both scoring Significant (4/5). These suggest your cells were under meaningful stress during or after the mission. Your immune system and inflammation levels are mild and consistent with normal post-flight recovery. The DNA and mitochondrial signals were picked up across multiple data sources, which makes them more reliable findings.`,
    followup: [
      'Follow up with your flight surgeon to review the DNA damage and mitochondrial findings in more detail',
      'Consider additional blood work at your next scheduled appointment — these signals are worth tracking over time',
      'Prioritize sleep and recovery — mitochondrial stress responds well to rest and reduced physical load in the short term',
    ]
  },
  C003: {
    summary: `Your profile shows the most variation across the crew. DNA damage response and mitochondrial function are both Significant (4/5). You also showed the highest immune regulation score of any crew member — Moderate (3/5) — meaning your immune system showed more shifts from baseline than others. This does not necessarily indicate illness, but suggests your body mounted a stronger response to spaceflight. Worth monitoring in follow-up.`,
    followup: [
      'Discuss your immune regulation scores with your flight surgeon — the elevated response is worth monitoring across the next few timepoints',
      'Follow up on DNA damage and mitochondrial markers at your next scheduled blood draw',
      'If you notice unusual fatigue or difficulty recovering from exercise, flag this early — it may be consistent with what the data is showing',
    ]
  },
  C004: {
    summary: `Your profile is one of the more stable in the crew. Immune regulation is Low (1/5) — the best score across all four crew members in that category. DNA damage response is Moderate (3/5) and mitochondrial function is Significant (4/5), consistent with the rest of the crew in those areas. Overall your inflammatory and immune markers look close to your pre-flight baseline.`,
    followup: [
      'Routine post-flight monitoring is sufficient given your overall stable profile',
      'Follow up on mitochondrial function markers at your next check-up — this was elevated across the whole crew and worth tracking',
      'No urgent action required — your immune and inflammatory markers look close to your pre-flight baseline',
    ]
  },
}

async function loadScores(): Promise<void> {
  const res = await fetch('./scores.json')
  scores = await res.json()
  render()
}

function getRiskColor(score: number): string {
  return RISK_COLORS[score] ?? '#94a3b8'
}

function getRiskLabel(score: number): string {
  return RISK_LABELS[score] ?? 'Unknown'
}

function buildHeader(): string {
  return `
    <header>
      <div class="header-inner">
        <div class="header-left">
          <div class="mission-badge">INSPIRATION4</div>
          <h1>AstronautMD</h1>
          <p class="tagline">Individualized Health Risk Dashboard</p>
        </div>
        <div class="header-right">
          <div class="stat">
            <span class="stat-num">5</span>
            <span class="stat-label">Domains</span>
          </div>
          <div class="stat">
            <span class="stat-num">5</span>
            <span class="stat-label">Datasets</span>
          </div>
          <div class="stat">
            <span class="stat-num">4</span>
            <span class="stat-label">Crew</span>
          </div>
        </div>
      </div>
    </header>
  `
}

function buildTabs(): string {
  return `
    <div class="tab-bar">
      ${SUBJECTS.map(s => `
        <button class="tab-btn ${s === activeSubject ? 'active' : ''}"
                data-subject="${s}">
          <span class="tab-id">${s}</span>
        </button>
      `).join('')}
    </div>
  `
}

function buildScoreCards(astro: AstronautScores): string {
  return `
    <div class="score-grid">
      ${DOMAINS.map(d => {
        const score = astro[d]
        const color = getRiskColor(score.risk_score)
        const label = getRiskLabel(score.risk_score)
        return `
          <div class="score-card" style="--accent: ${color}">
            <div class="score-card-domain">${DOMAIN_LABELS[d]}</div>
            <div class="score-card-number">${score.risk_score}<span class="score-denom">/5</span></div>
            <div class="score-card-label" style="color: ${color}">${label}</div>
            <div class="score-card-sources">${score.sources}</div>
          </div>
        `
      }).join('')}
    </div>
  `
}

function buildNotification(subject: string): string {
  const data = ASTRONAUT_SUMMARIES[subject]
  if (!data) return ''

  if (!anemiaVisible) return `
    <button class="anemia-show-btn" id="anemia-show-btn">
      View your personal health summary
    </button>
  `

  return `
    <div class="anemia-card" id="anemia-card">
      <div class="anemia-header">
        <div class="anemia-title-group">
          <span class="anemia-icon">📋</span>
          <div>
            <div class="anemia-title">Your Post-Flight Health Summary</div>
            <div class="anemia-subtitle">Based on your personal pre-flight baseline across 5 biological domains</div>
          </div>
        </div>
        <button class="anemia-dismiss" id="anemia-dismiss">Dismiss</button>
      </div>
      <div class="anemia-body">
        <p class="anemia-text">${data.summary}</p>
        <div class="followup-section">
          <div class="followup-title">Recommended Follow-up</div>
          <ul class="followup-list">
            ${data.followup.map(item => `<li>${item}</li>`).join('')}
          </ul>
        </div>
        <div class="anemia-finding">
          <span class="anemia-finding-label">⚠ Key Finding:</span>
          Red blood cell production genes were suppressed at landing across the crew, confirmed by blood counts, gene expression, and epigenetic data.
        </div>
      </div>
    </div>
  `
}

function buildPanel(subject: string): string {
  const astro = scores[subject]
  if (!astro) return ''

  return `
    <div class="panel">
      ${buildNotification(subject)}
      <div class="panel-top">
        <div class="panel-left">
          <div class="section-label">Domain Risk Scores</div>
          ${buildScoreCards(astro)}
          <div class="caveat">
            ⚠ n=4 limitation: All scores are within-person comparisons.
            Each astronaut serves as their own baseline. No group-level
            p-values are reported — only deviation magnitudes.
          </div>
        </div>
        <div class="panel-right">
          <div class="section-label">Risk Radar</div>
          <div class="radar-wrap">
            <canvas id="radar-chart"></canvas>
          </div>
        </div>
      </div>

      <div class="panel-bottom">
        <div class="section-label">Data Sources per Domain</div>
        <div class="sources-grid">
          ${DOMAINS.map(d => `
            <div class="source-row">
              <span class="source-domain">${DOMAIN_LABELS[d]}</span>
              <span class="source-tags">
                ${astro[d].sources.split(', ').map(s =>
                  `<span class="tag">${s}</span>`
                ).join('')}
              </span>
              <span class="source-frac">${(astro[d].fraction * 100).toFixed(1)}% markers affected</span>
            </div>
          `).join('')}
        </div>
      </div>
    </div>
  `
}

function buildHeatmap(): string {
  return `
    <div class="heatmap-section">
      <div class="section-label">Crew-Wide Comparison</div>
      <div class="heatmap-grid">
        <div class="heatmap-header-row">
          <div class="heatmap-corner"></div>
          ${DOMAINS.map(d => `
            <div class="heatmap-col-label">${DOMAIN_LABELS[d]}</div>
          `).join('')}
        </div>
        ${SUBJECTS.map(s => {
          const astro = scores[s]
          if (!astro) return ''
          return `
            <div class="heatmap-row">
              <div class="heatmap-row-label">${s}</div>
              ${DOMAINS.map(d => {
                const score = astro[d].risk_score
                const color = getRiskColor(score)
                const label = getRiskLabel(score)
                return `
                  <div class="heatmap-cell" style="background: ${color}22; border: 1px solid ${color}44;"
                       title="${DOMAIN_LABELS[d]}: ${score}/5 (${label})">
                    <span class="heatmap-score" style="color: ${color}">${score}</span>
                    <span class="heatmap-label">${label}</span>
                  </div>
                `
              }).join('')}
            </div>
          `
        }).join('')}
      </div>
    </div>
  `
}

function buildFooter(): string {
  return `
    <footer>
      <span>AstronautMD — Jeremy John, Rohan Pandit, Kai Mayberry</span>
      <span>Data: NASA OSDR OSD-569, OSD-575, OSD-656, GLDS-561</span>
      <span>Torchlight Biosovereignty Hackathon 2026</span>
    </footer>
  `
}

function drawRadar(subject: string): void {
  if (radarChart) {
    radarChart.destroy()
    radarChart = null
  }

  const canvas = document.getElementById('radar-chart') as HTMLCanvasElement
  if (!canvas) return

  const astro = scores[subject]
  if (!astro) return

  const data = DOMAINS.map(d => astro[d].risk_score)
  const labels = DOMAINS.map(d => {
    const words = DOMAIN_LABELS[d].split(' ')
    if (words.length <= 2) return words
    return [
      words.slice(0, Math.ceil(words.length / 2)).join(' '),
      words.slice(Math.ceil(words.length / 2)).join(' ')
    ]
  })

  radarChart = new Chart(canvas, {
    type: 'radar',
    data: {
      labels,
      datasets: [{
        data,
        backgroundColor: 'rgba(96, 165, 250, 0.15)',
        borderColor: 'rgba(96, 165, 250, 0.9)',
        borderWidth: 2,
        pointBackgroundColor: data.map(v => getRiskColor(v)),
        pointBorderColor: '#fff',
        pointBorderWidth: 1.5,
        pointRadius: 5,
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: true,
      scales: {
        r: {
          min: 0,
          max: 5,
          ticks: {
            stepSize: 1,
            color: '#64748b',
            font: { size: 10, family: 'Space Mono' },
            backdropColor: 'transparent',
          },
          grid: { color: 'rgba(255,255,255,0.08)' },
          angleLines: { color: 'rgba(255,255,255,0.08)' },
          pointLabels: {
            color: '#94a3b8',
            font: { size: 11, family: 'DM Sans', weight: 500 },
          }
        }
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: ctx => `Risk: ${ctx.raw}/5 — ${getRiskLabel(ctx.raw as number)}`
          },
          backgroundColor: '#1e293b',
          titleColor: '#e2e8f0',
          bodyColor: '#94a3b8',
          borderColor: '#334155',
          borderWidth: 1,
        }
      }
    }
  })
}

function wireAnemia(): void {
  const dismiss = document.getElementById('anemia-dismiss')
  const show = document.getElementById('anemia-show-btn')
  if (dismiss) {
    dismiss.addEventListener('click', () => {
      anemiaVisible = false
      const container = document.getElementById('panel-container')!
      container.innerHTML = buildPanel(activeSubject)
      drawRadar(activeSubject)
      wireAnemia()
    })
  }
  if (show) {
    show.addEventListener('click', () => {
      anemiaVisible = true
      const container = document.getElementById('panel-container')!
      container.innerHTML = buildPanel(activeSubject)
      drawRadar(activeSubject)
      wireAnemia()
    })
  }
}

function render(): void {
  const app = document.getElementById('app')!
  app.innerHTML = `
    ${buildHeader()}
    <main>
      ${buildTabs()}
      <div id="panel-container">
        ${buildPanel(activeSubject)}
      </div>
      ${buildHeatmap()}
    </main>
    ${buildFooter()}
  `

  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      activeSubject = (btn as HTMLElement).dataset.subject!
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'))
      btn.classList.add('active')
      const container = document.getElementById('panel-container')!
      container.innerHTML = buildPanel(activeSubject)
      drawRadar(activeSubject)
      wireAnemia()
    })
  })

  drawRadar(activeSubject)
  wireAnemia()
}

loadScores()