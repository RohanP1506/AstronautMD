import { Chart, RadarController, RadialLinearScale, PointElement,
         LineElement, Filler, Tooltip, Legend } from 'chart.js'
import type { ScoresJSON, AstronautScores } from './types'
import { DOMAINS, DOMAIN_LABELS, RISK_LABELS, RISK_COLORS, SUBJECTS } from './types'

Chart.register(RadarController, RadialLinearScale, PointElement,
               LineElement, Filler, Tooltip, Legend)

let scores: ScoresJSON = {}
let activeSubject = SUBJECTS[0]
let radarChart: Chart | null = null

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

function buildPanel(subject: string): string {
  const astro = scores[subject]
  if (!astro) return ''

  return `
    <div class="panel">
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
  return [words.slice(0, Math.ceil(words.length / 2)).join(' '),
          words.slice(Math.ceil(words.length / 2)).join(' ')]
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

  // Tab click handlers
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      activeSubject = (btn as HTMLElement).dataset.subject!
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'))
      btn.classList.add('active')
      const container = document.getElementById('panel-container')!
      container.innerHTML = buildPanel(activeSubject)
      drawRadar(activeSubject)
    })
  })

  drawRadar(activeSubject)
}

loadScores()
