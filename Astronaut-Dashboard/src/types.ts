export interface DomainScore {
  risk_score: number
  fraction: number
  n_sources: number
  sources: string
}

export interface AstronautScores {
  immune_regulation: DomainScore
  inflammation: DomainScore
  oxidative_stress: DomainScore
  dna_damage: DomainScore
  mitochondrial: DomainScore
}

export type ScoresJSON = Record<string, AstronautScores>

export const DOMAINS = [
  'immune_regulation',
  'inflammation',
  'oxidative_stress',
  'dna_damage',
  'mitochondrial',
] as const

export type Domain = typeof DOMAINS[number]

export const DOMAIN_LABELS: Record<Domain, string> = {
  immune_regulation: 'Immune Regulation',
  inflammation:      'Inflammation',
  oxidative_stress:  'Oxidative Stress',
  dna_damage:        'DNA Damage Response',
  mitochondrial:     'Mitochondrial Function',
}

export const RISK_LABELS: Record<number, string> = {
  1: 'Low',
  2: 'Mild',
  3: 'Moderate',
  4: 'Significant',
  5: 'High',
}

export const RISK_COLORS: Record<number, string> = {
  1: '#4ade80',
  2: '#a3e635',
  3: '#facc15',
  4: '#fb923c',
  5: '#f87171',
}

export const SUBJECTS = ['C001', 'C002', 'C003', 'C004']