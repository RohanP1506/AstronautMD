import pandas as pd

# Specific markers found in the OSD-656 Alamar Panel
# Use the short gene names; the function will find them in the long strings
DOMAIN_MAP = {
    'inflammation': [
        'il6', 'il1b', 'tnf', 'ccl2', 'ccl25', 'cxcl2', 'mmp3', 'mmp8', 'ptx3'
    ],
    'oxidative_stress': [
        'ager', 'gzmb', '8ohdg', 'hgf'
    ]
}

def calculate_composite_scores(processed_df, domain_mapping):
    # 1. Create a flattened lookup for case-insensitive matching
    def find_domain(marker_name):
        marker_clean = str(marker_name).lower()
        for domain, keywords in domain_mapping.items():
            for kw in keywords:
                # This matches 'il6' inside 'il6_concentration_npq'
                if kw.lower() in marker_clean:
                    return domain
        return None

    df_domains = processed_df.copy()
    df_domains['domain'] = df_domains['marker'].apply(find_domain)
    
    # 2. Filter out non-domain markers
    df_domains = df_domains.dropna(subset=['domain'])
    
    if df_domains.empty:
        print("!!! Still no matches. Check if keywords match the marker names below:")
        print(processed_df['marker'].unique()[:5])
        return pd.DataFrame()

    # 3. Aggregate: Mean of marker_score, then scale 0-3 -> 0-5
    summary = df_domains.groupby(['astronaut', 'timepoint', 'domain'])['marker_score'].mean().reset_index()
    summary['composite_score'] = (summary['marker_score'] / 3) * 5
    
    return summary.round(2)