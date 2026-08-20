import streamlit as st
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import subprocess
import tempfile
import os
import re
import time
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, Select

st.set_page_config(
    page_title="M.tb Pocket Analyser",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');
html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }
.stApp { background-color: #0D1117; color: #E6EDF3; }
.hero { background: linear-gradient(135deg, #0D1117 0%, #161B22 100%); border-bottom: 1px solid #21262D; padding: 2rem 0 1.5rem 0; margin-bottom: 2rem; }
.hero-title { font-family: 'IBM Plex Mono', monospace; font-size: 1.6rem; font-weight: 500; color: #58A6FF; letter-spacing: -0.02em; margin: 0; }
.hero-sub { font-size: 0.9rem; color: #8B949E; margin-top: 0.4rem; font-weight: 300; }
.card { background: #161B22; border: 1px solid #21262D; border-radius: 8px; padding: 1.25rem 1.5rem; margin-bottom: 1rem; }
.card-title { font-family: 'IBM Plex Mono', monospace; font-size: 0.75rem; color: #58A6FF; text-transform: uppercase; letter-spacing: 0.08em; margin-bottom: 0.75rem; }
.metric-row { display: flex; gap: 1.5rem; flex-wrap: wrap; }
.metric { flex: 1; min-width: 100px; }
.metric-val { font-family: 'IBM Plex Mono', monospace; font-size: 1.6rem; font-weight: 500; color: #E6EDF3; line-height: 1; }
.metric-label { font-size: 0.75rem; color: #8B949E; margin-top: 0.25rem; }
.badge { display: inline-block; padding: 2px 10px; border-radius: 99px; font-size: 0.72rem; font-weight: 500; margin-right: 6px; }
.badge-good { background: #0D4429; color: #3FB950; border: 1px solid #238636; }
.badge-warn { background: #2D1F00; color: #E3B341; border: 1px solid #9E6A03; }
.badge-bad { background: #3D0B0B; color: #F85149; border: 1px solid #DA3633; }
.badge-pathway { background: #0C2D6B; color: #58A6FF; border: 1px solid #1F6FEB; }
.stTextInput > div > div > input { background-color: #161B22 !important; border: 1px solid #30363D !important; border-radius: 6px !important; color: #E6EDF3 !important; font-family: 'IBM Plex Mono', monospace !important; font-size: 0.9rem !important; }
.stButton > button { background: #238636 !important; color: #fff !important; border: none !important; border-radius: 6px !important; font-weight: 500 !important; padding: 0.5rem 1.5rem !important; }
.stButton > button:hover { background: #2EA043 !important; }
.status-line { font-family: 'IBM Plex Mono', monospace; font-size: 0.8rem; color: #8B949E; padding: 0.2rem 0; }
.status-ok { color: #3FB950; }
.status-warn { color: #E3B341; }
.status-err { color: #F85149; }
.preset-label { font-size: 0.75rem; color: #8B949E; margin-bottom: 0.4rem; font-family: 'IBM Plex Mono', monospace; }
[data-testid="stSidebar"] { background: #161B22 !important; border-right: 1px solid #21262D; }
</style>
""", unsafe_allow_html=True)

KNOWN = {
    'P9WIA3': ('inhA',  'FAS-II',   'Enoyl-ACP reductase — primary isoniazid target'),
    'P9WJA3': ('kasA',  'FAS-II',   'Beta-ketoacyl-ACP synthase I'),
    'P9WJA1': ('kasB',  'FAS-II',   'Beta-ketoacyl-ACP synthase II — top ranked target'),
    'P9WIS3': ('mabA',  'FAS-II',   '3-ketoacyl reductase'),
    'P9WI77': ('hadA',  'FAS-II',   'Hydroxyacyl-ACP dehydratase subunit A'),
    'P9WI75': ('hadB',  'FAS-II',   'Hydroxyacyl-ACP dehydratase subunit B'),
    'P9WKG3': ('fabH',  'FAS-II',   '3-oxoacyl-ACP synthase III'),
    'P9WKR5': ('fabD',  'FAS-II',   'Malonyl-CoA ACP transacylase'),
    'P9WNE9': ('gyrA',  'DNA-rep',  'DNA gyrase subunit A'),
    'P9WNF1': ('gyrB',  'DNA-rep',  'DNA gyrase subunit B'),
    'P0A5R3': ('dnaN',  'DNA-rep',  'DNA polymerase III beta clamp'),
    'P9WPZ9': ('echA6', 'lipid',    'Enoyl-CoA hydratase'),
    'P9WK71': ('accA3', 'lipid',    'Acyl-CoA carboxylase alpha'),
    'P00374': ('DHFR',  'human',    'Human dihydrofolate reductase — control'),
}

PATHWAY_COLORS = {
    'FAS-II':  '#1F6FEB',
    'DNA-rep': '#D29922',
    'lipid':   '#238636',
    'human':   '#DA3633',
    'unknown': '#8B949E',
}

def fetch_alphafold(accession):
    api = f'https://alphafold.ebi.ac.uk/api/prediction/{accession}'
    r = requests.get(api, timeout=15)
    if r.status_code != 200:
        return None, f'AlphaFold API returned {r.status_code}'
    data = r.json()
    if not data:
        return None, 'No structure found for this ID'
    pdb_url = data[0].get('pdbUrl')
    if not pdb_url:
        return None, 'No PDB URL in response'
    r2 = requests.get(pdb_url, timeout=30)
    if r2.status_code != 200:
        return None, f'PDB download failed: {r2.status_code}'
    return r2.content, None

def fetch_uniprot(accession):
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        gene = ''
        if data.get('genes'):
            gene = data['genes'][0].get('geneName', {}).get('value', '')
        protein_name = (data.get('proteinDescription', {})
                           .get('recommendedName', {})
                           .get('fullName', {})
                           .get('value', 'Unknown'))
        organism = data.get('organism', {}).get('scientificName', 'Unknown')
        length = data.get('sequence', {}).get('length', 0)
        return {'gene': gene, 'protein_name': protein_name,
                'organism': organism, 'length': length}
    except:
        return {'gene': '', 'protein_name': 'Unknown',
                'organism': 'Unknown', 'length': 0}

class HighConfidenceSelect(Select):
    def __init__(self, confident_ids):
        self.confident_ids = confident_ids
    def accept_residue(self, residue):
        return residue.get_id()[1] in self.confident_ids

def filter_pdb(pdb_bytes, accession, threshold=70, window=5):
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='wb') as f:
        f.write(pdb_bytes)
        raw_path = f.name
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(accession, raw_path)
    scores = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    scores[residue.get_id()[1]] = residue['CA'].get_bfactor()
    if not scores:
        return None, 0, 0, []
    resnums = sorted(scores.keys())
    plddt = np.array([scores[r] for r in resnums])
    n_start = 0
    for i in range(len(plddt) - window + 1):
        if np.mean(plddt[i:i+window]) >= threshold:
            n_start = i
            break
    c_end = len(plddt)
    for i in range(len(plddt) - window, -1, -1):
        if np.mean(plddt[i:i+window]) >= threshold:
            c_end = i + window
            break
    confident = {r for r in resnums[n_start:c_end] if scores[r] >= threshold}
    filtered_path = raw_path.replace('.pdb', '_filtered.pdb')
    io = PDBIO()
    io.set_structure(structure)
    io.save(filtered_path, HighConfidenceSelect(confident))
    os.unlink(raw_path)
    return filtered_path, len(scores), len(confident), list(plddt)

def run_fpocket(pdb_path):
    stem = Path(pdb_path).stem
    tmp_pdb = f'/tmp/fp_{stem}.pdb'
    os.system(f'cp "{pdb_path}" "{tmp_pdb}"')
    subprocess.run(['fpocket', '-f', tmp_pdb], capture_output=True, text=True, timeout=120)
    out_dir = f'/tmp/fp_{stem}_out'
    if Path(out_dir).exists():
        return out_dir
    return None

def parse_pockets(out_dir):
    folder = Path(out_dir)
    stem = folder.stem.replace('_out', '')
    info = folder / f'{stem}_info.txt'
    if not info.exists():
        return []
    pockets, cur = [], {}
    with open(info) as f:
        for line in f:
            line = line.strip()
            m = re.match(r'Pocket\s+(\d+)\s*:', line)
            if m:
                if cur: pockets.append(cur)
                cur = {'pocket_num': int(m.group(1))}
                continue
            for key, patterns in [
                ('drug_score',     ['Drug Score']),
                ('volume',         ['Real volume', 'Volume']),
                ('hydrophobicity', ['Hydrophobicity score']),
                ('polarity',       ['Polarity score']),
                ('n_spheres',      ['Number of alpha sphere']),
            ]:
                if any(p in line for p in patterns) and key not in cur:
                    v = re.search(r'([-+]?[0-9]*\.?[0-9]+)', line.split(':')[-1])
                    if v:
                        cur[key] = float(v.group(1))
    if cur: pockets.append(cur)
    pockets.sort(key=lambda x: x.get('drug_score', 0), reverse=True)
    return pockets

def install_fpocket():
    if subprocess.run(['which', 'fpocket'], capture_output=True).returncode == 0:
        return True
    os.system('apt-get install -y -q libnetcdf-dev 2>/dev/null')
    os.system('git clone -q https://github.com/Discngine/fpocket.git /tmp/fpocket_src 2>/dev/null')
    os.system('cd /tmp/fpocket_src && make -s 2>/dev/null && cp bin/fpocket /usr/local/bin/ 2>/dev/null')
    os.system('chmod +x /usr/local/bin/fpocket 2>/dev/null')
    return subprocess.run(['which', 'fpocket'], capture_output=True).returncode == 0

def plddt_plot(plddt_list):
    fig, ax = plt.subplots(figsize=(10, 2.8))
    fig.patch.set_facecolor('#161B22')
    ax.set_facecolor('#161B22')
    colors = []
    for s in plddt_list:
        if s >= 90:   colors.append('#0053D6')
        elif s >= 70: colors.append('#65CBF3')
        elif s >= 50: colors.append('#FFDB13')
        else:         colors.append('#FF7D45')
    ax.bar(range(len(plddt_list)), plddt_list, color=colors, width=1.0, linewidth=0)
    ax.axhline(70, color='#F85149', linestyle='--', linewidth=1.2)
    ax.axhspan(0, 70, alpha=0.06, color='#F85149')
    ax.set_xlabel('Residue position', fontsize=9, color='#8B949E')
    ax.set_ylabel('pLDDT', fontsize=9, color='#8B949E')
    ax.set_ylim(0, 100)
    ax.tick_params(colors='#8B949E', labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor('#30363D')
    from matplotlib.patches import Patch
    patches = [
        Patch(color='#0053D6', label='≥90 very high'),
        Patch(color='#65CBF3', label='70–90 confident'),
        Patch(color='#FFDB13', label='50–70 low'),
        Patch(color='#FF7D45', label='<50 very low'),
    ]
    ax.legend(handles=patches, loc='lower right', fontsize=7,
              facecolor='#21262D', edgecolor='#30363D', labelcolor='#8B949E')
    plt.tight_layout()
    return fig

def pocket_bar_plot(pockets):
    if not pockets:
        return None
    top5 = pockets[:min(5, len(pockets))]
    nums = [f"#{p['pocket_num']}" for p in top5]
    scores = [p.get('drug_score', 0) for p in top5]
    vols = [p.get('volume', 0) for p in top5]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3))
    fig.patch.set_facecolor('#161B22')
    for ax in [ax1, ax2]:
        ax.set_facecolor('#161B22')
        ax.tick_params(colors='#8B949E', labelsize=8)
        for spine in ax.spines.values():
            spine.set_edgecolor('#30363D')
    bar_colors = ['#238636' if s >= 0.5 else '#E3B341' if s >= 0.3 else '#8B949E' for s in scores]
    ax1.bar(nums, scores, color=bar_colors, edgecolor='#0D1117', linewidth=0.5)
    ax1.set_title('Druggability score', fontsize=9, color='#8B949E', pad=6)
    ax1.set_ylim(0, 1)
    ax1.axhline(0.5, color='#238636', linestyle=':', linewidth=1, alpha=0.6)
    ax1.set_ylabel('Score (0–1)', fontsize=8, color='#8B949E')
    ax2.bar(nums, vols, color='#1F6FEB', edgecolor='#0D1117', linewidth=0.5)
    ax2.set_title('Pocket volume (Å³)', fontsize=9, color='#8B949E', pad=6)
    ax2.set_ylabel('Volume', fontsize=8, color='#8B949E')
    plt.tight_layout()
    return fig

# Header
st.markdown("""
<div class="hero">
  <p class="hero-title">🧬 M.tb Pocket Analyser</p>
  <p class="hero-sub">AlphaFold structure → pLDDT filtering → binding pocket detection — in one step</p>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.markdown("### About this tool")
    st.markdown("""
This tool implements the three-phase pipeline from the **AlphaFold Pocket Comparison** project:

1. Fetch predicted structure from AlphaFold DB
2. Filter low-confidence regions (pLDDT < 70)
3. Detect binding pockets with fpocket

Built as part of a BSc computational biology portfolio project analysing *M. tuberculosis* drug targets at IISER Bhopal (2026).
    """)
    st.markdown("---")
    st.markdown("**Top targets from full analysis:**")
    st.markdown("🥇 kasB — selectivity 0.711")
    st.markdown("🥈 fabH — selectivity 0.546")
    st.markdown("🥉 fabD — selectivity 0.422")

col_input, col_info = st.columns([1.2, 1])

with col_input:
    st.markdown('<div class="card"><div class="card-title">Enter UniProt Accession ID</div>', unsafe_allow_html=True)
    st.markdown('<div class="preset-label">Quick select from study proteins:</div>', unsafe_allow_html=True)
    preset_cols = st.columns(4)
    presets = [('kasB','P9WJA1'), ('fabH','P9WKG3'), ('inhA','P9WIA3'), ('DHFR','P00374')]
    chosen_preset = None
    for i, (label, acc) in enumerate(presets):
        with preset_cols[i]:
            if st.button(label, key=f'preset_{acc}'):
                chosen_preset = acc
    accession_input = st.text_input('UniProt ID', value=chosen_preset or '',
                                     placeholder='e.g. P9WJA1', label_visibility='collapsed')
    run = st.button('Analyse structure →', use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)

with col_info:
    if accession_input and accession_input.strip():
        acc = accession_input.strip().upper()
        if acc in KNOWN:
            gene, pathway, desc = KNOWN[acc]
            st.markdown(f"""
<div class="card">
  <div class="card-title">Known protein</div>
  <div style="font-size:1.3rem;font-weight:600;color:#E6EDF3;font-family:'IBM Plex Mono',monospace">{gene}</div>
  <div style="color:#8B949E;font-size:0.85rem;margin:4px 0 8px 0">{desc}</div>
  <span class="badge badge-pathway">{pathway}</span>
</div>
""", unsafe_allow_html=True)
        else:
            st.markdown(f"""
<div class="card">
  <div class="card-title">Custom protein</div>
  <div style="font-family:'IBM Plex Mono',monospace;color:#58A6FF">{acc}</div>
  <div style="color:#8B949E;font-size:0.82rem;margin-top:6px">Will fetch from UniProt + AlphaFold</div>
</div>
""", unsafe_allow_html=True)

if run and accession_input.strip():
    acc = accession_input.strip().upper()
    st.markdown("---")
    log = st.empty()
    def status(msg, kind=''):
        log.markdown(f'<div class="status-line status-{kind}">{msg}</div>', unsafe_allow_html=True)

    status('⚙ Checking fpocket installation (first run takes ~2 min)...')
    fp_ok = install_fpocket()
    if not fp_ok:
        st.error('fpocket could not be installed.')
        st.stop()

    status('① Fetching metadata from UniProt...')
    meta = fetch_uniprot(acc)
    if acc in KNOWN:
        gene, pathway, desc = KNOWN[acc]
        meta['gene'] = gene

    status('② Downloading AlphaFold predicted structure...')
    pdb_bytes, err = fetch_alphafold(acc)
    if err:
        st.error(f'Could not fetch structure: {err}')
        st.stop()
    status('② Structure downloaded.', 'ok')

    status('③ Filtering low-confidence residues (pLDDT < 70)...')
    filtered_path, n_total, n_kept, plddt_list = filter_pdb(pdb_bytes, acc)
    if not filtered_path:
        st.error('pLDDT filtering failed.')
        st.stop()
    pct = int(100 * n_kept / n_total) if n_total else 0
    status(f'③ Filtered: {n_total} → {n_kept} residues kept ({pct}%)', 'ok')

    status('④ Running fpocket pocket detection...')
    out_dir = run_fpocket(filtered_path)
    pockets = parse_pockets(out_dir) if out_dir else []
    os.unlink(filtered_path)
    if not pockets:
        status('④ No pockets detected.', 'warn')
    else:
        status(f'④ {len(pockets)} pockets found. Top druggability: {pockets[0].get("drug_score","?"):.3f}', 'ok')

    log.empty()

    st.markdown(f"""
<div class="card">
  <div class="card-title">Protein overview</div>
  <div style="font-size:1.4rem;font-weight:600;color:#E6EDF3;font-family:'IBM Plex Mono',monospace;margin-bottom:4px">
    {meta.get('gene') or acc} &nbsp;<span style="font-size:0.9rem;color:#8B949E;font-weight:400">{acc}</span>
  </div>
  <div style="color:#8B949E;font-size:0.85rem;margin-bottom:10px">{meta.get('protein_name','')}</div>
  <div style="color:#8B949E;font-size:0.82rem;font-style:italic">{meta.get('organism','')}</div>
  <div class="metric-row" style="margin-top:14px">
    <div class="metric"><div class="metric-val">{meta.get('length',0)}</div><div class="metric-label">amino acids</div></div>
    <div class="metric"><div class="metric-val">{n_total}</div><div class="metric-label">total residues</div></div>
    <div class="metric"><div class="metric-val">{n_kept}</div><div class="metric-label">after pLDDT filter</div></div>
    <div class="metric"><div class="metric-val">{len(pockets)}</div><div class="metric-label">pockets found</div></div>
  </div>
</div>
""", unsafe_allow_html=True)

    st.markdown('<div class="card"><div class="card-title">pLDDT confidence profile</div>', unsafe_allow_html=True)
    st.pyplot(plddt_plot(plddt_list), use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)

    if pockets:
        top = pockets[0]
        ds = top.get('drug_score', 0)
        vol = top.get('volume', 0)
        hyd = top.get('hydrophobicity', 0)
        pol = top.get('polarity', 0)
        if ds >= 0.5:
            badge = '<span class="badge badge-good">Druggable</span>'
        elif ds >= 0.3:
            badge = '<span class="badge badge-warn">Moderately druggable</span>'
        else:
            badge = '<span class="badge badge-bad">Low druggability</span>'

        st.markdown(f"""
<div class="card">
  <div class="card-title">Top binding pocket (pocket #{top.get('pocket_num','?')})</div>
  {badge}
  <div class="metric-row" style="margin-top:12px">
    <div class="metric"><div class="metric-val">{ds:.3f}</div><div class="metric-label">druggability score</div></div>
    <div class="metric"><div class="metric-val">{vol:.0f}</div><div class="metric-label">volume (Å³)</div></div>
    <div class="metric"><div class="metric-val">{hyd:.1f}</div><div class="metric-label">hydrophobicity</div></div>
    <div class="metric"><div class="metric-val">{pol:.0f}</div><div class="metric-label">polarity score</div></div>
  </div>
</div>
""", unsafe_allow_html=True)

        st.markdown('<div class="card"><div class="card-title">All detected pockets</div>', unsafe_allow_html=True)
        fig2 = pocket_bar_plot(pockets)
        if fig2:
            st.pyplot(fig2, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)

        with st.expander('Full pocket table'):
            df = pd.DataFrame(pockets)
            df = df.rename(columns={
                'pocket_num': 'Pocket', 'drug_score': 'Drug score',
                'volume': 'Volume (Å³)', 'hydrophobicity': 'Hydrophobicity',
                'polarity': 'Polarity', 'n_spheres': 'Alpha spheres'
            })
            st.dataframe(df, use_container_width=True)
    else:
        st.warning('No pockets detected. The protein may be too small or too disordered after filtering.')

elif run:
    st.warning('Please enter a UniProt accession ID first.')
