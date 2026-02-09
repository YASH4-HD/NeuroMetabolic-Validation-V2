import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# --- PAGE CONFIG ---
st.set_page_config(page_title="NeuroMetabolic Validation", page_icon="🔬", layout="wide")

st.title("🔬 Clinical Validation & PPI Interactome")
st.markdown("Phase 3: STRING-DB Physical Interactions | Phase 4: GEO Patient Data Overlay")

# --- FUNCTIONS ---

@st.cache_data
def get_kegg_genes(pathway_id):
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    response = requests.get(url)
    genes = []
    if response.status_code == 200:
        lines = response.text.split('\n')
        is_gene_section = False
        for line in lines:
            if line.startswith('GENE'):
                is_gene_section = True
                line = line.replace('GENE', '').strip()
            elif line.startswith('COMPOUND') or line.startswith('REFERENCE'):
                is_gene_section = False
            if is_gene_section and line and ';' in line:
                parts = line.split('; ')
                id_sym = parts[0].split(None, 1)
                if len(id_sym) >= 2:
                    clean_symbol = id_sym[1].split(',')[0].strip()
                    genes.append({'Symbol': clean_symbol, 'Description': parts[1].strip()})
    return pd.DataFrame(genes)

def get_string_interactions(gene_list):
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": 9606, 
        "caller_identity": "research_app"
    }
    try:
        response = requests.post(url, data=params)
        return response.json()
    except:
        return []

# --- SIDEBAR ---
st.sidebar.header("🧬 Study Parameters")
pathway_map = {
    "Alzheimer's": "hsa05010",
    "Huntington's": "hsa05016", 
    "Parkinson's": "hsa05012",
    "Type II Diabetes": "hsa04930"
}
disease_choice = st.sidebar.selectbox("Target Pathology:", list(pathway_map.keys()))
pathway_id = pathway_map[disease_choice]

st.sidebar.header("📊 Phase 4: Clinical Data")
uploaded_file = st.sidebar.file_uploader("Upload GEO Patient CSV")

st.sidebar.header("⚙️ Network Settings")
confidence = st.sidebar.slider("STRING Confidence Score", 0, 1000, 400)
node_spread = st.sidebar.slider("Node Spacing", 1.0, 5.0, 2.0)

# --- DATA PROCESSING ---
df_kegg = get_kegg_genes(pathway_id)

if not df_kegg.empty:
    gene_list = df_kegg['Symbol'].unique().tolist()[:40]
    
    with st.spinner('Building Interactome...'):
        interactions = get_string_interactions(gene_list)
    
    if uploaded_file:
        geo_df = pd.read_csv(uploaded_file)
        df_kegg = pd.merge(df_kegg, geo_df[['Symbol', 'LogFC']], on='Symbol', how='left')
        df_kegg['LogFC'] = df_kegg['LogFC'].fillna(0)
    else:
        df_kegg['LogFC'] = 0

    # --- NETWORK ---
    G = nx.Graph()
    for _, row in df_kegg.iterrows():
        if row['Symbol'] in gene_list:
            G.add_node(row['Symbol'], logfc=row['LogFC'])

    edges_found = 0
    if isinstance(interactions, list):
        for edge in interactions:
            if edge['score'] >= (confidence / 1000):
                n_a, n_b = edge['preferredName_A'], edge['preferredName_B']
                if n_a in G.nodes() and n_b in G.nodes():
                    G.add_edge(n_a, n_b)
                    edges_found += 1

    # --- VISUALS ---
    col1, col2 = st.columns([3, 1])

    with col1:
        fig, ax = plt.subplots(figsize=(12, 10))
        pos = nx.spring_layout(G, k=node_spread/np.sqrt(len(G.nodes())), iterations=50)
        
        node_colors = []
        for node in G.nodes():
            val = df_kegg.loc[df_kegg['Symbol'] == node, 'LogFC'].values
            val = val[0] if len(val) > 0 else 0
            if val > 1.0: node_colors.append('#FF4B4B')
            elif val < -1.0: node_colors.append('#4B4BFF')
            else: node_colors.append('#D5D8DC')

        nx.draw_networkx_edges(G, pos, alpha=0.3)
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1000)
        nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
        
        plt.axis('off')
        st.pyplot(fig)

    with col2:
        st.metric("PPI Edges", edges_found)
        st.metric("Nodes", G.number_of_nodes())
        st.write("**Top Hubs**")
        degrees = dict(G.degree())
        for hub, deg in sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:5]:
            st.write(f"• {hub}: {deg} links")
else:
    st.error("Check connection.")
