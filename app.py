import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# --- PAGE CONFIG ---
st.set_page_config(page_title="NeuroMetabolic Validation v2.5", page_icon="🔬", layout="wide")

st.title("🔬 Clinical Validation & PPI Interactome")
st.markdown("### Systems Biology Pipeline: Phase 3 (Interactome) & Phase 4 (Clinical Overlay)")

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
                    clean_symbol = id_sym[1].split(',')[0].strip().upper()
                    genes.append({'Symbol': clean_symbol, 'Description': parts[1].strip()})
    return pd.DataFrame(genes)

def get_string_interactions(gene_list):
    url = "https://string-db.org/api/json/network"
    params = {"identifiers": "%0d".join(gene_list), "species": 9606, "caller_identity": "research_app"}
    try:
        response = requests.post(url, data=params)
        data = response.json()
        return data if isinstance(data, list) else []
    except:
        return []

# --- SIDEBAR ---
st.sidebar.header("🧬 Study Parameters")
pathway_map = {"Alzheimer's": "hsa05010", "Huntington's": "hsa05016", "Parkinson's": "hsa05012", "Type II Diabetes": "hsa04930"}
disease_choice = st.sidebar.selectbox("Target Pathology:", list(pathway_map.keys()))
pathway_id = pathway_map[disease_choice]

st.sidebar.header("📊 Phase 4: Clinical Validation")
uploaded_file = st.sidebar.file_uploader("Upload Patient Data (GEO CSV)")

st.sidebar.header("⚙️ Network Rigor")
confidence = st.sidebar.slider("STRING Interaction Confidence Threshold", 0, 1000, 400)
node_spread = st.sidebar.slider("Node Spacing (Layout Force)", 1.0, 5.0, 3.4)

st.sidebar.header("🎨 Visualization Polish")
# REFINEMENT C: Defaulting to "Hubs Only" for better UX/readability
label_mode = st.sidebar.radio("Show Labels for:", ["Hubs Only (Degree > 2)", "All Nodes", "None"], index=0)

# REFINEMENT B: Global Disclaimer Line
st.sidebar.info("💡 *Network topology reflects functional coupling and pathway co-occurrence, not direct molecular causality.*")

# --- DATA PROCESSING ---
df_kegg = get_kegg_genes(pathway_id)

if not df_kegg.empty:
    gene_list = df_kegg['Symbol'].unique().tolist()[:45]
    
    with st.spinner('Calculating Interactome Topology...'):
        interactions = get_string_interactions(gene_list)
    
    if uploaded_file:
        try:
            geo_df = pd.read_csv(uploaded_file)
            geo_df['Symbol'] = geo_df['Symbol'].astype(str).str.strip().str.upper()
            geo_df['LogFC'] = pd.to_numeric(geo_df['LogFC'], errors='coerce').fillna(0)
            df_kegg = pd.merge(df_kegg, geo_df[['Symbol', 'LogFC']], on='Symbol', how='left')
            df_kegg['LogFC'] = df_kegg['LogFC'].fillna(0)
        except Exception as e:
            st.error(f"CSV Error: {e}")
            df_kegg['LogFC'] = 0
    else:
        df_kegg['LogFC'] = 0

    # --- NETWORK CONSTRUCTION ---
    G = nx.Graph()
    for _, row in df_kegg.iterrows():
        if row['Symbol'] in gene_list:
            G.add_node(row['Symbol'], logfc=row['LogFC'])

    edges_found = 0
    if interactions:
        for edge in interactions:
            if 'preferredName_A' in edge and 'preferredName_B' in edge:
                if edge['score'] >= (confidence / 1000):
                    n_a, n_b = edge['preferredName_A'].upper(), edge['preferredName_B'].upper()
                    if n_a in G.nodes() and n_b in G.nodes():
                        G.add_edge(n_a, n_b)
                        edges_found += 1

    # --- VISUALIZATION ---
    col1, col2 = st.columns([3, 1])

    with col1:
        fig, ax = plt.subplots(figsize=(12, 10))
        pos = nx.spring_layout(G, k=(node_spread)/np.sqrt(len(G.nodes())), iterations=50)
        
        node_colors = []
        for node in G.nodes():
            val = df_kegg.loc[df_kegg['Symbol'] == node, 'LogFC'].max()
            if val > 1.0: node_colors.append('#FF4B4B')
            elif val < -1.0: node_colors.append('#4B4BFF')
            else: node_colors.append('#D5D8DC')

        nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='grey')
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1000, edgecolors='white')
        
        # UX Logic for Labels
        labels = {n: n for n in G.nodes() if (label_mode == "All Nodes" or (label_mode == "Hubs Only (Degree > 2)" and G.degree(n) > 2))}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold')
        
        plt.axis('off')
        st.pyplot(fig)
        
        focus_area = "Insulin signaling" if "Diabetes" in disease_choice else "Mitochondrial ETC pathways"
        
        # REFINEMENT A: Updated Interpretation and Legend Wording
        st.markdown(f"**Analysis Interpretation:** *Highlighted hubs represent high-connectivity genes emerging under STRING confidence ≥ {confidence/1000}; when expression data is available, nodes are additionally colored by differential regulation, suggesting **{focus_area}** as a key driver of pathology in {disease_choice}.*")
        st.markdown(f"""
        **Node Legend:**
        - 🔴 = **High-centrality hubs** (upregulated when expression data is available)
        - 🔵 = **Downregulated** (LogFC < -1.0)
        - ⚪ = **Background / No expression data**
        """)

    with col2:
        st.subheader("Network Metrics")
        st.metric("Validated Edges", edges_found)
        st.metric("Total Nodes", G.number_of_nodes())
        st.write("---")
        st.write("**Top Centrality Hubs**")
        degrees = dict(G.degree())
        for hub, deg in sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:8]:
            if deg > 0: st.write(f"• **{hub}**: {deg} interactions")
else:
    st.error("Data fetch failed.")
