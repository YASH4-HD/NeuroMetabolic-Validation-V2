import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# --- PAGE CONFIG ---
st.set_page_config(page_title="NeuroMetabolic Phase 3 & 4", page_icon="🔬", layout="wide")

# Custom CSS to improve look
st.markdown("""
    <style>
    .main { background-color: #f5f7f9; }
    .stMetric { background-color: #ffffff; padding: 15px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_allow_stdio=True)

st.title("🔬 Clinical Validation & PPI Interactome")
st.markdown("### Phase 3: STRING-DB Physical Interactions | Phase 4: GEO Patient Data Overlay")

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
    "Type II Diabetes": "hsa04930",
    "ALS": "hsa05014"
}
disease_choice = st.sidebar.selectbox("Target Pathology:", list(pathway_map.keys()))
pathway_id = pathway_map[disease_choice]

st.sidebar.header("📊 Phase 4: Clinical Data")
uploaded_file = st.sidebar.file_uploader("Upload GEO Patient CSV", help="Requires 'Symbol' and 'LogFC' columns")

st.sidebar.header("⚙️ Network Settings")
confidence = st.sidebar.slider("STRING Confidence Score", 0, 1000, 400)
node_spread = st.sidebar.slider("Node Spacing (Layout)", 0.5, 5.0, 2.5)

# --- DATA PROCESSING ---
df_kegg = get_kegg_genes(pathway_id)

if not df_kegg.empty:
    # Use top 50 genes
    gene_list = df_kegg['Symbol'].unique().tolist()[:50]
    
    with st.spinner('Building Physical Interactome...'):
        interactions = get_string_interactions(gene_list)
    
    # Merge GEO Data
    if uploaded_file:
        geo_df = pd.read_csv(uploaded_file)
        df_kegg = pd.merge(df_kegg, geo_df[['Symbol', 'LogFC']], on='Symbol', how='left')
        df_kegg['LogFC'] = df_kegg['LogFC'].fillna(0)
    else:
        df_kegg['LogFC'] = 0

    # --- NETWORK CONSTRUCTION ---
    G = nx.Graph()
    for _, row in df_kegg.iterrows():
        if row['Symbol'] in gene_list:
            G.add_node(row['Symbol'], logfc=row['LogFC'])

    edges_found = 0
    if isinstance(interactions, list):
        for edge in interactions:
            if edge['score'] >= (confidence / 1000):
                node_a, node_b = edge['preferredName_A'], edge['preferredName_B']
                if node_a in G.nodes() and node_b in G.nodes():
                    G.add_edge(node_a, node_b)
                    edges_found += 1

    # --- VISUALIZATION ---
    col1, col2 = st.columns([3, 1])

    with col1:
        st.subheader(f"Validated Interactome: {disease_choice}")
        
        # IMPROVED LAYOUT FOR READABILITY
        fig, ax = plt.subplots(figsize=(14, 11))
        pos = nx.spring_layout(G, k=node_spread/np.sqrt(len(G.nodes())), iterations=100)
        
        # Color nodes by LogFC
        node_colors = []
        for node in G.nodes():
            val = df_kegg.loc[df_kegg['Symbol'] == node, 'LogFC'].values
            val = val[0] if len(val) > 0 else 0
            if val > 1.0: node_colors.append('#ef5350') # Significant Up
            elif val < -1.0: node_colors.append('#42a5f5') # Significant Down
            else: node_colors.append('#cfd8dc') # Neutral

        # Draw elements
        nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color="#90a4ae", width=1.5)
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1200, 
                               edgecolors="white", linewidths=2)
        
        # Optimized Labels (White background for readability)
        for node, (x, y) in pos.items():
            plt.text(x, y, s=node, fontsize=9, fontweight="bold", 
                     ha='center', va='center',
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
        
        plt.axis('off')
        st.pyplot(fig)
        st.info("💡 **Visual Guide:** Red nodes are upregulated in patients. Blue nodes are downregulated. Lines indicate physical protein-protein binding.")

    with col2:
        st.subheader("Network Analytics")
        st.metric("PPI Edges", edges_found)
        st.metric("Clinical Nodes", G.number_of_nodes())
        
        st.write("---")
        st.write("**Top Pathological Hubs**")
        degrees = dict(G.degree())
        sorted_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:8]
        
        hub_data = []
        for hub, deg in sorted_hubs:
            lfc = df_kegg.loc[df_kegg['Symbol'] == hub, 'LogFC'].values[0]
            status = "🔴" if lfc > 1 else ("🔵" if lfc < -1 else "⚪")
            st.write(f"{status} **{hub}**: {deg} links")
            hub_data.append({"Gene": hub, "Links": deg, "LogFC": lfc})

        st.write("---")
        # Download Results
        csv = pd.DataFrame(hub_data).to_csv(index=False)
        st.download_button("📩 Download Hub Analysis", csv, "hubs.csv", "text/csv")

else:
    st.error("Failed to fetch biological data. Please refresh.")
