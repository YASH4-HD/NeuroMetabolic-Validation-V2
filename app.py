import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

# --- PAGE CONFIG ---
st.set_page_config(page_title="NeuroMetabolic Phase 3 & 4", page_icon="🔬", layout="wide")

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
                    # Clean the symbol (remove commas or extra spaces)
                    clean_symbol = id_sym[1].split(',')[0].strip()
                    genes.append({'Symbol': clean_symbol, 'Description': parts[1].strip()})
    return pd.DataFrame(genes)

def get_string_interactions(gene_list):
    """PHASE 3: Fetches interactions and maps names correctly"""
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": 9606, # Human
        "caller_identity": "research_app"
    }
    response = requests.post(url, data=params)
    return response.json()

# --- SIDEBAR ---
st.sidebar.header("Step 1: Select Disease")
pathway_map = {
    "Alzheimer's": "hsa05010",
    "Huntington's": "hsa05016", 
    "Parkinson's": "hsa05012",
    "Type II Diabetes": "hsa04930"
}
disease_choice = st.sidebar.selectbox("Target Pathology:", list(pathway_map.keys()))
pathway_id = pathway_map[disease_choice]

st.sidebar.header("Step 2: Phase 4 Validation")
uploaded_file = st.sidebar.file_uploader("Upload GEO Patient Data (CSV)")
confidence = st.sidebar.slider("STRING Confidence Score", 0, 1000, 400)

# --- DATA PROCESSING ---
df_kegg = get_kegg_genes(pathway_id)

if not df_kegg.empty:
    # Use the first 40 genes for a cleaner network
    gene_list = df_kegg['Symbol'].unique().tolist()[:40]
    
    with st.spinner('Connecting to STRING-DB...'):
        interactions = get_string_interactions(gene_list)
    
    # Merge GEO Data
    if uploaded_file:
        geo_df = pd.read_csv(uploaded_file)
        df_kegg = pd.merge(df_kegg, geo_df[['Symbol', 'LogFC']], on='Symbol', how='left')
    else:
        df_kegg['LogFC'] = 0

    # --- NETWORK CONSTRUCTION ---
    G = nx.Graph()
    
    # Add Nodes
    for _, row in df_kegg.iterrows():
        if row['Symbol'] in gene_list:
            G.add_node(row['Symbol'], logfc=row.get('LogFC', 0))

    # Add Edges from STRING
    edges_found = 0
    if isinstance(interactions, list):
        for edge in interactions:
            if edge['score'] >= (confidence / 1000): # STRING API uses 0-1 scale
                # STRING often returns slightly different names, we map them back
                node_a = edge['preferredName_A']
                node_b = edge['preferredName_B']
                if node_a in gene_list and node_b in gene_list:
                    G.add_edge(node_a, node_b)
                    edges_found += 1

    # --- VISUALIZATION ---
    col1, col2 = st.columns([3, 1])

    with col1:
        st.subheader(f"Physical Interactome: {disease_choice}")
        if edges_found == 0:
            st.warning("No physical interactions found at this confidence level. Try lowering the Confidence Score.")
        
        fig, ax = plt.subplots(figsize=(12, 10))
        # Use spring_layout to show clusters
        pos = nx.spring_layout(G, k=0.4, iterations=50)
        
        # Color nodes by LogFC
        node_colors = []
        for node in G.nodes():
            val = df_kegg.loc[df_kegg['Symbol'] == node, 'LogFC'].values
            val = val[0] if len(val) > 0 else 0
            if val > 1: node_colors.append('#FF4B4B') # Red
            elif val < -1: node_colors.append('#4B4BFF') # Blue
            else: node_colors.append('#D5D8DC') # Grey

        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=800, alpha=0.9)
        nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color="grey", width=2)
        nx.draw_networkx_labels(G, pos, font_size=9, font_weight="bold")
        
        plt.axis('off')
        st.pyplot(fig)

    with col2:
        st.subheader("Metrics")
        st.metric("Physical Interactions", edges_found)
        st.metric("Active Nodes", G.number_of_nodes())
        
        st.write("**Central Hubs**")
        degrees = dict(G.degree())
        for hub, deg in sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:5]:
            if deg > 0:
                st.write(f"• **{hub}**: {deg} connections")

else:
    st.error("Data fetch failed.")
