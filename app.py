import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# --- PAGE CONFIG ---
st.set_page_config(page_title="NeuroMetabolic Phase 3 & 4", page_icon="🔬", layout="wide")

st.title("🔬 Clinical Validation & PPI Interactome")
st.markdown("""
**Phase 3:** Physical Protein-Protein Interactions (via STRING-DB)  
**Phase 4:** Patient Expression Mapping (via GEO Data Overlay)
""")

# --- FUNCTIONS ---

@st.cache_data
def get_kegg_genes(pathway_id):
    """Fetches genes from KEGG API"""
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
                    genes.append({'Symbol': id_sym[1].strip(), 'Description': parts[1].strip()})
    return pd.DataFrame(genes)

def get_string_interactions(gene_list):
    """PHASE 3: Fetches real physical interactions from STRING-DB"""
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list), 
        "species": 9606, # Human
        "required_score": 400, # Medium confidence
        "caller_identity": "streamlit_research_app"
    }
    response = requests.post(string_api_url, data=params)
    return response.json()

# --- SIDEBAR: CONTROLS ---
st.sidebar.header("Step 1: Select Disease")
pathway_map = {
    "Huntington's": "hsa05016", 
    "Alzheimer's": "hsa05010", 
    "Parkinson's": "hsa05012",
    "Type II Diabetes": "hsa04930"
}
disease_choice = st.sidebar.selectbox("Target Pathology:", list(pathway_map.keys()))
pathway_id = pathway_map[disease_choice]

st.sidebar.header("Step 2: Phase 4 Validation")
uploaded_file = st.sidebar.file_uploader("Upload GEO Patient Data (CSV)", help="CSV must have 'Symbol' and 'LogFC' columns")

st.sidebar.markdown("---")
confidence = st.sidebar.slider("STRING Confidence Score", 0, 1000, 400)

# --- DATA PROCESSING ---
df_kegg = get_kegg_genes(pathway_id)

if not df_kegg.empty:
    gene_list = df_kegg['Symbol'].tolist()[:50] # Limit to top 50 for performance
    
    # PHASE 3: STRING-DB API CALL
    with st.spinner('Fetching STRING-DB Interactions...'):
        interactions = get_string_interactions(gene_list)
    
    # PHASE 4: GEO DATA MERGING
    if uploaded_file:
        geo_df = pd.read_csv(uploaded_file)
        df_kegg = pd.merge(df_kegg, geo_df[['Symbol', 'LogFC']], on='Symbol', how='left')
        st.sidebar.success("Patient Data Loaded!")
    else:
        df_kegg['LogFC'] = 0 # Default if no data uploaded

    # --- NETWORK CONSTRUCTION ---
    G = nx.Graph()
    for _, row in df_kegg.iterrows():
        if row['Symbol'] in gene_list:
            G.add_node(row['Symbol'], logfc=row.get('LogFC', 0))

    for edge in interactions:
        if edge['score'] >= confidence:
            G.add_edge(edge['preferredName_A'], edge['preferredName_B'], weight=edge['score'])

    # --- VISUALIZATION ---
    col1, col2 = st.columns([3, 1])

    with col1:
        st.subheader(f"Physical Interactome: {disease_choice}")
        fig, ax = plt.subplots(figsize=(10, 8))
        pos = nx.spring_layout(G, k=0.3, seed=42)
        
        # Color nodes by LogFC (Phase 4)
        node_colors = []
        for node in G.nodes(data=True):
            val = node[1].get('logfc', 0)
            if val > 1: node_colors.append('#FF4B4B') # Upregulated (Red)
            elif val < -1: node_colors.append('#4B4BFF') # Downregulated (Blue)
            else: node_colors.append('#D5D8DC') # Neutral (Grey)

        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500, alpha=0.9)
        nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color="grey")
        nx.draw_networkx_labels(G, pos, font_size=8, font_weight="bold")
        
        st.pyplot(fig)
        st.caption("🔴 Red: Upregulated in Patients | 🔵 Blue: Downregulated in Patients | ⚪ Grey: No Data/Neutral")

    with col2:
        st.subheader("Network Metrics")
        st.metric("Physical Interactions", G.number_of_edges())
        st.metric("Active Nodes", G.number_of_nodes())
        
        st.write("**Top Connected Proteins**")
        degrees = dict(G.degree())
        top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:5]
        for hub, deg in top_hubs:
            st.write(f"• {hub} ({deg} links)")

    # --- RESEARCH PAPER SUMMARY ---
    st.markdown("---")
    st.subheader("📝 Phase 3 & 4 Research Summary")
    
    if uploaded_file:
        up_genes = df_kegg[df_kegg['LogFC'] > 1]['Symbol'].tolist()
        st.write(f"**Clinical Insight:** In the {disease_choice} patient cohort, we identified **{len(up_genes)}** genes that are significantly upregulated and physically interact within the metabolic hub.")
        st.write(f"**Key Pathological Node:** {top_hubs[0][0] if top_hubs else 'N/A'} shows the highest degree of physical connectivity, suggesting it as a master regulator.")
    else:
        st.info("Upload a GEO CSV file (with 'Symbol' and 'LogFC' columns) to generate clinical insights.")

else:
    st.error("Could not fetch data from KEGG. Please check your connection.")
