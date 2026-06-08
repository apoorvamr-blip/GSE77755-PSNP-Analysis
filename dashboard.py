#!/usr/bin/env python3
"""
GSE77755 Zebrafish PSNP Immune Response Dashboard
Streamlit-ready dashboard for hackathon presentation
Run: streamlit run dashboard.py
"""

import streamlit as st
import pandas as pd
import plotly.express as px

st.set_page_config(
    page_title="PSNP Immune Response | GSE77755",
    page_icon="🐟",
    layout="wide"
)

# ── CSS Styling ──
st.markdown("""
<style>
.big-font { font-size:20px !important; font-weight: bold; }
.metric-card { background-color: #1e3a5f; padding: 15px; border-radius: 10px; color: white; }
.highlight { background-color: #ff4b4b20; border-left: 4px solid #ff4b4b; padding: 10px; }
</style>
""", unsafe_allow_html=True)

# ── Header ──
st.header("Immune-Related Differentially Expressed Genes")
st.markdown("**Dataset:** GSE77755 | **Organism:** *Danio rerio* | **Analysis:** DESeq2 RNA-seq")
st.markdown("---")

# ── DATA ──
# DEGs 1dpi
deg1 = pd.read_csv("data/DEG_1dpi.csv")

# DEGs 3dpi
deg3 = pd.read_csv("data/DEG_3dpi.csv")

deg1 = deg1.rename(columns={
    "X": "Gene",
    "log2FoldChange": "log2FC"
})

deg3 = deg3.rename(columns={
    "X": "Gene",
    "log2FoldChange": "log2FC"
})

# Add timepoint labels
deg1["Timepoint"] = "1 dpi"
deg3["Timepoint"] = "3 dpi"

gene_map = pd.read_csv("data/Final_Immune_Gene_Table.csv")

symbol_dict = dict(
    zip(
        gene_map["SYMBOL"],
        gene_map["GENENAME"]
    )
)

immune_top20 = pd.read_csv(
    "data/Immune_Genes_GSE77755.csv"
)
# Combined DEG table
all_degs = pd.concat([deg1, deg3], ignore_index=True)
# ── TABS ──
tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
    "Overview",
    "DEGs",
    "Immune Genes",
    "Functional Interpretation",
    "Proposed Biological Mechanism",
    "Summary & Discussion",
    "Figures"
])

with tab1:
    st.header("Study Overview")

    c1, c2, c3, c4 = st.columns(4)

    c1.metric("Total Genes Tested", "~31,000")
    c2.metric("DEGs at 1 dpi", "4")
    c3.metric("DEGs at 3 dpi", "4")
    c4.metric("Immune Genes Identified", "2")

    st.markdown("---")

    st.subheader("Temporal Response to Polystyrene Nanoparticle Exposure")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        ###  1 dpi — Early Stress Response

        - **npas4a** was the most significantly upregulated gene
        - Additional zebrafish-specific genes showed altered expression
        - No strong immune-related genes were detected
        - Indicates an early cellular stress/adaptation response

        **Biological theme:** Initial stress signaling
        """)

    with col2:
        st.markdown("""
        ###  3 dpi — Innate Immune Activation

        - **c3a.3** significantly upregulated
        - **cfb** significantly upregulated
        - **cetp** associated with lipid metabolism changes
        - **igfbp1a** linked to stress and growth regulation

        **Biological theme:** Complement-mediated innate immunity
        """)

    st.markdown("---")

    st.subheader("Log2 Fold Change of Significant Differentially Expressed Genes")

    fig = px.bar(
        all_degs,
        x="Gene",
        y="log2FC",
        color="Timepoint",
        facet_col="Timepoint",
        text="log2FC",
        title="Differentially Expressed Genes Across Timepoints",
        height=500
    )

    fig.update_traces(
        texttemplate="%{text:.2f}",
        textposition="outside"
    )

    st.plotly_chart(fig, use_container_width=True)

    st.info(
        "The transcriptomic response was limited at 1 dpi but shifted toward complement-associated immune activation by 3 dpi, particularly through c3a.3 and cfb."
    )
    
# ──────────────────────────────────────────
# TAB 2: DEGs
# ──────────────────────────────────────────
with tab2:

    st.header("Differentially Expressed Genes")

    col1, col2 = st.columns(2)

    with col1:

        st.subheader("1 dpi Significant Genes")

        st.dataframe(
            deg1[["Gene","log2FC","padj"]]
            .sort_values("padj"),
            use_container_width=True
        )

        fig1 = px.bar(
            deg1,
            x="Gene",
            y="log2FC",
            text="log2FC",
            title="1 dpi Differentially Expressed Genes"
        )

        fig1.update_traces(
            texttemplate="%{text:.2f}",
            textposition="outside"
        )

        st.plotly_chart(fig1, use_container_width=True)

    with col2:

        st.subheader("3 dpi Significant Genes")

        st.dataframe(
            deg3[["Gene","log2FC","padj"]]
            .sort_values("padj"),
            use_container_width=True
        )

        fig3 = px.bar(
            deg3,
            x="Gene",
            y="log2FC",
            text="log2FC",
            title="3 dpi Differentially Expressed Genes"
        )

        fig3.update_traces(
            texttemplate="%{text:.2f}",
            textposition="outside"
        )

        st.plotly_chart(fig3, use_container_width=True)

    st.markdown("---")

    st.subheader("Combined DEG Table")

    st.dataframe(
        all_degs[
            ["Gene","Timepoint","log2FC","padj"]
        ].sort_values("padj"),
        use_container_width=True
    )


# ──────────────────────────────────────────
# TAB 3: Immune Genes
# ──────────────────────────────────────────
with tab3:

    st.header("Immune-Related Differentially Expressed Genes")

    st.dataframe(
        immune_top20,
        use_container_width=True
    )

    fig = px.bar(
        immune_top20,
        x="Gene_Symbol",
        y="log2FC",
        color="Gene_Symbol",
        text="log2FC",
        title="Immune Gene Expression Changes"
    )

    fig.update_traces(
        texttemplate="%{text:.2f}",
        textposition="outside"
    )

    st.plotly_chart(fig, use_container_width=True)

    st.info(
        "Two immune-associated genes were significantly upregulated at 3 dpi: c3a.3 (Complement component C3a) and cfb (Complement factor B), suggesting activation of complement-mediated innate immunity."
    )

# ──────────────────────────────────────────
# TAB 4: Functional Interpretation
# ──────────────────────────────────────────
with tab4:

    st.header("Functional Interpretation of Significant Genes")

    interpretation = pd.DataFrame([
        {
            "Gene":"npas4a",
            "Timepoint":"1 dpi",
            "Function":"Neuronal stress response and transcriptional regulation",
            "Biological Significance":"Early cellular stress response following PSNP exposure"
        },
        {
            "Gene":"cetp",
            "Timepoint":"3 dpi",
            "Function":"Cholesterol and lipid transport",
            "Biological Significance":"Suggests metabolic disruption following exposure"
        },
        {
            "Gene":"c3a.3",
            "Timepoint":"3 dpi",
            "Function":"Complement component C3a",
            "Biological Significance":"Activation of innate immune signaling"
        },
        {
            "Gene":"cfb",
            "Timepoint":"3 dpi",
            "Function":"Complement factor B",
            "Biological Significance":"Activation of alternative complement pathway"
        },
        {
            "Gene":"igfbp1a",
            "Timepoint":"3 dpi",
            "Function":"Growth and stress regulation",
            "Biological Significance":"Adaptive response to environmental stress"
        }
    ])

    st.dataframe(
        interpretation,
        use_container_width=True
    )

    st.markdown("---")

    st.subheader("Major Biological Findings")

    st.success("""
    1. A limited transcriptional response was observed at 1 dpi.

    2. npas4a was the strongest early-response gene, suggesting
       activation of cellular stress pathways.

    3. By 3 dpi, the response shifted toward innate immunity.

    4. Significant upregulation of c3a.3 and cfb indicates
       complement-mediated immune activation.

    5. cetp and igfbp1a suggest additional metabolic and
       stress-related adaptations.
    """)

    st.markdown("---")

    st.subheader("Immune Genes Identified")

    st.markdown("""
    **c3a.3** and **cfb** were the only clearly immune-related
    genes that passed the differential expression threshold.

    Together these genes support activation of the complement
    cascade, one of the primary components of the innate immune
    system in zebrafish.
    """)

# ──────────────────────────────────────────
# TAB 5: Proposed Biological Mechanism
# ──────────────────────────────────────────
with tab5:

    st.header("Proposed Biological Mechanism")

    st.markdown("""
    ### Polystyrene Nanoparticle Exposure Response

    PSNP Exposure
    ↓
    Cellular Stress Response
    ↓
    npas4a Activation (1 dpi)
    ↓
    Continued Exposure
    ↓
    Innate Immune Activation
    ↓
    c3a.3 + cfb Upregulation (3 dpi)
    ↓
    Complement Cascade Activation
    ↓
    Inflammatory Response
    """)

    st.markdown("---")

    st.subheader("Mechanistic Interpretation")

    mechanism_df = pd.DataFrame([
        {
            "Stage":"Early Response",
            "Gene":"npas4a",
            "Role":"Stress-responsive transcription factor"
        },
        {
            "Stage":"Immune Activation",
            "Gene":"c3a.3",
            "Role":"Complement-mediated immune signaling"
        },
        {
            "Stage":"Complement Amplification",
            "Gene":"cfb",
            "Role":"Alternative complement pathway activation"
        },
        {
            "Stage":"Metabolic Response",
            "Gene":"cetp",
            "Role":"Lipid metabolism and transport"
        },
        {
            "Stage":"Adaptive Response",
            "Gene":"igfbp1a",
            "Role":"Stress and growth regulation"
        }
    ])

    st.dataframe(
        mechanism_df,
        use_container_width=True
    )

    st.markdown("---")

    st.success("""
    Proposed mechanism:

    1. Polystyrene nanoparticles trigger an initial cellular stress response.

    2. Early transcriptional changes are dominated by npas4a.

    3. Continued exposure activates innate immune pathways.

    4. Significant upregulation of c3a.3 and cfb indicates
       activation of the complement system.

    5. cetp and igfbp1a suggest additional metabolic and
       physiological adaptation to nanoparticle exposure.
    """)

# ──────────────────────────────────────────
# TAB 6: Summary
# ──────────────────────────────────────────
with tab6:

    st.header("Summary and Discussion")

    st.subheader("Key Findings")

    findings = [
        (
            "1",
            "Limited Early Response",
            "Only four genes were significantly differentially expressed at 1 dpi, indicating a modest early transcriptional response to PSNP exposure."
        ),
        (
            "2",
            "npas4a as the Major Early Response Gene",
            "npas4a showed the strongest statistical signal at 1 dpi and is associated with cellular stress and transcriptional regulation."
        ),
        (
            "3",
            "Complement Activation at 3 dpi",
            "c3a.3 and cfb were significantly upregulated at 3 dpi, suggesting activation of complement-mediated innate immunity."
        ),
        (
            "4",
            "Metabolic and Physiological Adaptation",
            "cetp and igfbp1a indicate additional metabolic and stress-related responses following prolonged exposure."
        ),
        (
            "5",
            "Evidence for Delayed Immune Response",
            "The transcriptomic profile shifted from stress-associated signaling at 1 dpi to immune-associated pathways at 3 dpi."
        )
    ]

    for num, title, detail in findings:
        with st.expander(f"Finding {num}: {title}"):
            st.write(detail)

    st.markdown("---")

    st.subheader("Executive Summary")

    st.markdown("""
    Polystyrene nanoparticle (PSNP) exposure was investigated using the
    zebrafish RNA-seq dataset GSE77755.

    Differential expression analysis identified four significant genes
    at 1 dpi and four significant genes at 3 dpi.

    The strongest early response gene was **npas4a**, suggesting
    activation of cellular stress pathways shortly after exposure.

    By 3 dpi, significant upregulation of **c3a.3** and **cfb**
    indicated activation of complement-mediated innate immunity.
    Additional genes (**cetp** and **igfbp1a**) suggested metabolic
    and physiological adaptation to nanoparticle exposure.

    Overall, the results support a model in which PSNP exposure
    triggers an early stress response followed by a delayed innate
    immune response involving complement pathway activation.
    """)

    st.markdown("---")

    st.subheader("Discussion Points")

    discussion = [
        (
            "Why were only eight genes significant?",
            "Stringent multiple-testing correction and biological variability reduced the number of statistically significant genes."
        ),
        (
            "What is the strongest immune signal?",
            "The complement-associated genes c3a.3 and cfb showed significant upregulation at 3 dpi."
        ),
        (
            "What does npas4a suggest?",
            "npas4a indicates an early stress-response program activated shortly after exposure."
        ),
        (
            "Why are cetp and igfbp1a important?",
            "These genes suggest that PSNP exposure may affect metabolic and physiological regulation in addition to immunity."
        ),
        (
            "What is the overall biological interpretation?",
            "PSNP exposure appears to induce an initial stress response followed by complement-mediated innate immune activation."
        )
    ]

    for q, a in discussion:
        with st.expander(q):
            st.info(a)

    st.markdown("---")

    st.success("""
    Final Conclusion:

    Polystyrene nanoparticle exposure produced a limited but biologically
    meaningful transcriptional response in zebrafish larvae. The most
    notable finding was the upregulation of complement-associated genes
    c3a.3 and cfb at 3 dpi, indicating activation of innate immune
    pathways following sustained exposure.
    """)

    st.caption(
        "Dataset: GSE77755 | Analysis: DESeq2 | Organism: Danio rerio | Dashboard: Streamlit"
    )

# ──────────────────────────────────────────
# TAB 7: Figures
# ──────────────────────────────────────────
with tab7:

    st.header("Analysis Figures")

    st.subheader("PCA Plot")
    st.image("figures/PCA_plot.png")

    st.subheader("Volcano Plot (1 dpi)")
    st.image("figures/Volcano_1dpi.png")

    st.subheader("Volcano Plot (3 dpi)")
    st.image("figures/Volcano_3dpi.png")

    st.subheader("Heatmap")
    st.image("figures/Heatmap_DEGs.png")
