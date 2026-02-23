import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import base64
from utils.parsing import parse_stamp_file, extract_tissue_name
from components.downloads import create_download_button
from components.styling import apply_plot_style

def show():
    """Upload and Single Tissue Analysis Page"""
    
    st.header("📤 Upload & Single Tissue Analysis")
    st.markdown("Upload STAMP gene switching files for individual tissue analysis.")
    
    # Global age groups
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]
    
    # Upload section
    st.markdown("""
    <div class="analysis-section">
        <h3>📂 File Upload</h3>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    with col1:
        file1 = st.file_uploader(
            "📂 Upload the first file", 
            type=["txt"], 
            key="file1",
            help="Upload a .txt file with gene switching data"
        )
    with col2:
        file2 = st.file_uploader(
            "📂 Upload the second file (optional)", 
            type=["txt"], 
            key="file2",
            help="Optional: Upload a second file for comparison"
        )
    
    # Age group selection
    st.markdown("### 🎯 Age Group Selection")
    selected_ages = st.multiselect(
        "Select age groups to display:",
        age_groups, 
        default=age_groups,
        help="Choose which age groups to include in the analysis"
    )
    
    if not selected_ages:
        st.warning("⚠️ Please select at least one age group.")
        return
    
    # Parse files
    fasce1, counts1, df1 = parse_stamp_file(file1, age_groups)
    fasce2, counts2, df2 = parse_stamp_file(file2, age_groups)
    
    # === SINGLE TISSUE ANALYSIS ===
    if file1 and not file2:
        st.markdown("""
        <div class="analysis-section">
            <h2>📊 Single Tissue Analysis Results</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # NOME PULITO USANDO LA FUNZIONE
        tissue_name = extract_tissue_name(file1.name)
        
        # Filter data
        df1_filtered = df1[df1["Age"].isin(selected_ages)]
        counts1_filtered = [
            df1_filtered[df1_filtered["Age"] == age].shape[0]
            for age in age_groups if age in selected_ages
        ]
        
        # Display metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(df1)}</h3>
                <p>Total Genes</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(selected_ages)}</h3>
                <p>Age Groups</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            avg_genes = np.mean(counts1_filtered) if counts1_filtered else 0
            st.markdown(f"""
            <div class="metric-card">
                <h3>{avg_genes:.1f}</h3>
                <p>Avg Genes/Group</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Bar chart
        st.markdown("### 📊 Gene Distribution by Age Group")
        
        # Create enhanced bar plot
        fig, ax = plt.subplots(figsize=(12, 6))
        apply_plot_style()
        
        bars = ax.bar(selected_ages, counts1_filtered, 
                     color='skyblue', alpha=0.8, edgecolor='navy', linewidth=1.2)
        
        # Add value labels on bars
        for i, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{int(height)}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_xlabel("Age Group", fontsize=12, fontweight='bold')
        ax.set_ylabel("Number of Switching Genes", fontsize=12, fontweight='bold')
        # TITOLO CON NOME PULITO
        ax.set_title(f"Distribution of Switching Genes by Age Group\n{tissue_name}", 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        st.pyplot(fig)
        
        # Download button for plot - NOME PULITO
        create_download_button(fig, f"gene_distribution_{tissue_name.replace(' ', '_').replace('-', '_')}.png")
        
        # Download CSV
        st.markdown("""
        <div class="download-section">
            <h4>📥 Download Data</h4>
        </div>
        """, unsafe_allow_html=True)
        
        csv1 = df1.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="⬇️ Download Complete Gene List (CSV)",
            data=csv1,
            file_name=f"switching_genes_{tissue_name.replace(' ', '_').replace('-', '_')}.csv",
            mime='text/csv'
        )
        
        # Gene lists by age group
        st.markdown("### 📋 Switching Genes by Age Group")
        
        for age in selected_ages:
            age_genes = df1[df1["Age"] == age]["Gene"].tolist()
            if age_genes:
                with st.expander(f"🔹 Age Group {age} – {len(age_genes)} genes"):
                    st.write(", ".join(age_genes))
                    
                    # Download button for this age group - NOME PULITO
                    age_df = pd.DataFrame(age_genes, columns=["Gene"])
                    csv_age = age_df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label=f"⬇️ Download {age} genes (CSV)",
                        data=csv_age,
                        file_name=f"genes_{age.replace('–','_')}_{tissue_name.replace(' ', '_').replace('-', '_')}.csv",
                        mime='text/csv',
                        key=f"download_{age}"
                    )
        
        # === INTRA-TISSUE HEATMAP ===
        if not df1.empty:
            st.markdown("""
            <div class="analysis-section">
                <h2>Intra-Tissue Heatmap - Age Group Overlap</h2>
            </div>
            """, unsafe_allow_html=True)
            
            # Create gene sets
            gene_sets = {
                age: set(df1[df1["Age"] == age]["Gene"])
                for age in age_groups
            }
            
            # Calculate overlap matrix — force int dtype
            matrix = np.zeros((len(age_groups), len(age_groups)), dtype=np.int64)
            for i, age_i in enumerate(age_groups):
                for j, age_j in enumerate(age_groups):
                    matrix[i, j] = int(len(gene_sets[age_i] & gene_sets[age_j]))

            # Build DataFrame for download
            df_plot = pd.DataFrame(matrix, index=age_groups, columns=age_groups)

            # FIX 2 — NON usare apply_plot_style() qui (interferisce con heatmap su mpl 3.10)

            # Dedicated figure
            fig, ax = plt.subplots(figsize=(12, 8))

            # FIX 1 — annot=False, annotazione MANUALE con ax.text
            sns.heatmap(
                matrix,
                annot=False,
                xticklabels=age_groups,
                yticklabels=age_groups,
                cmap="YlGnBu",
                ax=ax,
                cbar_kws={"label": "Number of Shared Genes"},
                square=True,
                linewidths=0.5,
                linecolor="white"
            )

            # Fix clipping (matplotlib 3.10 bug)
            bottom, top = ax.get_ylim()
            ax.set_ylim(bottom + 0.5, top - 0.5)

            # Annotazioni manuali — sempre visibili
            n = matrix.shape[0]
            for i in range(n):
                for j in range(n):
                    ax.text(j + 0.5, i + 0.5, f"{int(matrix[i, j])}",
                            ha="center", va="center", fontsize=14,
                            fontweight="bold", color="black")

            # FIX 4 — Titolo senza emoji
            ax.set_title(f"Gene Overlap Between Age Groups - {tissue_name}",
                        fontsize=16, fontweight='bold', pad=20)
            ax.set_xlabel("Age Group", fontsize=12, fontweight='bold')
            ax.set_ylabel("Age Group", fontsize=12, fontweight='bold')

            # FIX 3 — Layout/salvataggio stabile
            fig.tight_layout()

            # Download PRIMA di clear_figure (salva bytes in memoria)
            create_download_button(fig, f"age_overlap_heatmap_{tissue_name.replace(' ', '_').replace('-', '_')}.png")

            st.pyplot(fig, clear_figure=True)
            # Download overlap matrix
            df_matrix = pd.DataFrame(matrix, index=age_groups, columns=age_groups)
            csv_matrix = df_matrix.to_csv(index=True).encode('utf-8')
            st.download_button(
                label="⬇️ Download Overlap Matrix (CSV)",
                data=csv_matrix,
                file_name=f"overlap_matrix_{tissue_name.replace(' ', '_').replace('-', '_')}.csv",
                mime='text/csv'
            )
    
    # === TWO TISSUE COMPARISON ===
    elif file1 and file2:
        st.markdown("""
        <div class="analysis-section">
            <h2>📊 Comparison Between Two Tissues</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # NOMI PULITI USANDO LA FUNZIONE
        tissue1_name = extract_tissue_name(file1.name)
        tissue2_name = extract_tissue_name(file2.name)
        
        # Filter data
        df1_filtered = df1[df1["Age"].isin(selected_ages)]
        df2_filtered = df2[df2["Age"].isin(selected_ages)]
        
        counts1_filtered = [
            df1_filtered[df1_filtered["Age"] == age].shape[0]
            for age in age_groups if age in selected_ages
        ]
        counts2_filtered = [
            df2_filtered[df2_filtered["Age"] == age].shape[0]
            for age in age_groups if age in selected_ages
        ]
        
        # Comparison metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(df1)}</h3>
                <p>{tissue1_name}</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(df2)}</h3>
                <p>{tissue2_name}</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            total_unique = len(set(df1["Gene"]) | set(df2["Gene"]))
            st.markdown(f"""
            <div class="metric-card">
                <h3>{total_unique}</h3>
                <p>Unique Genes</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Comparison bar chart
        st.markdown("### 📊 Gene Count Comparison")
        
        fig, ax = plt.subplots(figsize=(12, 6))
        apply_plot_style()
        
        x = np.arange(len(selected_ages))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, counts1_filtered, width,
                      label=tissue1_name, 
                      color="skyblue", alpha=0.8, edgecolor='navy')
        bars2 = ax.bar(x + width/2, counts2_filtered, width,
                      label=tissue2_name, 
                      color="salmon", alpha=0.8, edgecolor='darkred')
        
        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                       f'{int(height)}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_xticks(x)
        ax.set_xticklabels(selected_ages)
        ax.set_title("Comparison of Switching Genes Between Tissues", 
                    fontsize=14, fontweight='bold', pad=20)
        ax.set_ylabel("Number of Switching Genes", fontsize=12, fontweight='bold')
        ax.set_xlabel("Age Group", fontsize=12, fontweight='bold')
        ax.legend(frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        st.pyplot(fig)
        
        # Download comparison plot
        create_download_button(fig, "tissue_comparison.png")
        
        # Download both datasets
        st.markdown("""
        <div class="download-section">
            <h4>📥 Download Datasets</h4>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        with col1:
            csv1 = df1.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=f"⬇️ Download {tissue1_name} (CSV)",
                data=csv1,
                file_name=f"dataset1_{tissue1_name.replace(' ', '_').replace('-', '_')}.csv",
                mime='text/csv'
            )
        
        with col2:
            csv2 = df2.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=f"⬇️ Download {tissue2_name} (CSV)",
                data=csv2,
                file_name=f"dataset2_{tissue2_name.replace(' ', '_').replace('-', '_')}.csv",
                mime='text/csv'
            )
    
    elif not file1:
        st.info("👆 Please upload at least one file to start the analysis.")
        
        # Show example format
        st.markdown("### 📋 Expected File Format")
        st.code("""
GENE1 GENE2 GENE3 GENE4
GENE5 GENE6 GENE7
GENE8 GENE9 GENE10 GENE11 GENE12
GENE13 GENE14
GENE15 GENE16 GENE17 GENE18
        """)
        st.caption("Each line represents an age group (30-39, 40-49, 50-59, 60-69, 70-79)")