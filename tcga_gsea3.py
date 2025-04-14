import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gseapy import enrichr
from collections import defaultdict, Counter
import time
import numpy as np

# === CONFIGURATION === #
BASE_FOLDER = "./grade_generalised"
STAGES = ["StageI", "StageII", "StageIII", "StageIV"]
GENE_SETS = {
    "BP": "GO_Biological_Process_2021",
    "MF": "GO_Molecular_Function_2021",
    "CC": "GO_Cellular_Component_2021"
}
ALLOWED_EXTENSION = ".maf"
RESULT_DIR = "gsea_result_all_genes"
TOP_N = 10
RETRY_ATTEMPTS = 3
RETRY_DELAY = 5  # Seconds
IMPACTFUL_VARIANTS = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site']
TOP_TERMS_HEATMAP = 15  # For heatmap

# === INITIALIZATION === #
os.makedirs(RESULT_DIR, exist_ok=True)
writer = pd.ExcelWriter(os.path.join(RESULT_DIR, "gsea_all_stages.xlsx"), engine='openpyxl')
summary_matrix = defaultdict(dict)

def load_gene_list_from_file(file_path):
    """
    Extract unique genes from Hugo_Symbol in .maf files, filtering for impactful variants.
    """
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', low_memory=False)
        if 'Hugo_Symbol' not in df.columns:
            print(f"Error: 'Hugo_Symbol' column missing in {file_path}.")
            return []
        if 'Variant_Classification' in df.columns:
            df = df[df['Variant_Classification'].isin(IMPACTFUL_VARIANTS)]
        genes = df['Hugo_Symbol'].dropna().astype(str).str.strip()
        valid_genes = [g for g in genes if g and g != "Unknown" and g.isalnum()]
        if not valid_genes:
            print(f"No valid genes extracted from {file_path}.")
        return list(set(valid_genes))
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return []

def run_gsea(gene_list, gene_set, stage, ontology_key):
    """
    Run GSEA with retries and handle large gene lists.
    """
    if not gene_list:
        print(f"Skipping GSEA for {stage} ({ontology_key}): No genes found.")
        return pd.DataFrame()

    print(f"Running GSEA for {stage} ({ontology_key}) with {len(gene_list)} genes: {gene_list[:5]}...")

    for attempt in range(RETRY_ATTEMPTS):
        try:
            enr = enrichr(
                gene_list=gene_list,
                gene_sets=gene_set,
                organism="Human",
                outdir=None,
                cutoff=0.05,
                verbose=False
            )
            if enr.results is not None and not enr.results.empty:
                df = enr.results
                df["-log10(p-value)"] = df["P-value"].apply(lambda p: -np.log10(p) if p > 0 else 0)
                df["Ontology"] = ontology_key
                return df
            print(f"No significant results for {stage} ({ontology_key}).")
            return pd.DataFrame()
        except Exception as e:
            if attempt < RETRY_ATTEMPTS - 1:
                print(f"GSEA attempt {attempt + 1} failed for {stage} ({ontology_key}): {e}. Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
            else:
                print(f"GSEA failed for {stage} ({ontology_key}) after {RETRY_ATTEMPTS} attempts: {e}")
                return pd.DataFrame()

# === MAIN PIPELINE === #
for stage in STAGES:
    altered_genes = set()
    print(f"\n📂 Processing {stage}")
    stage_path = os.path.join(BASE_FOLDER, stage)
    if not os.path.exists(stage_path):
        print(f"Warning: Directory {stage_path} does not exist. Skipping.")
        continue

    for case_id in os.listdir(stage_path):
        case_path = os.path.join(stage_path, case_id)
        if not os.path.isdir(case_path):
            print(f"Warning: {case_path} is not a directory. Skipping.")
            continue
        for file in os.listdir(case_path):
            if file.endswith(ALLOWED_EXTENSION):
                file_path = os.path.join(case_path, file)
                genes = load_gene_list_from_file(file_path)
                if genes:
                    altered_genes.update(genes)
                else:
                    print(f"No genes extracted from {file_path}.")

    print(f"🔍 Found {len(altered_genes)} unique genes for {stage}")
    stage_dfs = []

    for key, gene_set in GENE_SETS.items():
        gsea_df = run_gsea(list(altered_genes), gene_set, stage, key)
        if not gsea_df.empty:
            top_df = gsea_df.sort_values("Adjusted P-value").head(TOP_N)
            top_df["Stage"] = stage
            stage_dfs.append(top_df)

            # Visualize Top N Enriched Terms
            plt.figure(figsize=(8, 5))
            sns.barplot(
                x="-log10(p-value)",
                y="Term",
                hue="Term",
                data=top_df,
                palette="viridis",
                legend=False
            )
            plt.title(f"{stage} - Top {TOP_N} {key} Enriched Terms")
            plt.tight_layout()
            plt.savefig(os.path.join(RESULT_DIR, f"{stage}_{key}_top_terms.png"))
            plt.close()

            # Populate matrix for heatmap
            for _, row in top_df.iterrows():
                term = row["Term"]
                score = row["Combined Score"]
                summary_matrix[term][stage] = score

    # Save stage-specific enrichment in Excel
    if stage_dfs:
        combined_df = pd.concat(stage_dfs)
        combined_df.to_excel(writer, sheet_name=stage, index=False)
    else:
        print(f"No GSEA results for {stage}. Skipping Excel sheet.")

# === HEATMAP OF TOP ENRICHED GO TERMS ACROSS STAGES === #
if summary_matrix:
    heatmap_df = pd.DataFrame(summary_matrix).T.fillna(0)
    # Select top terms by average Combined Score
    if len(heatmap_df) > TOP_TERMS_HEATMAP:
        heatmap_df['mean_score'] = heatmap_df.mean(axis=1)
        heatmap_df = heatmap_df.sort_values('mean_score', ascending=False).head(TOP_TERMS_HEATMAP)
        heatmap_df = heatmap_df.drop(columns='mean_score')
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        heatmap_df,
        cmap="YlGnBu",
        annot=True,
        fmt=".0f",
        linewidths=0.5,
        cbar_kws={'label': 'Combined Score'},
        annot_kws={'size': 8},
        xticklabels=True,
        yticklabels=True
    )
    plt.title("Top GO Terms Across Stages", fontsize=12)
    plt.xlabel("Stage", fontsize=10)
    plt.ylabel("GO Term", fontsize=10)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULT_DIR, "gsea_stage_comparison_heatmap.png"), dpi=500)
    plt.close()
else:
    print("No enrichment results for heatmap.")

# === SAVE EXCEL WORKBOOK === #
if writer.sheets:
    writer.close()
else:
    print("No sheets to save. Creating placeholder sheet.")
    pd.DataFrame({"Note": ["No GSEA results generated"]}).to_excel(writer, sheet_name="Summary")
    writer.close()

print("\n✅ GSEA analysis complete. Results saved to:", os.path.abspath(RESULT_DIR))