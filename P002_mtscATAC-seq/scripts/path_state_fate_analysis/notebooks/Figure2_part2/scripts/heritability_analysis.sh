#!/bin/bash
# Heritability analysis using phylogenetic trees
# Trees from: /lab-share/Hem-Sankaran-e2/Public/projects/tgao/plarry/phylosignal/final/trees

eval "$(mamba shell hook --shell bash)"
mamba activate base

home="/lab-share/Hem-Sankaran-e2/Public/projects/tgao"
wd="$home/plarry/phylosignal/final"
logdir="$wd/logs"
tree_dir="$wd/trees"
feature_mat="$wd/tf_auc_mat.csv"
n_cores=30

# Create log directory if it doesn't exist
mkdir -p "$logdir"

# Run heritability analysis for each tree type
for sample in \
    pL1000_CMP_phy_tru pL1000_CMP_phy_nj pL1000_CMP_phy_trim \
    pL1000_MEG_phy_tru pL1000_MEG_phy_nj pL1000_MEG_phy_trim \
    pL1000_MEP_phy_tru pL1000_MEP_phy_nj pL1000_MEP_phy_trim \
    pL1000_promyelo_phy_tru pL1000_promyelo_phy_nj pL1000_promyelo_phy_trim; do
    tree_file="$tree_dir/${sample}.newick"
    
    echo "Processing tree: $sample"
    
    sbatch --job-name="herit_${sample}" \
        --output="$logdir/herit_${sample}_%j.out" \
        --error="$logdir/herit_${sample}_%j.err" \
        --cpus-per-task=$n_cores \
        --mem=24G \
        --time=00:20:00 \
        $home/mitodrift_analysis/scripts/run_phylo_signal_wrapper.sh \
            --trait_mat $feature_mat \
            --tree $tree_file \
            --outfile $wd/${sample}_tf_cmeans.csv \
            --reps 9999 \
            --ncores $n_cores
done
