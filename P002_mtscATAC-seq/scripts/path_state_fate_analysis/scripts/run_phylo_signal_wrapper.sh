#!/bin/bash
# Wrapper script to run phylo signal analysis with proper environment

eval "$(mamba shell hook --shell bash)"
mamba activate base

Rscript /lab-share/Hem-Sankaran-e2/Public/projects/tgao/mitodrift_analysis/scripts/run_phylo_signal.R "$@"
