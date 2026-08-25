#!/bin/bash
set -e
mkdir -p /mnt/e/genotype_data/liftover
if [ ! -x "$HOME/miniconda3/bin/conda" ]; then
  echo "[setup] downloading miniconda..."
  curl -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
  echo "[setup] installing miniconda to $HOME/miniconda3..."
  bash /tmp/miniconda.sh -b -p "$HOME/miniconda3" > /tmp/miniconda_install.log 2>&1
fi
export PATH="$HOME/miniconda3/bin:$PATH"
echo "[setup] conda version:"
conda --version
echo "[setup] installing mamba..."
conda install -y -n base -c conda-forge mamba > /tmp/mamba_install.log 2>&1 || echo "[setup] mamba install failed, will fall back to conda"
echo "[setup] creating liftover env (minimap2 samtools bcftools crossmap)..."
if [ -x "$HOME/miniconda3/bin/mamba" ]; then
  "$HOME/miniconda3/bin/mamba" create -y -n liftover -c conda-forge -c bioconda minimap2 samtools bcftools crossmap > /tmp/liftover_env.log 2>&1
else
  "$HOME/miniconda3/bin/conda" create -y -n liftover -c conda-forge -c bioconda minimap2 samtools bcftools crossmap > /tmp/liftover_env.log 2>&1
fi
echo "[setup] DONE"
for t in minimap2 samtools bcftools CrossMap CrossMap.py; do
  p="$HOME/miniconda3/envs/liftover/bin/$t"
  if [ -x "$p" ]; then echo "[setup] OK $t"; else echo "[setup] MISSING $t"; fi
done
