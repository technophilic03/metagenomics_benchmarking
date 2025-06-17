# Load dependencies
module load samtools
module load bwa/0.7.17
module load htslib
module load boost
module load miniconda

# Set location of conda env
loc=/restricted/projectnb/pathoscope/work/aubrey/newMB_02_25_meta_benchmark/my_conda_env

# Install metaenv (DEPENDENCY)
conda install -p $loc -c bioconda -c conda-forge metasnv

# Test metasnv install
metaSNV.py --help
metaSNV_Filtering.py --help
metaSNV_DistDiv.py --help
metaSNV_subpopr.R --help

# Install mOTUs
mamba activate $loc
pip install motu-profiler

# Download database
motus downloadDB

# Test motus install
motus profile --test