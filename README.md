# Visium HD Spatial Transcriptomics Analysis

A analysis pipeline for Visium HD spatial transcriptomics data, featuring sample processing, quality control, batch correction, cell annotation, spatial visualization, and differential expression analysis.

## Environment Setup

### Option 1: Using Conda (Recommended)

This project contains both Python and R code. Use the conda environment file:

```bash
# Create environment from YAML file
conda env create -f environment.yml

# Activate the environment
conda activate visium_hd

# Install Python dependencies
pip install -r requirements.txt

# Install additional R packages in R session
R
> install.packages(c(
    "Seurat", "Signac", "monocle3",
    "patchwork", "ggpubr", "circlize",
    "clusterProfiler", "DESeq2"
  ))
> quit()
```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/visium_hd.git
cd visium_hd

# Follow environment setup above, then:

# Install Python package
pip install -e .
```

## Quick Start

### Using Python Script

```python
from pathlib import Path
from visium_hd import AnalysisConfig, create_sample_config, SpatialTranscriptomicsPipeline

# Create configuration
config = AnalysisConfig(
    base_dir=Path('.'),
    data_dir=Path('data'),
    output_dir=Path('outputs'),
    max_mt_pct=30.0,
    integration_method='harmony',
    default_resolution=0.5
)

# Create sample configuration
samples = create_sample_config(config.data_dir)

# Initialize pipeline
pipeline = SpatialTranscriptomicsPipeline(config, samples)

# Run full pipeline
pipeline.run_full_pipeline(
    batch_method='harmony',
    manual_annotation_mapping={
        '0': 'Epithelial',
        '1': 'Fibroblast',
        '2': 'T_Cells',
        # ... add more mappings
    },
    de_celltypes=['Epithelial', 'Macrophages']
)
```


## Project Structure

```
visium_hd/
├── src/
│   └── visium_hd/
│       ├── __init__.py
│       ├── config.py          # Configuration management
│       ├── data_io.py         # Data I/O utilities
│       ├── qc_processor.py    # QC and preprocessing
│       ├── batch_corrector.py # Batch correction
│       ├── cell_annotator.py  # Cell annotation
│       ├── spatial_visualizer.py  # Spatial visualization
│       ├── de_analyzer.py     # Differential expression
│       └── pipeline.py        # Main pipeline orchestrator
├── scripts/
│   └── run_analysis.py        # Command-line entry point
├── examples/
│   └── V1X.ipynb              # Jupyter notebook workflow
├── tests/
│   └── test_*.py              # Unit tests
├── data/                      # Input data directory
├── outputs/                   # Output directory
├── requirements.txt
├── pyproject.toml
└── README.md
```

## Configuration

Key parameters can be configured in `AnalysisConfig`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_cells_per_gene` | 10 | Minimum cells per gene for filtering |
| `min_counts_per_cell` | 50 | Minimum counts per cell |
| `max_mt_pct` | 30.0 | Maximum mitochondrial percentage |
| `n_top_genes` | 3000 | Number of highly variable genes |
| `integration_method` | 'harmony' | Batch correction method |
| `clustering_resolutions` | [0.1, 0.2, 0.4, 0.5, 0.6, 0.8] | Clustering resolutions |

## Dependencies

- Python 3.12
- scanpy
- spatialdata
- pydeseq2
- geosketch
- geopandas
- matplotlib

