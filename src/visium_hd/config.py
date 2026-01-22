from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
from typing_extensions import Literal

import scanpy as sc


@dataclass
class SampleInfo:
    name: str
    count_matrix: Path
    image: Path
    scalefactors: Path
    geojson: Path
    group: str


@dataclass
class AnalysisConfig:
    base_dir: Path = Path('.')
    data_dir: Path = field(default_factory=lambda: Path('data'))
    output_dir: Path = field(default_factory=lambda: Path('outputs'))

    min_cells_per_gene: int = 10
    min_counts_per_cell: int = 50
    max_counts_per_cell: Optional[int] = 6000
    max_mt_pct: float = 30.0
    min_genes_per_cell: int = 200
    max_genes_per_cell: Optional[int] = 5000

    n_top_genes: int = 3000
    n_pcs: int = 20
    n_neighbors: int = 20

    integration_method: Literal['harmony', 'bbknn', 'scvi', 'none'] = 'harmony'
    harmony_max_iter: int = 20

    clustering_resolutions: List[float] = field(default_factory=lambda: [0.1, 0.2, 0.4, 0.5, 0.6, 0.8])
    default_resolution: float = 0.5

    umap_min_dist: float = 0.5
    umap_spread: float = 2.0

    sketch_fraction: float = 0.25
    sketch_min_cells: int = 100

    random_state: int = 42

    def __post_init__(self):
        self.raw_dir = self.output_dir / "01_raw"
        self.filtered_dir = self.output_dir / "02_filtered"
        self.merged_dir = self.output_dir / "03_merged"
        self.plot_dir = self.output_dir / "04_plots"
        self.spatial_plot_dir = self.plot_dir / "05_spatial_plots"

        for d in [self.raw_dir, self.filtered_dir, self.merged_dir,
                  self.plot_dir, self.spatial_plot_dir]:
            d.mkdir(parents=True, exist_ok=True)

        sc.settings.verbosity = 2
        sc.settings.figdir = str(self.plot_dir)
        sc.set_figure_params(
            dpi=120, dpi_save=300, format='png',
            facecolor='white', frameon=False, figsize=(8, 8)
        )


def create_sample_config(data_dir: Path) -> Dict[str, SampleInfo]:
    samples = {
        "Colon_Cancer_P1": SampleInfo(
            name="Colon_Cancer_P1",
            count_matrix=data_dir / "Cancer_P1_filtered_feature_cell_matrix.h5",
            image=data_dir / "Cancer_P1_tissue_hires_image.png",
            scalefactors=data_dir / "Cancer_P1_scalefactors_json.json",
            geojson=data_dir / "Cancer_P1_cell_segmentations.geojson",
            group="Cancer"
        ),
        "Colon_Cancer_P2": SampleInfo(
            name="Colon_Cancer_P2",
            count_matrix=data_dir / "Cancer_P2_filtered_feature_cell_matrix.h5",
            image=data_dir / "Cancer_P2_tissue_hires_image.png",
            scalefactors=data_dir / "Cancer_P2_scalefactors_json.json",
            geojson=data_dir / "Cancer_P2_cell_segmentations.geojson",
            group="Cancer"
        ),
        "Colon_Normal_P3": SampleInfo(
            name="Colon_Normal_P3",
            count_matrix=data_dir / "Norm_P3_filtered_feature_cell_matrix.h5",
            image=data_dir / "Norm_P3_tissue_hires_image.png",
            scalefactors=data_dir / "Norm_P3_scalefactors_json.json",
            geojson=data_dir / "Norm_P3_cell_segmentations.geojson",
            group="Normal"
        ),
        "Colon_Normal_P5": SampleInfo(
            name="Colon_Normal_P5",
            count_matrix=data_dir / "Norm_P5_filtered_feature_cell_matrix.h5",
            image=data_dir / "Norm_P5_tissue_hires_image.png",
            scalefactors=data_dir / "Norm_P5_scalefactors_json.json",
            geojson=data_dir / "Norm_P5_cell_segmentations.geojson",
            group="Normal"
        ),
    }
    return samples
