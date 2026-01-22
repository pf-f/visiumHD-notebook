from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import scanpy as sc

from .batch_corrector import BatchCorrector
from .cell_annotator import CellAnnotator
from .config import AnalysisConfig, SampleInfo, create_sample_config
from .data_io import DataIO
from .de_analyzer import DEAnalyzer
from .qc_processor import QCProcessor
from .spatial_visualizer import SpatialVisualizer


class SpatialTranscriptomicsPipeline:
    def __init__(self, config: AnalysisConfig, samples: Dict[str, SampleInfo]):
        self.config = config
        self.samples = samples

        self.qc = QCProcessor(config)
        self.batch = BatchCorrector(config)
        self.annotator = CellAnnotator(config)
        self.visualizer = SpatialVisualizer(config)
        self.de = DEAnalyzer(config)

        self.sdata = None
        self.adata = None

    def run_data_preparation(self):
        print("\n" + "="*60)
        print("STEP 1: Data Preparation")
        print("="*60)

        zarr_paths = []
        for name, sample in self.samples.items():
            path = DataIO.build_single_sample_zarr(
                sample, self.config.raw_dir
            )
            zarr_paths.append(path)

        merged_path = self.config.raw_dir / 'merged_sdata.zarr'
        self.sdata = DataIO.concatenate_samples(
            zarr_paths, self.samples, merged_path
        )

        self.adata = self.sdata.tables['table'].copy()

        return self

    def run_qc_and_filtering(self, plot_per_sample: bool = True):
        print("\n" + "="*60)
        print("STEP 2: QC and Filtering")
        print("="*60)

        self.adata = self.qc.calculate_qc_metrics(self.adata)
        self.qc.plot_qc_metrics(self.adata, suffix='raw', per_sample=plot_per_sample)
        self.adata = self.qc.filter_cells_and_genes(self.adata)
        self.qc.plot_qc_metrics(self.adata, suffix='filtered', per_sample=plot_per_sample)

        self.adata.write(self.config.filtered_dir / 'filtered_adata.h5ad')

        return self

    def run_preprocessing(self):
        print("\n" + "="*60)
        print("STEP 3: Preprocessing")
        print("="*60)

        self.adata = self.qc.preprocess(self.adata)
        self.adata.write(self.config.merged_dir / 'preprocessed_adata.h5ad')

        return self

    def run_compare_batch_methods(
        self,
        methods: Optional[List[str]] = None
    ) -> Dict[str, sc.AnnData]:
        print("\n" + "="*60)
        print("STEP 4a: Compare Batch Correction Methods (Sketched Data)")
        print("="*60)
        print(f"Using single resolution: {self.config.default_resolution}")

        sketched = self.batch.geosketch_subsample(self.adata)
        results = self.batch.compare_methods(sketched, methods=methods)

        sketched.write(self.config.merged_dir / 'sketched_adata.h5ad')

        print("\n" + "="*60)
        print("COMPARISON COMPLETE!")
        print("="*60)
        print(f"\nPlease review UMAP plots in: {self.config.plot_dir}")
        print("\nAfter deciding the best method, run:")
        print("  pipeline.run_final_integration(method='your_chosen_method')")

        return results

    def run_final_integration(self, method: str = 'harmony'):
        print("\n" + "="*60)
        print("STEP 4b: Final Integration on Full Data")
        print("="*60)
        print(f"Selected method: {method}")
        print(f"Using multi-resolution: {self.config.clustering_resolutions}")

        self.adata = self.batch.run_integration(self.adata, method=method)
        self.adata = self.batch.cluster_and_visualize(
            self.adata,
            method_name=f'final_{method}',
            multi_resolution=True
        )

        self.adata.write(self.config.merged_dir / f'final_{method}_integrated_adata.h5ad')

        return self

    def run_marker_gene_analysis(self, groupby: str = 'clusters'):
        print("\n" + "="*60)
        print(f"Running Marker Gene Analysis (Groupby: {groupby})")
        print("="*60)

        self.annotator.find_marker_genes(self.adata, groupby=groupby)
        self.annotator.plot_known_markers(self.adata, groupby=groupby)

        print(f"✅ Marker gene analysis completed")
        return self

    def run_celltypist_annotation(
        self,
        celltypist_model: str = 'Human_Colorectal_Cancer.pkl'
    ):
        print("\n" + "="*60)
        print(f"Running CellTypist Annotation (Model: {celltypist_model})")
        print("="*60)

        try:
            self.adata = self.annotator.run_celltypist(self.adata, model_name=celltypist_model)
        except Exception as e:
            print(f"❌ CellTypist failed: {e}")

        print("✅ CellTypist annotation completed")
        return self

    def run_singler_annotation(
        self,
        singler_ref: str = 'blueprint_encode',
        ref_version: str = '2024-02-26',
        groupby: str = 'clusters'
    ):
        print("\n" + "="*60)
        print(f"Running SingleR Annotation (Ref: {singler_ref})")
        print("="*60)

        try:
            self.adata = self.annotator.run_singler(
                self.adata,
                singler_ref=singler_ref,
                ref_version=ref_version,
                groupby=groupby
            )
        except Exception as e:
            print(f"❌ SingleR failed: {e}")

        print("✅ SingleR annotation completed")
        return self

    def run_manual_annotation(
        self,
        manual_mapping: Dict[str, str],
        cluster_key: str = 'clusters',
        annotation_key: str = 'celltype'
    ):
        print("\n" + "="*60)
        print("STEP 5b: Manual Annotation")
        print("="*60)

        if cluster_key not in self.adata.obs.columns:
            available = [c for c in self.adata.obs.columns if 'leiden_0.2' in c.lower() or 'cluster' in c.lower()]
            if available:
                cluster_key = available[0]
                print(f"  Using cluster column: {cluster_key}")
            else:
                raise ValueError(f"Cluster column '{cluster_key}' not found.")

        clusters_in_data = set(self.adata.obs[cluster_key].astype(str).unique())
        clusters_in_mapping = set(manual_mapping.keys())

        missing = clusters_in_data - clusters_in_mapping
        extra = clusters_in_mapping - clusters_in_data

        if missing:
            print(f"  ⚠️  Warning: Missing clusters in mapping: {missing}")
            print(f"      These will be labeled as 'Unknown'")
        if extra:
            print(f"  ⚠️  Warning: Extra clusters in mapping (not in data): {extra}")

        self.adata = self.annotator.manual_annotation(
            self.adata,
            cluster_to_celltype=manual_mapping,
            cluster_key=cluster_key,
            annotation_key=annotation_key
        )

        print(f"\n注释结果 ('{annotation_key}')：")
        print(self.adata.obs[annotation_key].value_counts())

        self.adata.write(self.config.merged_dir / 'annotated_adata.h5ad')
        print(f"\n  Saved: {self.config.merged_dir / 'annotated_adata.h5ad'}")

        return self

    def run_spatial_visualization(self, colors: Optional[List[str]] = None):
        if colors is None:
            colors = ['clusters', 'celltype']

        print("\n" + "="*60)
        print("STEP 6: Spatial Visualization")
        print("="*60)

        self.sdata.tables['table'] = self.adata

        for color in colors:
            if color in self.adata.obs.columns or color in self.adata.var_names:
                self.visualizer.plot_all_samples(self.sdata, color=color)
            else:
                print(f"  Warning: {color} not found, skipping")

        return self

    def run_differential_expression(
        self,
        celltypes: Optional[List[str]] = None,
        contrast: Tuple[str, str, str] = ('group', 'Cancer', 'Normal')
    ):
        print("\n" + "="*60)
        print("STEP 7: Differential Expression Analysis")
        print("="*60)

        if celltypes is None and 'celltype' in self.adata.obs.columns:
            celltypes = self.adata.obs['celltype'].unique().tolist()
        elif celltypes is None:
            print("No celltype column found, skipping DE")
            return self

        for ct in celltypes:
            try:
                results = self.de.pseudobulk_de(
                    self.adata, celltype=ct, contrast=contrast
                )
                self.de.visualize_de_genes(
                    self.adata, results,
                    visualizer=self.visualizer, sdata=self.sdata
                )
            except Exception as e:
                print(f"  DE failed for {ct}: {e}")

        return self

    def save_final_results(self):
        print("\n" + "="*60)
        print("Saving Final Results")
        print("="*60)

        self.adata.write(self.config.merged_dir / 'final_adata.h5ad')

        final_sdata_path = self.config.merged_dir / 'final_sdata.zarr'
        self.sdata.tables['table'] = self.adata
        self.sdata.write(final_sdata_path, overwrite=True)

        for sample in self.adata.obs['sample'].unique():
            subset = self.adata[self.adata.obs['sample'] == sample]

            export_df = pd.DataFrame({
                'Barcode': subset.obs.index.str.split('_').str[-1],
                'Cluster': subset.obs['clusters'],
            })

            if 'celltype' in subset.obs.columns:
                export_df['CellType'] = subset.obs['celltype']

            export_df.to_csv(
                self.config.merged_dir / f'{sample}_annotations.csv',
                index=False
            )

        print("All results saved successfully!")

        return self

    def run_full_pipeline(
        self,
        manual_annotation_mapping: Optional[Dict[str, str]] = None,
        de_celltypes: Optional[List[str]] = None,
        batch_method: str = 'harmony',
        skip_manual_annotation: bool = False,
        use_celltypist: bool = True,
        celltypist_model: str = 'Human_Colorectal_Cancer.pkl',
        use_singler: bool = True,
        singler_ref: str = 'blueprint_encode',
        ref_version: str = '2024-02-26',
    ):
        (self
            .run_data_preparation()
            .run_qc_and_filtering()
            .run_preprocessing()
            .run_final_integration(method=batch_method)
        )

        self.run_marker_gene_analysis(groupby='leiden_0.2')

        if use_celltypist:
            self.run_celltypist_annotation(celltypist_model=celltypist_model)

        if use_singler:
            self.run_singler_annotation(
                singler_ref=singler_ref,
                ref_version=ref_version
            )

        if manual_annotation_mapping:
            self.run_manual_annotation(manual_annotation_mapping)
        elif not skip_manual_annotation:
            print("\n" + "!"*60)
            print("PIPELINE PAUSED - Manual Annotation Required")
            print("!"*60)
            print("\n查看自动注释结果，准备manual_mapping后运行：")
            print("""
            # 继续流程：
            manual_mapping = {
            '0': 'CellType_A',
            '1': 'CellType_B',
            # ...
            }
            pipeline.run_manual_annotation(manual_mapping)
            pipeline.run_spatial_visualization()
            pipeline.run_differential_expression()
            pipeline.save_final_results()
            """)
            return self

        (self
            .run_spatial_visualization()
            .run_differential_expression(celltypes=de_celltypes)
            .save_final_results()
        )

        print("\n" + "="*60)
        print("Pipeline completed successfully!")
        print("="*60)

        return self
