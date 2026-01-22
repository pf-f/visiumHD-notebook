from typing import Dict, List, Optional

import scanpy as sc

from .config import AnalysisConfig


class CellAnnotator:
    MARKER_GENES = {
        "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
        "Malignant_Epithelial": ["EPCAM", "KRT8", "CEACAM5", "CEACAM6", "CLDN7"],
        "T_Cells": ["CD3D", "CD3E", "CD3G", "TRAC"],
        "B_Cells": ["CD79A", "CD79B", "MS4A1"],
        "Plasma_Cells": ["SDC1", "JCHAIN", "MZB1"],
        "NK_Cells": ["NCAM1", "NCR1", "KLRD1"],
        "Macrophages": ["CD68", "CD163", "MRC1", "C1QA"],
        "Monocytes": ["CD14", "FCGR3A", "S100A8", "S100A9"],
        "Dendritic_Cells": ["CD1C", "CLEC9A", "FLT3"],
        "Mast_Cells": ["KIT", "TPSAB1", "MS4A2"],
        "Fibroblasts": ["DCN", "LUM", "COL1A1", "COL1A2", "PDPN"],
        "CAFs": ["FAP", "PDPN", "ACTA2", "PDGFRB"],
        "Endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5"],
        "Cycling": ["MKI67", "TOP2A", "CENPF"]
    }

    def __init__(self, config: AnalysisConfig):
        self.config = config

    def find_marker_genes(
        self,
        adata: sc.AnnData,
        groupby: str = 'clusters',
        method: str = 'wilcoxon'
    ) -> sc.AnnData:
        print(f"Finding marker genes for {groupby}...")

        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method=method,
            tie_correct=True
        )

        sc.pl.rank_genes_groups_dotplot(
            adata,
            groupby=groupby,
            standard_scale='var',
            n_genes=10,
            save='_marker_genes_dotplot.png'
        )

        df = sc.get.rank_genes_groups_df(adata, group=None, pval_cutoff=0.05)
        df.to_csv(self.config.merged_dir / 'marker_genes.csv', index=False)

        return adata

    def plot_known_markers(
        self,
        adata: sc.AnnData,
        groupby: str = 'clusters',
        markers: Optional[Dict[str, List[str]]] = None
    ):
        if markers is None:
            markers = self.MARKER_GENES

        filtered_markers = {}
        missing_genes = []

        for cell_type, genes in markers.items():
            valid = [g for g in genes if g in adata.var_names]
            if valid:
                filtered_markers[cell_type] = valid
            missing = [g for g in genes if g not in adata.var_names]
            missing_genes.extend(missing)

        if missing_genes:
            print(f"Missing {len(set(missing_genes))} genes: {list(set(missing_genes))[:10]}...")

        sc.pl.dotplot(
            adata,
            var_names=filtered_markers,
            groupby=groupby,
            standard_scale='var',
            save='_known_markers_dotplot.png'
        )

    def run_celltypist(
        self,
        adata: sc.AnnData,
        model_name: str = 'Human_Colorectal_Cancer.pkl',
        color_cols=None
    ) -> sc.AnnData:
        import celltypist
        from celltypist import models

        print(f"Running CellTypist with model: {model_name}")

        model = models.Model.load(model=model_name)
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=True
        )

        adata = predictions.to_adata()

        if color_cols is None:
            color_cols = ['clusters', 'leiden_0.1', 'leiden_0.2']
        if isinstance(color_cols, str):
            color_cols = [color_cols]
        color = ['predicted_labels', 'majority_voting'] + color_cols

        sc.pl.umap(
            adata,
            color=color,
            frameon=True,
            ncols=2,
            legend_loc='on data',
            legend_fontsize=8,
            save='_celltypist_annotation.png'
        )

        return adata

    def run_singler(
        self,
        adata: sc.AnnData,
        singler_ref: str = 'blueprint_encode',
        ref_version: str = '2024-02-26',
        groupby: str = 'clusters',
        color_cols=None
    ) -> sc.AnnData:
        import celldex
        import singler
        import scranpy
        from summarizedexperiment import SummarizedExperiment

        print(f"Running SingleR with reference: {singler_ref}")

        if adata.raw is not None:
            print("  Using adata.raw.X for SingleR annotation.")
            mat = adata.raw.X.T
            features = list(adata.raw.var_names)
        else:
            raise ValueError("Error: adata.raw does not exist!")

        test = scranpy.aggregate_across_cells(
            mat,
            [adata.obs[groupby]]
        ).sum

        test = SummarizedExperiment(
            assays={'counts': test},
            row_names=features,
            column_names=list(adata.obs[groupby].cat.categories)
        )

        ref_data = celldex.fetch_reference(
            singler_ref, ref_version, realize_assays=True
        )

        results = singler.annotate_single(
            test_data=test,
            test_features=features,
            ref_data=ref_data,
            ref_labels='label.main'
        )

        cluster_to_label = dict(zip(
            list(test.column_data.row_names),
            results['best']
        ))
        adata.obs['singler'] = adata.obs[groupby].map(cluster_to_label)

        if color_cols is None:
            color_cols = ['clusters', 'leiden_0.1', 'leiden_0.2']
        if isinstance(color_cols, str):
            color_cols = [color_cols]
        color = ['singler'] + color_cols

        sc.pl.umap(
            adata,
            color=color,
            legend_loc='on data',
            legend_fontsize=8,
            save='_singler_annotation.png'
        )

        return adata

    def manual_annotation(
        self,
        adata: sc.AnnData,
        cluster_to_celltype: Dict[str, str],
        cluster_key: str = 'clusters',
        annotation_key: str = 'celltype'
    ) -> sc.AnnData:
        adata.obs[annotation_key] = (
            adata.obs[cluster_key]
            .astype(str)
            .map(cluster_to_celltype)
            .astype('category')
        )

        sc.pl.umap(
            adata,
            color=[cluster_key, annotation_key],
            legend_loc='on data'
        )

        return adata
