from typing import Dict, List, Optional, Tuple

import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from .config import AnalysisConfig


class DEAnalyzer:
    def __init__(self, config: AnalysisConfig):
        self.config = config

    def pseudobulk_de(
        self,
        adata: sc.AnnData,
        celltype: str,
        condition_col: str = 'group',
        contrast: Tuple[str, str, str] = ('group', 'Cancer', 'Normal'),
        celltype_col: str = 'celltype',
        sample_col: str = 'sample'
    ) -> pd.DataFrame:
        print(f"Running pseudobulk DE for {celltype}...")

        aggregated = sc.get.aggregate(
            adata,
            by=[celltype_col, sample_col],
            func=['sum'],
            layer='counts'
        )

        subset = aggregated[aggregated.obs[celltype_col] == celltype].copy()

        if subset.n_obs < 4:
            print(f"  Warning: Only {subset.n_obs} samples, DE may be unreliable")

        metadata = subset.obs[[sample_col]].copy()
        metadata[condition_col] = (
            metadata[sample_col]
            .str.split('_')
            .str[0]
        )

        counts_df = pd.DataFrame(
            subset.layers['sum'],
            index=metadata.index,
            columns=subset.var_names
        ).astype(int)

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors=condition_col,
            refit_cooks=True
        )
        dds.deseq2()

        stats = DeseqStats(dds, contrast=list(contrast))
        stats.summary()

        results = stats.results_df.sort_values('padj')

        output_path = self.config.merged_dir / f'DE_{celltype}_{contrast[1]}_vs_{contrast[2]}.csv'
        results.to_csv(output_path)
        print(f"  Saved to {output_path}")

        return results

    def visualize_de_genes(
        self,
        adata: sc.AnnData,
        de_results: pd.DataFrame,
        top_n: int = 10,
        visualizer=None,
        sdata=None
    ):
        top_up = (
            de_results[de_results['log2FoldChange'] > 0]
            .sort_values('padj')
            .head(top_n)
        )

        top_down = (
            de_results[de_results['log2FoldChange'] < 0]
            .sort_values('padj')
            .head(top_n)
        )

        print(f"\nTop {top_n} upregulated genes:")
        print(top_up[['log2FoldChange', 'padj']])

        print(f"\nTop {top_n} downregulated genes:")
        print(top_down[['log2FoldChange', 'padj']])

        top_genes = list(top_up.index[:5]) + list(top_down.index[:5])
        top_genes = [g for g in top_genes if g in adata.var_names]

        if top_genes:
            sc.pl.umap(
                adata,
                color=top_genes,
                ncols=3,
                save='_de_top_genes.png'
            )

        if visualizer and sdata and top_genes:
            for gene in top_genes[:3]:
                visualizer.plot_all_samples(sdata, color=gene, suffix='_DE')
