import scanpy as sc

from .config import AnalysisConfig


class QCProcessor:
    def __init__(self, config: AnalysisConfig):
        self.config = config

    def calculate_qc_metrics(self, adata: sc.AnnData) -> sc.AnnData:
        print("Calculating QC metrics...")

        adata.var['mt'] = adata.var_names.str.match(r'^MT[-]')
        adata.var['ribo'] = adata.var_names.str.match(r'^RP[SL]')
        adata.var['hb'] = adata.var_names.str.match(r'^HB[^Pp]')

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mt', 'ribo', 'hb'],
            percent_top=None,
            log1p=True,
            inplace=True
        )

        return adata

    def plot_qc_metrics(
        self,
        adata: sc.AnnData,
        suffix: str = "",
        per_sample: bool = False
    ):
        metrics = [
            'n_genes_by_counts', 'total_counts',
            'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb'
        ]

        sc.pl.violin(
            adata,
            keys=metrics,
            groupby='sample',
            stripplot=False,
            inner='box',
            rotation=90,
            multi_panel=False,
            save=f'_qc_violin_{suffix}.png',
            show=True
        )

        if per_sample:
            for sample_id in adata.obs['sample'].unique():
                adata_sub = adata[adata.obs['sample'] == sample_id]

                sc.pl.violin(
                    adata_sub,
                    keys=metrics,
                    jitter=0.4,
                    multi_panel=True,
                    save=f'_{sample_id}_{suffix}_violin.png',
                    show=True
                )

                sc.pl.highest_expr_genes(
                    adata_sub,
                    n_top=20,
                    save=f'_{sample_id}_{suffix}_top20_genes.png',
                    show=True
                )

    def filter_cells_and_genes(self, adata: sc.AnnData) -> sc.AnnData:
        print(f"Before filtering: {adata.shape[0]} cells, {adata.shape[1]} genes")

        sc.pp.filter_genes(adata, min_cells=self.config.min_cells_per_gene)
        sc.pp.filter_cells(adata, min_counts=self.config.min_counts_per_cell)
        sc.pp.filter_cells(adata, min_genes=self.config.min_genes_per_cell)

        if self.config.max_counts_per_cell:
            adata = adata[
                adata.obs['total_counts'] <= self.config.max_counts_per_cell
            ].copy()

        adata = adata[
            adata.obs['pct_counts_mt'] <= self.config.max_mt_pct
        ].copy()

        print(f"After filtering: {adata.shape[0]} cells, {adata.shape[1]} genes")

        return adata

    def preprocess(self, adata: sc.AnnData) -> sc.AnnData:
        print("Running preprocessing...")

        adata.layers['counts'] = adata.X.copy()
        adata.raw = adata.copy()

        sc.pp.normalize_total(adata, target_sum=None)
        sc.pp.log1p(adata)

        sc.tl.pca(adata, n_comps=self.config.n_pcs, svd_solver='arpack')

        sc.pl.pca_variance_ratio(
            adata, log=True, n_pcs=self.config.n_pcs,
            save='_pca_elbow.png'
        )

        return adata
