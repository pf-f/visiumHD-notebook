from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scanpy.external as sce

from .config import AnalysisConfig


class BatchCorrector:
    def __init__(self, config: AnalysisConfig):
        self.config = config

    def geosketch_subsample(self, adata: sc.AnnData) -> sc.AnnData:
        import geosketch

        print("Running Geosketch subsampling...")

        sketched_list = []
        for sample in adata.obs['sample'].cat.categories:
            subset = adata[adata.obs['sample'] == sample].copy()
            n_cells = subset.n_obs

            n_sample = max(
                self.config.sketch_min_cells,
                int(n_cells * self.config.sketch_fraction)
            )
            n_sample = min(n_sample, n_cells)

            print(f"  {sample}: {n_sample}/{n_cells} cells")

            sketch_idx = geosketch.gs(
                X=subset.obsm['X_pca'],
                N=n_sample,
                seed=self.config.random_state
            )
            sketched_list.append(subset[sketch_idx].copy())

        sketched_adata = sc.concat(sketched_list)
        print(f"Sketched data: {sketched_adata.n_obs} cells")

        return sketched_adata

    def run_integration(
        self,
        adata: sc.AnnData,
        method: str = 'harmony',
        batch_key: str = 'sample',
        save: bool = True
    ) -> sc.AnnData:
        print(f"Running batch correction: {method}")
        result = adata.copy()

        if method == 'harmony':
            result = self._run_harmony(result, batch_key)
        elif method == 'bbknn':
            result = self._run_bbknn(result, batch_key)
        elif method == 'scvi':
            result = self._run_scvi(result, batch_key)
        elif method == 'none':
            result = self._run_uncorrected(result)
        else:
            raise ValueError(f"Unknown integration method: {method}")

        if save:
            save_dir = self.config.merged_dir / "integrated_methods"
            save_dir.mkdir(parents=True, exist_ok=True)
            save_path = save_dir / f"sketched_adata_{method}.h5ad"
            print(f"Saving integrated result to: {save_path}")
            result.write(save_path, compression='gzip')

        return result

    def _run_harmony(self, adata: sc.AnnData, batch_key: str) -> sc.AnnData:
        sce.pp.harmony_integrate(
            adata,
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_pca_harmony',
            max_iter_harmony=self.config.harmony_max_iter
        )

        adata.obsm['X_pca_orig'] = adata.obsm['X_pca'].copy()

        sc.pp.neighbors(
            adata,
            n_neighbors=self.config.n_neighbors,
            n_pcs=self.config.n_pcs,
            use_rep='X_pca_harmony',
            metric='correlation'
        )

        return adata

    def _run_bbknn(self, adata: sc.AnnData, batch_key: str) -> sc.AnnData:
        sce.pp.bbknn(
            adata,
            batch_key=batch_key,
            n_pcs=self.config.n_pcs,
            neighbors_within_batch=15
        )
        return adata

    def _run_scvi(self, adata: sc.AnnData, batch_key: str) -> sc.AnnData:
        import scvi

        scvi.model.SCVI.setup_anndata(
            adata,
            layer='counts',
            batch_key=batch_key
        )

        model = scvi.model.SCVI(
            adata,
            n_layers=2,
            n_latent=30,
            gene_likelihood='nb'
        )
        model.train(datasplitter_kwargs={"num_workers": 30}, early_stopping=True)

        adata.obsm['X_scVI'] = model.get_latent_representation()

        sc.pp.neighbors(
            adata,
            n_neighbors=self.config.n_neighbors,
            use_rep='X_scVI'
        )

        return adata

    def _run_uncorrected(self, adata: sc.AnnData) -> sc.AnnData:
        sc.pp.neighbors(
            adata,
            n_neighbors=self.config.n_neighbors,
            n_pcs=self.config.n_pcs,
            use_rep='X_pca',
            metric='correlation'
        )
        return adata

    def cluster_and_visualize(
        self,
        adata: sc.AnnData,
        method_name: str = 'integrated',
        multi_resolution: bool = True
    ) -> sc.AnnData:
        sc.tl.umap(
            adata,
            min_dist=self.config.umap_min_dist,
            spread=self.config.umap_spread,
            random_state=self.config.random_state
        )

        if multi_resolution:
            resolutions = self.config.clustering_resolutions
            print(f"  Multi-resolution clustering: {resolutions}")
        else:
            resolutions = [self.config.default_resolution]
            print(f"  Single-resolution clustering for comparison: {resolutions}")

        for res in resolutions:
            sc.tl.leiden(
                adata,
                flavor='igraph',
                resolution=res,
                key_added=f'leiden_{res}',
                random_state=self.config.random_state
            )

        adata.obs['clusters'] = adata.obs[f'leiden_{self.config.default_resolution}']

        sc.pl.umap(
            adata,
            color=['sample', 'group', 'clusters'],
            frameon=False,
            ncols=2,
            wspace=0.3,
            save=f'_{method_name}_umap.png'
        )

        self._plot_cluster_distribution(adata, method_name)

        return adata

    def _plot_cluster_distribution(self, adata: sc.AnnData, method_name: str):
        crosstab = pd.crosstab(adata.obs['sample'], adata.obs['clusters'])

        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(crosstab.values, cmap='hot', interpolation='nearest')

        ax.set_xticks(range(crosstab.shape[1]))
        ax.set_xticklabels(crosstab.columns)
        ax.set_yticks(range(crosstab.shape[0]))
        ax.set_yticklabels(crosstab.index)

        ax.set_xlabel('Cluster')
        ax.set_ylabel('Sample')
        ax.set_title(f'Cell Distribution - {method_name}')

        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.savefig(
            self.config.plot_dir / f'{method_name}_cluster_distribution.png',
            dpi=300
        )
        plt.close()

    def compare_methods(
        self,
        adata: sc.AnnData,
        methods: Optional[List[str]] = None
    ) -> Dict[str, sc.AnnData]:
        results = {}

        if methods is None:
            methods = ['none', 'harmony', 'bbknn', 'scvi']

        for method in methods:
            print(f"\n{'='*50}")
            print(f"Testing method: {method}")
            print('='*50)

            adata_method = self.run_integration(adata, method=method)
            adata_method = self.cluster_and_visualize(adata_method, method_name=method, multi_resolution=False)
            results[method] = adata_method

        return results
