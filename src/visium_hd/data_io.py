import gc
import json
from pathlib import Path
from typing import Dict, List, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as spd
from PIL import Image
from shapely.geometry import Polygon
from spatialdata.models import Image2DModel, ShapesModel, TableModel
from spatialdata.transformations import Identity, Scale

from .config import SampleInfo


class DataIO:
    @staticmethod
    def build_single_sample_zarr(
        sample: SampleInfo,
        output_dir: Path,
        overwrite: bool = True
    ) -> Path:
        print(f"Processing sample: {sample.name}")

        adata = sc.read_10x_h5(sample.count_matrix)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample.name
        adata.obs['group'] = sample.group
        adata.obs.index = f"{sample.name}_" + adata.obs.index.astype(str)

        img_array = DataIO._process_image(sample.image)

        with open(sample.scalefactors, 'r') as f:
            scale_data = json.load(f)
        hires_scale = scale_data.get('tissue_hires_scalef', 1.0)

        shapes_gdf, valid_indices = DataIO._process_geojson(
            sample.geojson, sample.name, adata.obs.index
        )

        n_original = len(adata)
        adata = adata[valid_indices].copy()
        print(f"  Matched {len(adata)}/{n_original} cells with shapes")

        adata.obs['cell_id'] = adata.obs.index
        adata.obs['region'] = f"{sample.name}_shapes"
        adata.obs['region'] = adata.obs['region'].astype('category')

        transform_shapes = {
            "downscale_to_hires": Scale([hires_scale, hires_scale], axes=("x", "y"))
        }
        transform_image = {"downscale_to_hires": Identity()}

        image_key = f"{sample.name}_image"
        shapes_key = f"{sample.name}_shapes"
        table_key = "table"

        sdata = spd.SpatialData(
            images={
                image_key: Image2DModel.parse(
                    img_array, transformations=transform_image
                )
            },
            shapes={
                shapes_key: ShapesModel.parse(
                    shapes_gdf, transformations=transform_shapes
                )
            },
            tables={
                table_key: TableModel.parse(
                    adata,
                    region=shapes_key,
                    region_key='region',
                    instance_key='cell_id'
                )
            }
        )

        zarr_path = output_dir / f"{sample.name}.zarr"
        sdata.write(zarr_path, overwrite=overwrite)
        print(f"  Saved to {zarr_path}")

        del sdata, adata, shapes_gdf, img_array
        gc.collect()

        return zarr_path

    @staticmethod
    def _process_image(image_path: Path) -> np.ndarray:
        img = np.array(Image.open(image_path))
        if img.ndim == 2:
            img = img[np.newaxis, :, :]
        elif img.ndim == 3:
            img = np.transpose(img, (2, 0, 1))
        return img

    @staticmethod
    def _process_geojson(
        geojson_path: Path,
        sample_name: str,
        obs_index: pd.Index
    ) -> Tuple[gpd.GeoDataFrame, List[str]]:
        with open(geojson_path, 'r') as f:
            geojson_data = json.load(f)

        geojson_map = {}
        for feature in geojson_data['features']:
            cell_id = feature['properties']['cell_id']
            formatted_id = f"{sample_name}_cellid_{cell_id:09d}-1"
            coords = feature['geometry']['coordinates'][0]
            geojson_map[formatted_id] = Polygon(coords)

        valid_indices = [idx for idx in obs_index if idx in geojson_map]
        geometries = [geojson_map[idx] for idx in valid_indices]

        shapes_gdf = gpd.GeoDataFrame(
            {'geometry': geometries},
            index=valid_indices
        )

        return shapes_gdf, valid_indices

    @staticmethod
    def concatenate_samples(
        zarr_paths: List[Path],
        samples: Dict[str, SampleInfo],
        output_path: Path
    ) -> spd.SpatialData:
        print("Reading and concatenating samples...")

        sdata_list = []
        for zarr_path in zarr_paths:
            sdata = spd.read_zarr(zarr_path)
            sample_name = zarr_path.stem.replace('.zarr', '')

            for table in sdata.tables.values():
                table.var_names_make_unique()
                if sample_name in samples:
                    table.obs['group'] = samples[sample_name].group

            sdata_list.append(sdata)

        merged_sdata = spd.concatenate(sdata_list, concatenate_tables=True)
        merged_sdata.write(output_path, overwrite=True)
        print(f"Merged data saved to {output_path}")

        del sdata_list
        gc.collect()

        return merged_sdata
