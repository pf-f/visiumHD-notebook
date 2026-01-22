from typing import Dict, List, Optional, Tuple

import spatialdata as spd
import spatialdata_plot as splt

from .config import AnalysisConfig


class SpatialVisualizer:
    def __init__(self, config: AnalysisConfig):
        self.config = config

    @staticmethod
    def crop_to_extent(
        sdata: spd.SpatialData,
        shape_key: str,
        coordinate_system: str = 'downscale_to_hires'
    ) -> spd.SpatialData:
        extent = spd.get_extent(
            sdata,
            elements=[shape_key],
            coordinate_system=coordinate_system
        )

        cropped = spd.bounding_box_query(
            sdata,
            min_coordinate=[extent['x'][0], extent['y'][0]],
            max_coordinate=[extent['x'][1], extent['y'][1]],
            axes=('x', 'y'),
            target_coordinate_system=coordinate_system
        )

        return cropped

    def plot_sample(
        self,
        sdata: spd.SpatialData,
        image_key: str,
        shape_key: str,
        color: str,
        title: Optional[str] = None,
        save_name: Optional[str] = None,
        coordinate_system: str = 'downscale_to_hires'
    ):
        if title is None:
            title = image_key.replace('_image', '')

        cropped = self.crop_to_extent(sdata, shape_key, coordinate_system)

        (
            cropped
            .pl.render_images(image_key)
            .pl.render_shapes(shape_key, color=color)
            .pl.show(
                coordinate_systems=coordinate_system,
                title=title,
                save=save_name
            )
        )

    def plot_all_samples(
        self,
        sdata: spd.SpatialData,
        color: str,
        suffix: str = ''
    ):
        image_keys = list(sdata.images.keys())
        shape_keys = list(sdata.shapes.keys())

        if len(image_keys) != len(shape_keys):
            print("Warning: Number of images and shapes don't match")

        for img_key, shp_key in zip(image_keys, shape_keys):
            sample_name = img_key.replace('_image', '')
            print(f"Plotting {sample_name} - {color}")

            save_path = self.config.spatial_plot_dir / f'{sample_name}_{color}{suffix}.png'

            try:
                self.plot_sample(
                    sdata, img_key, shp_key, color,
                    title=sample_name,
                    save_name=str(save_path)
                )
            except Exception as e:
                print(f"  Error: {e}")

    def plot_region_of_interest(
        self,
        sdata: spd.SpatialData,
        image_key: str,
        shape_key: str,
        bbox: Dict[str, Tuple[float, float]],
        color: str,
        title: str = 'ROI',
        coordinate_system: str = 'downscale_to_hires'
    ):
        cropped = spd.bounding_box_query(
            sdata,
            min_coordinate=[bbox['x'][0], bbox['y'][0]],
            max_coordinate=[bbox['x'][1], bbox['y'][1]],
            axes=('x', 'y'),
            target_coordinate_system=coordinate_system
        )

        (
            cropped
            .pl.render_images(image_key)
            .pl.render_shapes(shape_key, color=color)
            .pl.show(
                coordinate_systems=coordinate_system,
                title=title
            )
        )
