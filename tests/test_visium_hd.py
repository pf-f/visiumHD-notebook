import sys
sys.path.insert(0, 'src')

import pytest
from pathlib import Path
from dataclasses import dataclass

from visium_hd import (
    AnalysisConfig,
    SampleInfo,
    create_sample_config,
)


class TestConfig:
    def test_analysis_config_creation(self):
        config = AnalysisConfig(
            base_dir=Path('.'),
            data_dir=Path('data'),
            output_dir=Path('outputs'),
            max_mt_pct=30.0
        )
        
        assert config.max_mt_pct == 30.0
        assert config.integration_method == 'harmony'
        assert config.default_resolution == 0.5

    def test_analysis_config_directories(self, tmp_path):
        config = AnalysisConfig(
            base_dir=tmp_path,
            data_dir=tmp_path / 'data',
            output_dir=tmp_path / 'outputs'
        )
        
        assert config.raw_dir.exists()
        assert config.filtered_dir.exists()
        assert config.merged_dir.exists()
        assert config.plot_dir.exists()
        assert config.spatial_plot_dir.exists()

    def test_sample_info_creation(self):
        sample = SampleInfo(
            name='Test_Sample',
            count_matrix=Path('data/test.h5'),
            image=Path('data/test.png'),
            scalefactors=Path('data/test.json'),
            geojson=Path('data/test.geojson'),
            group='Test'
        )
        
        assert sample.name == 'Test_Sample'
        assert sample.group == 'Test'

    def test_create_sample_config(self, tmp_path):
        samples = create_sample_config(tmp_path / 'data')
        assert isinstance(samples, dict)


class TestQCProcessor:
    def test_qc_processor_init(self):
        config = AnalysisConfig()
        from visium_hd import QCProcessor
        
        qc = QCProcessor(config)
        assert qc.config == config


class TestBatchCorrector:
    def test_batch_corrector_init(self):
        config = AnalysisConfig()
        from visium_hd import BatchCorrector
        
        bc = BatchCorrector(config)
        assert bc.config == config


class TestCellAnnotator:
    def test_cell_annotator_init(self):
        config = AnalysisConfig()
        from visium_hd import CellAnnotator
        
        ca = CellAnnotator(config)
        assert ca.config == config

    def test_marker_genes_dict_exists(self):
        from visium_hd import CellAnnotator
        
        assert hasattr(CellAnnotator, 'MARKER_GENES')
        assert isinstance(CellAnnotator.MARKER_GENES, dict)


class TestSpatialVisualizer:
    def test_spatial_visualizer_init(self):
        config = AnalysisConfig()
        from visium_hd import SpatialVisualizer
        
        sv = SpatialVisualizer(config)
        assert sv.config == config


class TestDEAnalyzer:
    def test_de_analyzer_init(self):
        config = AnalysisConfig()
        from visium_hd import DEAnalyzer
        
        de = DEAnalyzer(config)
        assert de.config == config


class TestPipeline:
    def test_pipeline_init(self):
        config = AnalysisConfig()
        samples = create_sample_config(Path('data'))
        from visium_hd import SpatialTranscriptomicsPipeline
        
        pipeline = SpatialTranscriptomicsPipeline(config, samples)
        assert pipeline.config == config
        assert pipeline.samples == samples


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
