#!/usr/bin/env python
"""
Visium HD Spatial Transcriptomics Analysis Pipeline

Usage:
    visium-hd run --data-dir DATA --output-dir OUTPUT [options]
    visium-hd run-full --data-dir DATA --output-dir OUTPUT [options]
    visium-hd --help

Options:
    --data-dir DIR          Input data directory [default: data]
    --output-dir DIR        Output directory [default: outputs]
    --batch-method METHOD   Batch correction method: harmony, bbknn, scvi, none [default: harmony]
    --max-mt-pct PCT        Maximum mitochondrial percentage [default: 30.0]
    --resolution RES        Clustering resolution [default: 0.5]
    --use-celltypist        Use CellTypist for annotation [default: True]
    --use-singler           Use SingleR for annotation [default: True]
    --skip-qc               Skip quality control step
    --help                  Show this help message
"""

import argparse
import sys
from pathlib import Path

from visium_hd import (
    AnalysisConfig,
    SpatialTranscriptomicsPipeline,
    create_sample_config,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Visium HD Spatial Transcriptomics Analysis Pipeline'
    )
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    run_parser = subparsers.add_parser('run', help='Run analysis pipeline')
    run_full_parser = subparsers.add_parser('run-full', help='Run full pipeline with all steps')

    for p in [run_parser, run_full_parser]:
        p.add_argument('--data-dir', type=Path, default=Path('data'),
                       help='Input data directory')
        p.add_argument('--output-dir', type=Path, default=Path('outputs'),
                       help='Output directory')
        p.add_argument('--batch-method', type=str, default='harmony',
                       choices=['harmony', 'bbknn', 'scvi', 'none'],
                       help='Batch correction method')
        p.add_argument('--max-mt-pct', type=float, default=30.0,
                       help='Maximum mitochondrial percentage')
        p.add_argument('--resolution', type=float, default=0.5,
                       help='Default clustering resolution')
        p.add_argument('--use-celltypist', action='store_true', default=True,
                       help='Use CellTypist for annotation')
        p.add_argument('--no-celltypist', dest='use_celltypist', action='store_false',
                       help='Skip CellTypist')
        p.add_argument('--use-singler', action='store_true', default=True,
                       help='Use SingleR for annotation')
        p.add_argument('--no-singler', dest='use_singler', action='store_false',
                       help='Skip SingleR')
        p.add_argument('--skip-qc', action='store_true',
                       help='Skip quality control step')

    return parser.parse_args()


def run_pipeline(args):
    config = AnalysisConfig(
        base_dir=Path('.'),
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        max_mt_pct=args.max_mt_pct,
        integration_method=args.batch_method,
        default_resolution=args.resolution
    )

    samples = create_sample_config(config.data_dir)
    pipeline = SpatialTranscriptomicsPipeline(config, samples)

    if args.skip_qc:
        (pipeline
            .run_data_preparation()
            .run_preprocessing()
            .run_final_integration(method=args.batch_method)
        )
    else:
        (pipeline
            .run_data_preparation()
            .run_qc_and_filtering()
            .run_preprocessing()
            .run_final_integration(method=args.batch_method)
        )

    pipeline.run_marker_gene_analysis(groupby='leiden_0.2')

    if args.use_celltypist:
        pipeline.run_celltypist_annotation()

    if args.use_singler:
        pipeline.run_singler_annotation()

    print("\n" + "!"*60)
    print("PIPELINE PAUSED - Manual Annotation Required")
    print("!"*60)
    print("\nPlease review outputs/04_plots for annotation results,")
    print("then provide manual_mapping and run:")
    print("  pipeline.run_manual_annotation(manual_mapping)")
    print("  pipeline.run_spatial_visualization()")
    print("  pipeline.run_differential_expression()")
    print("  pipeline.save_final_results()")


def run_full_pipeline(args):
    config = AnalysisConfig(
        base_dir=Path('.'),
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        max_mt_pct=args.max_mt_pct,
        integration_method=args.batch_method,
        default_resolution=args.resolution
    )

    samples = create_sample_config(config.data_dir)
    pipeline = SpatialTranscriptomicsPipeline(config, samples)

    pipeline.run_full_pipeline(
        batch_method=args.batch_method,
        use_celltypist=args.use_celltypist,
        use_singler=args.use_singler
    )


def main():
    args = parse_args()

    if args.command is None:
        print(__doc__)
        sys.exit(1)

    if args.command == 'run':
        run_pipeline(args)
    elif args.command == 'run-full':
        run_full_pipeline(args)


if __name__ == '__main__':
    main()
