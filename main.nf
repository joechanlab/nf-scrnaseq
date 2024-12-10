#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXTRACT_RNA } from './modules/extract_rna'
include { CELLBENDER } from './modules/cellbender'
include { DOUBLETDETECTION } from './modules/doubletdetection'
include { AGGREGATION } from './modules/aggregation'
include { OUTLIER_FILTER } from './modules/outlierfilter'
include { SCRAN } from './modules/scran'
include { SCVI } from './modules/scvi'
include { POSTPROCESSING } from './modules/postprocessing'
include { CELLTYPIST } from './modules/celltypist'
include { REPORT } from './modules/report'
include { SEACELLS } from './modules/seacells'

workflow {
    // Access the samplesheet
    sample_sheet = file(params.samplesheet)

    // Read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // Create a channel from the paths
    out = Channel.from(sample_sheet_data).map { row ->
        def (name, raw_path, filtered_path, demultiplexing, expected_droplets) = row
        return tuple(name, file(raw_path), filtered_path, demultiplexing.toLowerCase(), expected_droplets)
    }

    // Run Cellbender
    if (params.atac) {
        EXTRACT_RNA(out)
        out = EXTRACT_RNA.out.output
    }

    // run Cellbender
    if (!params.cellbender.skip) {
        CELLBENDER(out)
        out = CELLBENDER.out.output
    }

    // Run DoubletDetection
    DOUBLETDETECTION(out)

    // Optional: aggregate the outputs
    if (params.aggregation) {
        AGGREGATION(DOUBLETDETECTION.out.output_h5ad.collect())
        outlier_input = AGGREGATION.out
    } else {
        outlier_input = DOUBLETDETECTION.out
    }

    // Filter out outliers
    OUTLIER_FILTER(outlier_input.name, outlier_input.output_h5ad)

    // SCRAN normalization
    SCRAN(OUTLIER_FILTER.out.name, OUTLIER_FILTER.out.outlier_filtered_h5ad)

    // SCVI batch correction
    SCVI(SCRAN.out.name, SCRAN.out.scran_h5ad)

    // Postprocessing
    POSTPROCESSING(SCVI.out.name, SCVI.out.scvi_h5ad)

    // Celltypist and Report
    if (params.celltypist.skip) {
        REPORT(POSTPROCESSING.out.name, POSTPROCESSING.out.postprocessing_h5ad, POSTPROCESSING.out.postprocessing_scvi_h5ad)
    } else {
        CELLTYPIST(POSTPROCESSING.out.name, POSTPROCESSING.out.postprocessing_scvi_h5ad, POSTPROCESSING.out.postprocessing_h5ad)
        REPORT(CELLTYPIST.out.name, CELLTYPIST.out.postprocessing_h5ad, CELLTYPIST.out.celltypist_scvi_h5ad)
    }

    // SEACells
    if (params.atac) {
        SEACELLS(POSTPROCESSING.out.name, POSTPROCESSING.out.postprocessing_h5ad)
    }

}
