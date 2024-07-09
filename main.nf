#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {CELLBENDER} from './modules/cellbender'
include {DOUBLETDETECTION} from './modules/doubletdetection'
include {AGGREGATION} from './modules/aggregation'
include {OUTLIER_FILTER} from './modules/outlierfilter'
include {SCRAN} from './modules/scran'
include {SCVI} from './modules/scvi'
include {POSTPROCESSING} from './modules/postprocessing'
include {LEIDEN_RESOLUTION} from './modules/leiden_resolution'
include {CELLTYPIST} from './modules/celltypist'
include {REPORT} from './modules/report'

workflow {
    // access the samplesheet
    sample_sheet = file(params.samplesheet)

    // read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def name = row[0]
        def raw_path = file(row[1])
        def filtered_path = row[2]
        def demultiplexing = row[3].toBoolean()
        return tuple(name, raw_path, filtered_path, demultiplexing)
    }

    // run Cellbender
    CELLBENDER(ch_input)

    // run DoubletDetection with optional demultiplexing
    DOUBLETDETECTION(CELLBENDER.out.name, CELLBENDER.out.raw_h5, CELLBENDER.out.filtered_path, CELLBENDER.out.demultiplexing)

    // aggregate the outputs
    AGGREGATION(DOUBLETDETECTION.out.doublet_h5ad.collect())

    // filter out outliers
    OUTLIER_FILTER(AGGREGATION.out.aggregation_h5ad)

    // SCRAN normalization
    SCRAN(OUTLIER_FILTER.out.outlier_filtered_h5ad)

    // SCVI batch correction
    SCVI(SCRAN.out.scran_h5ad)

    // Postprocessing
    POSTPROCESSING(SCVI.out.scvi_h5ad)

    // Celltypist
    CELLTYPIST(POSTPROCESSING.out.postprocessing_scvi_h5ad)

    // Leiden resolution search
    LEIDEN_RESOLUTION(POSTPROCESSING.out.postprocessing_scvi_h5ad)

    // Report
    REPORT(POSTPROCESSING.out.postprocessing_h5ad,
    CELLTYPIST.out.celltypist_scvi_h5ad)
}
