#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {CELLBENDER} from './modules/cellbender'
include {DOUBLETDETECTION} from './modules/doubletdetection'
include {AGGREGATION} from './modules/aggregation'
include {SCRAN} from './modules/scran'
include {SCVI} from './modules/scvi'

workflow {
    // access the samplesheet
    sample_sheet = file(params.samplesheet)
    
    // read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }
    
    // create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def name = row[0]
        def raw_path = file(row[1])
        def filtered_path = file(row[2])
        return tuple(name, raw_path, filtered_path)
    }
    
    // run Cellbender
    CELLBENDER(ch_input)
    
    // run DoubletDetection
    DOUBLETDETECTION(ch_input, CELLBENDER.out.cellbender_h5)
    
    // aggregate the outputs
    AGGREGATION(DOUBLETDETECTION.out.doublet_h5ad.collect())
    
    // SCRAN normalization
    SCRAN(AGGREGATION.out.aggregation_h5ad)
    
    // SCVI batch correction
    SCVI(SCRAN.out.scran_h5ad)
}
