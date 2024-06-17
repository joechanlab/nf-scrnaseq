#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {CELLBENDER} from './modules/cellbender'

workflow {
    // access the samplesheet
    sample_sheet = file(params.samplesheet)
    
    // read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }
    
    // create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def name = row[0]
        def path = file(row[1])
        def expected_cells = row[2].trim().toInteger()
        return tuple(name, path, expected_cells)
    }
    
    // run Cellbender
    CELLBENDER(ch_input)
    
    // run DoubletDetection
    DOUBLETDETECTION(ch_input, CELLBENDER.out.cellbender_h5)
}
