#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {CELLBENDER} from './modules/cellbender'

workflow {
    main:
    ch_input = Channel.fromPath(params.input, checkIfExists: true)

    CELLBENDER(ch_input)
}
