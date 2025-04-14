#!/usr/bin/env nextflow

// Initialize the input channel with the tuple [A, [1, 2, 3]]
Channel
    .of( ['A', [1, 2, 3]] )
    .flatMap { item ->
        def key = item[0]
        def values = item[1]
        // Create a new list of tuples combining the key with each value
        return values.collect { value -> [key, value] }
    }
    .view()
