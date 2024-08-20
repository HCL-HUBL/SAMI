process splicing_filter {
    cpus { params.CPU_splicing }
    label 'multicore'
    label 'retriable'
    publishDir "${params.out}/splicing", mode: params.publish

    input:
    path(exons)
    path('*')
    path('*')
    path('*')
    path(targetGTF)
    val(dir)

    output:
    path("${dir}"), emit: splicing_output

    """
    Rscript --vanilla "${projectDir}/scripts/splicing_filter.R" ${params.CPU_splicing} "$targetGTF" "$exons" ${params.plot} ${params.fusions} ${params.min_I} ${params.min_PSI} "$params.symbols" "$params.classes" "$params.focus"
    """
}
