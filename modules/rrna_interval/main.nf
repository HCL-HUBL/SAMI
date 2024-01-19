process rrna_interval {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir params.store, mode: "copy"

    input:
    path(GTF)
    path(chrNameLength)

    output:
    path("${GTF}.rRNA"), emit: rRNAs

    """
    # Header (consider unsorted to be safe)
    echo -e "@HD\tVN:1.0\tSO:unsorted" > "./${GTF}.rRNA"

    # Chromosomes (from STAR genome)
    sed -r 's/^(.+)\t(.+)\$/@SQ\tSN:\\1\tLN:\\2/' "${chrNameLength}" >> "./${GTF}.rRNA"

    # BED-like content
    grep -E 'transcript_(bio)?type "rRNA"' "$GTF" | awk -F "\t" '\$3 == "transcript" { id=gensub(/^.+transcript_id \"([^\"]+)\";.+\$/, "\\\\1", "g", \$9); print \$1"\t"\$4"\t"\$5"\t"\$7"\t"id }' >> "./${GTF}.rRNA"
    """
}
