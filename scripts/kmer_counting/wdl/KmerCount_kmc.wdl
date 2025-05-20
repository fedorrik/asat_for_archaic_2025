version 1.0

workflow kmc_count_wf {

    call kmc_count

    output {
        File kmc_common_txt = kmc_count.kmc_common_txt_out
    }

    meta {
        author: "Fedor Ryabov"
        email: "fedorik1@gmail.com"
        description: "Fork from [kmcCount.wdl](https://github.com/kmiga/alphaAnnotation/blob/main/tools/kmcCount/kmcCount.wdl). [KMC](https://github.com/refresh-bio/KMC) DB creation and intersection with DB of interest. Outputs a text dump of intersection."
    }
}

task kmc_count {
    input {
        Array[File] input_reads
        Array[String] kmers
        String sample_name

        Int kmer_len    = 18
        Int min_cutoff  = 0
        Int max_count   = 10000000

        Int thread_count = 12
        Int mem_size_gb  = 36
        Int addl_disk    = 64
        Int preempts     = 2
    }

    parameter_meta {
        input_reads: "Reads to be merged and convertged to kmer db"
        kmers: "array of kmers"
        sample_name: "sample name to attach to output files"
        kmer_len: "length of kmer. Must match kmer length used to generate query DB."
        min_cutoff: "minimum count to include in output DB and text dump. kmers below this count are excluded."
        max_count: "maximum count to assign to a kmer. kmers above this count are assigned max_count counts."
    }

    Int input_size = ceil(size(input_reads, "GB"))
    Int final_disk_dize = input_size * 4 + addl_disk

    command <<<
        set -eux -o pipefail

        date +"%D %T"
        echo "Merging forward and reversed fastq files"
        cat ~{sep=' ' input_reads} > merged_reads.fastq.gz
        rm ~{sep=' ' input_reads}
        mkdir kmc_tmp_dir
        
        date +"%D %T"
        echo "Creating kmer_dbs"
        
        # create reads kmer db
        kmc \
            -t~{thread_count} \
            -k~{kmer_len} \
            -m~{mem_size_gb} \
            -ci~{min_cutoff} \
            -cs~{max_count} \
            merged_reads.fastq.gz \
            ~{sample_name} \
            kmc_tmp_dir/

        # kmer list to fasta
        for kmer in ~{sep=' ' kmers}; do
            echo ">kmer" >> query_kmers.fa
            echo $kmer >> query_kmers.fa
        done

        # create query kmer db
        kmc -fa -k18 -ci0 query_kmers.fa query_kmers .


        # interseect kmer dbs
        date +"%D %T"
        echo "Intersecting kmer_dbs"
        kmc_tools \
            -t~{thread_count} \
            simple \
            ~{sample_name} \
            query_kmers \
            intersect \
            common-kmers \
            -ocleft

        # dump the intersection to a text file
        kmc_tools \
            -t~{thread_count} \
            transform \
            common-kmers \
            dump \
            "~{sample_name}.txt"
        cat "~{sample_name}.txt"
    >>>

    output {
        File kmc_common_txt_out = "~{sample_name}.txt"
    }

    runtime {
        memory: mem_size_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "quay.io/biocontainers/kmc@sha256:41161e8b9da4a02a6e14ce1a9130474baba223c87213ba75d9092842d29400d1"
        preemptible: preempts
    }
}
