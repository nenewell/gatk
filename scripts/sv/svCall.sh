#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [5] absolute path to the reference fasta on the cluster (2bit is assumed accompanying with same basename and extension \".2bit\", skiplist with extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [6] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "Example:"
    echo -e " bash svDiscover.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-test-cluster \\"
    echo -e "      /test-sample \\"
    echo -e "      /data/NA12878_test.bam \\"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta.img"
    exit 1
fi

cd $(dirname -- $"0")
./scanBam.sh "$@"
./callVariants.sh "$1" "$2" "$3" "$5"
