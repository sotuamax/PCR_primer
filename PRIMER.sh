#!/bin/bash

function usage {
    printf "\nUsage: PRIMER.sh [options] \n"
    printf "\n"
    printf "OPTIONS\n"
    printf "  -ref|--reference <fasta>	reference fasta file\n"
    printf "  -gtf|--gtf <gtf>	gtf annotation file for the reference genome\n"
    printf "  -vcf|--vcf <vcf>	vcf file for variants candidate\n"
    printf "  -bam|--bam <bam>	bam file of DNA-mapping against the reference genome\n"
    printf "  -enzyme|--enzymte <enzyme.txt>	Enzyme table with name and recognizable site\n"
    printf "  -tig|--contig <tig_name>	contig name for primer design\n"
    printf "  -start|--search_start <int>	start searching region on the contig\n"
    printf "  -end|--search_end <int>	end searching region on the contig\n"
    printf "  -primer_range|--primer_range <int>	region around the variant for primer design\n"
    printf "  -out|--output <prefix>	prefix name for the output file\n"
    printf "  -cov|--coverage <int>	estimated coverage based on pre-knowledge\n"
    printf "  -CDS|--CDS_search <boolean> If only searching for primers using gene CDS region\n"
    printf "  -h|--help	Display help information\n"
    printf "  -v|--verbose	write log information to standard output rather than log file\n"
    printf "\n"
    exit 1
}

# Handle arguments

RANGE=300
COV=NULL
VERBOSE=true
CDS=true

while [[ $# -gt 0 ]]; do
    case $1 in
        -ref|--reference)
            REF=$2
            shift; shift ;;
        -gtf|--gtf) 
            GTF=$2
            shift; shift ;;
        -vcf|--vcf)
            VCF=$2
            shift; shift ;;
        -bam|--bam)
            BAM=$2
            shift; shift ;;
        -enzyme|--enzyme)
            ENZYME=$2
            shift; shift ;;
        -tig|--contig)
            TIG=$2
            shift; shift ;;
        -start|--search_start)
            START=$2
            shift; shift ;;
        -end|--search_end)
            END=$2
            shift; shift ;;
        -primer_range|--primer_range)
            RANGE=$2
            shift; shift ;;
        -out|--output)
            OUT=$2
            shift; shift ;;
        -cov|--coverage)
            COV=$2
            shift; shift ;;
        -CDS|--CDS_search)
            CDS=$2
            shift; shift ;;
        -v|--verbose)
            VERBOSE=$2
            shift; shift ;;
        -h|--help)
            usage; break ;;
        --) shift ; break ;;
    esac
done


# Run code
if [ "$COV" = NULL ]; then 
    python var_flanking_seq.py -ref $REF -gtf $GTF -vcf $VCF -bam $BAM -enzyme $ENZYME -tig $TIG -search_start $START -search_end $END -primer_range $RANGE -out $OUT -CDS_search $CDS 2>&1 
    python primer_explore.py -flanking_seq $OUT.txt -primer_range $RANGE -out $OUT 2>&1
else
    python var_flanking_seq.py -ref $REF -gtf $GTF -vcf $VCF -bam $BAM -enzyme $ENZYME -tig $TIG -search_start $START -search_end $END -primer_range $RANGE -out $OUT -cov $COV -CDS_search $CDS 2>&1
    python primer_explore.py -flanking_seq $OUT.txt -primer_range $RANGE -out $OUT 2>&1
fi
