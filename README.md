# PCR_primer

PCR (polymerase chain reaction) primer design is a commonly used and fundamental skill in laboratories, and is important for genotyping, variant screening and many other researches. Here, I developed a pipeline for the design of PCR primers as used in genotyping purposes. 

There are a lot of genetic variations in different populations. By taking advantage the modern high-throughput sequencing technology, we are able to identify the genetic variations within and between species. These identified variants can then be used as markers for following genotyping purposes of unknown individual. For example, there are two maize inbred lines, and based on sequencing and variant calling process we know the genetic variants between them. The crosses design using the two parental maize inbreds generated a large population we call it as recombinant inbred lines (RILs). Because of recombination events happened in the reproducing process, the genetic compositions in RILs are different from each other. 

The question is: how can we use available information and resource to know the genotype at a given loci for these RILs? One simple (and quick) and also easy-to-do method is PCR the target loci. 

Input files you have: 
1. Reference genome (in fasta format)
2. VCF file (a standard file format for the storage of variants information)
3. BED file (a standard file format for gene location information) 
4. 
