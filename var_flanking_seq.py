import primer3
import pandas as pd
import argparse
import pysam
from collections import Counter
from Bio import SeqIO
import random
import statistics
import sys

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="\n\
        Usage: ")
    parser.add_argument("-ref", "--reference", help = "reference genome. ")
    parser.add_argument("-gtf", "--gtf", required = False, help = "gtf annotation file matched to the reference genome. ")
    parser.add_argument("-vcf", "--vcf", help = "vcf variant file. ")
    parser.add_argument("-bam", "--bam", required = False, help = "mapping file in bam. ")
    parser.add_argument("-enzyme", "--enzyme", help="restriction enzyme table for selection.")
    parser.add_argument("-bed", "--bed", help = "BED file containing the targeted chromosome, start and end region. ")
    parser.add_argument("-primer_range", "--primer_range", default = 400, type = int, help = "primer distance from the cutting site (default: 400 bp). ")
    parser.add_argument("-CDS_search", "--CDS_search", default = False, action = 'store_true', help = "Search recognizable sites within gene CDS region. ")
    parser.add_argument("-out", "--output", help = "output prefix. ")
    parser.add_argument("-cov", "--coverage", required = False, type = int, help = "coverage estimation of bam file. ")
    args = parser.parse_args()
    return args

def enzyme_df_fun(args):
    '''parse restriction enzyme table'''
    enzyme = args.enzyme
    enzyme_df = pd.read_table(enzyme, sep = "\t", header = 0)
    enzyme_df = enzyme_df[["Enzyme", "Site", "label"]]
    enzyme_df["site_len"] = enzyme_df.apply(lambda row:len(row.Site), axis = 1)
    enzyme_df = enzyme_df.sort_values(by = "site_len")
    return enzyme_df 

def contig_seq(args):
    '''Parse fasta sequence file and build a dictionary'''
    ref_fasta = args.reference
    tig = args.contig
    seq_parse = SeqIO.parse(ref_fasta, "fasta")
    seq_dict = SeqIO.to_dict(seq_parse)
    seq_return = str(seq_dict[tig].seq) # all bases of the given contig
    return seq_return 

def vcf_df_fun(args):
    '''Parse vcf file and search for variant given the search region'''
    vcf = args.vcf
    vcf_handle = pysam.TabixFile(vcf, parser = pysam.asVCF())
    bed = args.bed
    bed_df = pd.read_table(bed, sep = "\t", header = None)
    for row in bed_df.itertuples():
        tig = row[1]
        pos_up = row[2]
        pos_down = row[3]
        tig_list = list()
        pos_list = list()
        ref_list = list()
        alt_list = list()
        # 
        for v in vcf_handle.fetch(tig, pos_up, pos_down):
            if v.filter == "PASS":
                tig = str(v.contig)
                ref = str(v.ref)
                alt = str(v.alt)
                pos = v.pos # 0-based variants
                tig_list.append(tig)
                pos_list.append(pos)
                ref_list.append(ref)
                alt_list.append(alt)
        # build a dataframe for enzyme
        df = pd.DataFrame({"contig":tig_list, "pos":pos_list, "ref":ref_list, "alt":alt_list})
        return df

def gff_df_fun(args):
    '''parse gff annotation file and return gene CDS range given the search region. '''
    gtf = args.gtf
    tig = args.contig
    pos_up = args.search_start
    pos_down = args.search_end
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    # 
    contig_list = list()
    gid_list = list()
    start_list = list()
    end_list = list()
    for gene in gtf_handle.fetch(tig, pos_up, pos_down):
        feature = gene.feature
        tig = gene.contig
        start = gene.start # 0-based
        end = gene.end # 1-based
        gid = gene.gene_id
        if feature == "CDS":
            contig_list.append(tig)
            gid_list.append(gid)
            start_list.append(start)
            end_list.append(end)
    df = pd.DataFrame({"contig":contig_list, "gid":gid_list, "start":start_list, "end":end_list})
    return df

def gff_gene_fun(args):
    '''parse gff annotation file and return gene region (including up- and down-stream 5000 bp) within the search region. '''
    gtf = args.gtf
    tig = args.contig
    pos_up = args.search_start
    pos_down = args.search_end
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    # 
    contig_list = list()
    gid_list = list()
    start_list = list()
    end_list = list()
    for gene in gtf_handle.fetch(tig, pos_up, pos_down):
        feature = gene.feature
        if feature == "gene":
            start = gene.start - 5000 # 0-based
            end = gene.end + 5000 # 1-based
            gid = gene.gene_id
            #
            contig_list.append(tig)
            gid_list.append(gid)
            start_list.append(start)
            end_list.append(end)
    df = pd.DataFrame({"contig":contig_list, "gid":gid_list, "start":start_list, "end":end_list})
    return df

def check_region(vcf_df, gff_df):
    '''Check the parsed variant and gff file and return if any variants within the gene coding region found'''
    select_vpos = list()
    for v in vcf_df.itertuples():
        vpos = v.pos # 0-based
        gff_sub = gff_df[(gff_df["start"] <= vpos) & (gff_df["end"] >= vpos)]
        if len(gff_sub) > 0:
            select_vpos.append(vpos)
    #
    if len(select_vpos) > 0:
        vcf_select = vcf_df[vcf_df.pos.isin(select_vpos)]
        #vcf_nselect = vcf_df[~vcf_df.pos.isin(select_vpos)]
        return vcf_select
    else:
        sys.exit("No variants found on gene coding region! ")

def ref_allele_validate(vcf, seq):
    '''Validate the ref-allele matched to the vcf file. '''
    for row in vcf.itertuples():
        vpos = row.pos
        ref = row.ref
        if seq[vpos] != ref:
            sys.exit("The reference sequence is not matched to the variant file. ")

def variant_seq(pos, tig_seq, recog_len, alt):
    '''Takes variant position and return sequences around it. '''
    ref_s = tig_seq[pos-recog_len+1:pos+recog_len]
    alt_s = tig_seq[pos-recog_len+1:pos] + alt + tig_seq[pos+1:pos+recog_len] # the variant always in the middel the extracted allele
    return ref_s, alt_s

def var_recog(args, vcf_df, enzyme_df, tig_seq):
    '''Identify recoginizing site by restriction enzyme. '''
    output = args.output
    recog_df = open(output+"_enzyme_recog_var.txt", "w")
    recog_df.write("pos\tref\talt\tenzyme\tsite\trecognize\n") # 0-based position
    for row in vcf_df.itertuples():
        vpos = row.pos
        palt = row.alt
        for row in enzyme_df.itertuples():
            enzyme = row.Enzyme
            site = row.Site
            site_len = row.site_len
            ref_s, alt_s = variant_seq(vpos, tig_seq, site_len, palt)
            if (site in ref_s) or (site in alt_s):
                if site in ref_s:
                    recog_df.write(f"{vpos}\t{ref_s}\t{alt_s}\t{enzyme}\t{site}\tref\n")
                else:
                    recog_df.write(f"{vpos}\t{ref_s}\t{alt_s}\t{enzyme}\t{site}\talt\n")
    recog_df.close()
    #
    recog_df = pd.read_table(output+"_enzyme_recog_var.txt", sep = "\t", header = 0)
    return recog_df

def var_distance(var_use):
    '''Variant sites recognized by cutting enzyme are distantly distributed.'''
    delete_pos_list = list()
    enzymes = sorted(set(var_use["enzyme"]))
    #################
    enzyme_clean = list()
    for en in enzymes:
        var_sub = var_use[var_use["enzyme"] == en]
        var_sort = var_sub.sort_values(by = ["pos"], ascending = True)
        var_sort["order"] = range(0, len(var_sort))
        init_pos = int(var_sort[var_sort["order"] == 0]["pos"])
        for i in range(1, len(var_sort)):
            var_pos = int(var_sort[var_sort["order"] == i]["pos"])
            if var_pos - init_pos < 300:
                delete_pos_list.append(init_pos)
                delete_pos_list.append(var_pos)
            init_pos = var_pos
        var_clean = var_sort[~var_sort["pos"].isin(delete_pos_list)]
        enzyme_clean.append(var_clean)
    var_distance_use = pd.concat(enzyme_clean, axis = 0)
    return var_distance_use

def bam_coverage(args, tig_seq):
    '''Estimate total coverage on the contig scale.'''
    if args.coverage == None:
        bam = args.bam
        tig = args.contig
        tig_len = len(tig_seq)
        bam_handle = pysam.AlignmentFile(bam, "rb")
        pos_rand = sorted(random.sample(range(0, tig_len, 10000), int(tig_len/10000-1)))
        coverage_random_check = list()
        for p in pos_rand:
            for pileupcolumn in bam_handle.pileup(tig, p, p+1):
                if pileupcolumn.pos == p:
                    coverage_random_check.append(pileupcolumn.n)
        coverage_soft_estimate = statistics.median(coverage_random_check)
        print("Estimated coverage value:", coverage_soft_estimate)
    else:
        coverage_soft_estimate = args.coverage
    return coverage_soft_estimate

def var_coverage(args, pos_s, coverage_estimate):
    '''Variant sites are ressured by its coverage. '''
    bam = args.bam
    tig = args.contig
    bam_handle = pysam.AlignmentFile(bam, "rb")
    pos_cov_validate = list()
    for p in pos_s:
        for pileupcolumn in bam_handle.pileup(tig, p, p+1):
            if pileupcolumn.pos == p:
                pos_cov = pileupcolumn.n
                if (pos_cov < int(coverage_estimate)*1.7) & (pos_cov > 5):
                    pos_cov_validate.append(p)
    return pos_cov_validate

def flanking_seq_info(recog_df, tig_seq, size, output):
    fw = open(output + ".txt", "w")
    fw.write("pos\tref\talt\tenzyme\tsite\trecognize\tsequence\n")
    for row in recog_df.itertuples():
        vpos = row.pos
        ref = row.ref
        alt = row.alt
        enzyme = row.enzyme
        site = row.site
        recognize = row.recognize
        start_site = vpos - size
        end_site = vpos + size + 1
        flank_s = tig_seq[start_site:end_site]
        fw.write(f"{vpos}\t{ref}\t{alt}\t{enzyme}\t{site}\t{recognize}\t{flank_s}\n")
    fw.close()

def main():
    args = args_parser()
    output = args.output
    tig_seq = contig_seq(args)
    primer_range = args.primer_range
    vcf_df = vcf_df_fun(args)
    print("Validate the reference fasta and vcf file matched ......")
    ref_allele_validate(vcf_df, tig_seq) # check the vcf file and ref-seq matched
    if args.CDS_search == True:
        print("CDS region ......")
        gff_df = gff_df_fun(args)
    else:
        print("gene region ......")
        gff_df = gff_gene_fun(args)
    ##
    enzyme_df = enzyme_df_fun(args)
    #
    print("Search variants on gene coding region ...... ")
    check_var = check_region(vcf_df, gff_df) # vcf of variants on gene-coding region
    #
    recog_df = var_recog(args, check_var, enzyme_df, tig_seq)
    recog_df = var_distance(recog_df) # variants distantly recognized
    #
    print("Check coverage at variant site ...... ")
    coverage = bam_coverage(args, tig_seq)
    pos_list = sorted(set(recog_df["pos"]))
    pos_new = var_coverage(args, pos_list, coverage) # proper coverage on the variant site
    recog_df_new = recog_df[recog_df.pos.isin(pos_new)]
    #
    print("Generate flanking sequences around the variant site ...... ")
    flanking_seq_info(recog_df_new, tig_seq, primer_range, output)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
