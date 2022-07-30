import primer3
import pandas as pd
import argparse
import sys

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="\n\
        Usage: ")
    parser.add_argument("-flanking_seq", "--flanking_sequence", help="flanking sequence for primer designation.")
    parser.add_argument("-primer_range", "--primer_range", default = 300, type = int, help = "primer distance from the cutting site. ")
    parser.add_argument("-out", "--output", help = "output prefix. ")
    args = parser.parse_args()
    return args

def GC_content(seq):
    GC_prop = (seq.count("C") + seq.count("G"))/len(seq)*100
    return GC_prop

def primer_info(args):
    '''Write primer design results into output file. '''
    seq = args.flanking_sequence
    primer_range = args.primer_range
    output = args.output
    seq_df = pd.read_table(seq, header = 0, sep = "\t")
    fw = open(output + "_PRIMER.txt", "w")
    for row in seq_df.itertuples():
        vpos = row.pos # 0-based
        vseq = row.sequence
        enzyme = row.enzyme
        recognize = row.recognize
        site = row.site
        label = enzyme + "_" + str(vpos)
        fw.write(f"#### {label}\n")
        design_candidate = primer_design(vseq, label, primer_range)
        number_return = design_candidate["PRIMER_PAIR_NUM_RETURNED"]
        if number_return > 0:
            for i in range(0, number_return):
                fw.write(f"######## Candidate_{i} {site} {recognize}\n")
                left_seq = "PRIMER_LEFT_" + str(i) + "_SEQUENCE"
                right_seq = "PRIMER_RIGHT_" + str(i) + "_SEQUENCE"
                left_pos = "PRIMER_LEFT_" + str(i)
                right_pos = "PRIMER_RIGHT_" + str(i)
                left_tm = "PRIMER_LEFT_" + str(i) + "_TM"
                right_tm = "PRIMER_RIGHT_" + str(i) + "_TM"
                product_size = "PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"
                #
                left_start = int(str(design_candidate[left_pos]).split(", ")[0].strip("("))
                right_end = int(str(design_candidate[right_pos]).split(", ")[0].strip("("))
                pcr_product = vseq[left_start:right_end+1]
                pcr_gc = GC_content(pcr_product)
                # 
                ref_cut_fragment_list = pcr_product.split(site)
                fw.write(f"\tCG_content_{i} : {pcr_gc}\n")
                fw.write(f"\tLEFT_{i} : {design_candidate[left_seq]}\n")
                fw.write(f"\tRIGHT_{i} : {design_candidate[right_seq]}\n")
                fw.write(f"\tPCR_product{i} : {pcr_product}\n")
                left_homodimer = primer3.calcHomodimer(design_candidate[left_seq])
                left_homotm = float(str(left_homodimer).split("tm=")[-1].split(", ")[0])
                right_homodimer = primer3.calcHomodimer(design_candidate[right_seq])
                right_homotm = float(str(right_homodimer).split("tm=")[-1].split(", ")[0])
                left_hairpin = primer3.calcHairpin(design_candidate[left_seq])
                left_hairtm = float(str(left_hairpin).split("tm=")[-1].split(", ")[0])
                right_hairpin = primer3.calcHairpin(design_candidate[right_seq])
                right_hairtm = float(str(right_hairpin).split("tm=")[-1].split(", ")[0])
                heterodimer = primer3.calcHeterodimer(design_candidate[left_seq], design_candidate[right_seq])
                heterotm = float(str(heterodimer).split("tm=")[-1].split(", ")[0])
                #
                if (recognize == "ref" and len(ref_cut_fragment_list) == 2) or (recognize == "alt" and len(ref_cut_fragment_list) == 1):
                    if left_homotm <= 45 and right_homotm <= 45 and left_hairtm <= 45 and right_hairtm <= 45 and heterotm <= 45:
                        cut_left = primer_range - left_start + 1
                        cut_right = right_end - primer_range
                        left_primer_coord = vpos - cut_left
                        right_primer_coord = vpos + cut_right
                        #
                        if cut_left >= 100 and cut_right >= 100:
                            fw.write(f"\tLEFT_p{i} : {design_candidate[left_pos]}\n")
                            fw.write(f"\tRIGHT_p{i} : {design_candidate[right_pos]}\n")
                            fw.write(f"\tLEFT_genomep{i} : {left_primer_coord}\n")
                            fw.write(f"\tRIGHT_genomep{i} : {right_primer_coord}\n")
                            fw.write(f"\tCut_{i} : {cut_left}, {cut_right}\n")
                            fw.write(f"\tLEFT_tm{i} : {round(design_candidate[left_tm],1)}\n")
                            fw.write(f"\tRIGHT_tm{i} : {round(design_candidate[right_tm],1)}\n")
                            fw.write(f"\tPRODUCT_{i} : {design_candidate[product_size]}\n")
                            #
                            lhomo_tm = f"\tleft_homodimer{i} : " + str(left_homodimer)
                            fw.write(f"{lhomo_tm}\n")
                            rhomo_tm = f"\tright_homodimer{i} : " + str(right_homodimer)
                            fw.write(f"{rhomo_tm}\n")
                            #
                            lhairp_tm = f"\tleft_hairpin{i} : " + str(left_hairpin)
                            fw.write(f"{lhairp_tm}\n")
                            rhairp_tm = f"\tright_hairpin{i} : " + str(right_hairpin)
                            fw.write(f"{rhairp_tm}\n")
                            #
                            hetero_tm = str(heterodimer)
                            fw.write(f"\theterodimer{i} : {hetero_tm}\n")
                        else:
                            fw.write(f"*** The fragment size after digestion is small (<120 bp)! {cut_left}, {cut_right}\n")
                    else:
                        fw.write("*** melting temperature is high!\n")
                        fw.write(f"*** {left_homotm}\t{right_homotm}\t{left_hairtm}\t{right_hairtm}\t{heterotm}\n")
                else:
                    fw.write(f"*** There are more than one recognizable sites!\n")
    fw.close()

def primer_design(fasta, label, size):
    '''return paired primer design'''
    fa_length = len(fasta)
    primer_return = primer3.bindings.designPrimers({"SEQUENCE_ID":label, "SEQUENCE_TEMPLATE":fasta, "SEQUENCE_INCLUDED_REGION":[0, fa_length], "SEQUENCE_TARGET":[size, 1]}, {"PRIMER_PICK_LEFT_PRIMER":1, "PRIMER_PICK_INTERNAL_OLIGO":0, "PRIMER_PICK_RIGHT_PRIMER":1, 'PRIMER_OPT_SIZE': 22, 'PRIMER_MIN_SIZE': 19, 'PRIMER_MAX_SIZE': 24, 'PRIMER_OPT_TM': 60.5, 'PRIMER_MIN_TM': 59.5, 'PRIMER_MAX_TM': 62.5, 'PRIMER_MIN_GC': 45.0, 'PRIMER_MAX_GC': 65.0, "PRIMER_NUM_RETURN":3, "PRIMER_PRODUCT_SIZE_RANGE":[[380, 450], [450, 550], [550, 650], [650,800], [800, 900], [900, fa_length]]})
    return primer_return

def main():
    args = args_parser()
    primer_info(args)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
