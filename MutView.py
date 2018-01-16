#-*- coding:utf-8 -*-
# author: liteng
# email: 707704459@qq.com
#!/usr/bin/env python


import re
import sys
import pysam
import argparse
import dear_seq

def print_html_header(title_name):
    '''
    this function  we can use jinja2 to instead but that will be slow 
    and a little consumption of memeroy.
    If in the future html code becomes complicated, you can consider using jinja2.
    '''
    #print '<header>'
    #print '<meta charset="utf-8">'
    print '<style type="text/css">'
    print 'table {border:1px solid gray;padding:2x;border-collapse:collapse;}'
    print 'td {border:1px solid gray;padding-left:2px;padding-right:2px;text-align:left}'
    print 'p {color:white;font-family: consolas}'
    print '</style>'
    print '<title>{0}</title>'.format(title_name)
    #print '</header>'
    print '<body style="background-color:black;font-family: consolas">'
    print '<table>'

def print_html_tail():
    print '</table>'
    print '</body>'

def insert_str(seq, pos_list, char):
    l = list(seq)
    for pos in pos_list:
        l.insert(pos, char)
    return l

def insert_list(l, pos_list, char):
    for pos in pos_list:
        l.insert(pos, char)

def show_sorted_seq(mut_base, D_record, max_n):
    n = 0
    order_rule = ["A", "a", "T", "t", "C", "c", "G", "g", "I", "i", "_", "N", "n"]
    order_rule.remove(mut_base)
    order_rule.remove(mut_base.lower())
    for allele in order_rule:
        for show_seq in D_record.get(allele, []):
            print show_seq
            n += 1
            if n >= max_n:
                return None
    for show_seq in D_record.get(mut_base, []):
        print show_seq
        n += 1
        if n >  max_n:
            return None
    for show_seq in D_record.get(mut_base.lower, []):
        print show_seq
        n += 1
        if n >  max_n:
            return None


def filter_base_qual(seq_l, qual_l, min_bq):
    return [ j if qual_l[i] >= min_bq else j.lower() for i,j in enumerate(seq_l) ]

def show_reads(sam_iter, chrom, var_pos, window_start, 
              window_end, genome_base_list, mut_base, show_genome_info, 
              min_bq,min_aq,max_record_num=100, show_mut=False):
    
    D_seq = {}
    D_record = {}
    D_base_count = {}
    depth_count = 0
    show_record_num = 0
    show_record_list = []
    for record in sam_iter:
        if show_record_num > max_record_num:
            break

        # the pysam alin pos is advance 1 bp
        spos = record.reference_start + 1
        show_start = spos
        show_end = spos - 1

        # get the read info 
        read_name = record.qname
        if record.is_read1:
            pair_info = "R1"
        elif record.is_read2:
            pair_info = "R2"
        else:
            pair_info = "*"
        if record.is_reverse:
            strand_info = "reverse"
        else:
            strand_info = "positive"

        read_info_html = "<td><p>{0}_{1}_{2}</p></td>".format(read_name,
                                                              pair_info,
                                                              strand_info)
        #raw_base_list = list(record.seq)
        #raw_qual_list = list(record.qual)
        new_base_list = []

        if show_start > var_pos:
            continue
        index = 0
        match_l = 0
        match_or_mismatch_l = 0
        md_index = 0



        if not record.cigartuples:
            continue
        
        for i, (flag, length) in enumerate(record.cigartuples):
            # softclip
            sub_list = []
            if flag == 4:
                #sub_list = [ dear_seq.color_seq(base) for base 
                #           in raw_base_list[index:index+length] ]
                sub_list = filter_base_qual(record.seq[index:index+length],
                                            record.qual[index:index+length],
                                            min_bq)
                sub_list = dear_seq.generate_color_html_seq(sub_list, 
                                                            "l")
                index += length
                if i == 0:
                    show_start = spos - length
                else:
                    show_end += length
            # match or mismatch
            if flag == 0:
                #sub_list = raw_base_list[index:index+length]
                sub_list = filter_base_qual(record.seq[index:index+length],
                                            record.qual[index:index+length],
                                            min_bq)
                index += length
                show_end += length
            # deletion
            if flag == 2:
                sub_list = ["_"]*length
                show_end += length
            # insert
            if flag == 1:
                var_base = "I" if record.qual[index] >= min_bq else "i"
                index += length
                
                if new_base_list:
                    new_base_list[-1] = var_base
                else:
                    new_base_list.append(var_base)
                    
            new_base_list.extend(sub_list)

        # reomove the record that not cover user pos 
        if show_end < var_pos or show_start > var_pos:
            continue
        # find the base in user_pos
        depth_count += 1
        read_pos = var_pos - show_start
        temp = -1
        for i, j in record.cigartuples:
            if i != 2:
                temp += j
                if temp >= read_pos:
                    read_pos_cigar = i
                    break
        else:
            continue

        if ">" in new_base_list[read_pos]:
            single_mut = re.search(r'>([ATCGNIatcgni])<', new_base_list[read_pos])
            if single_mut:
                record_allele = single_mut.group(1)
            else:
                if re.search(r'>(\_+)<', new_base_list[read_pos]):
                    record_allele = "_"
                else:
                    record_allele = None

 
            #record_allele = re.search(r'>([ATCGN\_I])<', new_base_list[read_pos]).group(1)
        else:
            record_allele = new_base_list[read_pos]
        if read_pos_cigar != 4 and read_pos_cigar != 3:
            D_base_count[record_allele] = D_base_count.get(record_allele, 0) + 1
        
        insert_list(new_base_list, [read_pos+1, read_pos], 
                    dear_seq.generate_color_html_seq("|", "s"))

        # cut the show seq if the show seq is longer than window size
        if show_end > window_end:
            new_base_list = new_base_list[:window_end-show_end]
          
        if show_start >= window_start:
            s_aln = show_start - window_start
        else:
            s_aln = 0
            new_base_list= new_base_list[window_start - show_start:]

        genome_base_list_cut = genome_base_list[s_aln:] 

        # compare reference and record
        for k  in xrange(len(new_base_list)):
            v = new_base_list[k]
            if ("<" not in v) and ( v.upper() != genome_base_list_cut[k])  \
                and (v.upper() != "N") and ("|" not in v):
                new_base_list[k] = dear_seq.generate_color_html_seq(v, "s")
        if s_aln:
            new_base_list = ["&nbsp;"]*(show_start - window_start) + new_base_list

        # print the result
        #new_base_list.append("</br>") #B9B9FF
        show_seq = "".join(new_base_list)
        if record.mapping_quality >= min_aq:
            seq_info_html = '<td><p>{0}</p></td>'.format(show_seq)
        else:
            seq_info_html = '<td bgcolor="#B9B9FF"><p>{0}</p></td>'.format(show_seq)
        show_info = '<tr>\n{0}\n{1}</tr>'.format(read_info_html, seq_info_html)
        if show_mut:
            D_record.setdefault(record_allele, []).append(show_info)
        else:
            #print show_info
            #show_record_num += 1
            if len(show_record_list) < max_record_num:
                show_record_list.append(show_info)
    show_depth_info = '<tr> \
                            <td><p>{0}:{1}</p></td> \
                            <td><p>A:{2}&nbsp;T:{3}&nbsp;C:{4}&nbsp;G:{5}&nbsp; \
                            N:{6}&nbsp;ins:{7}&nbsp;del:{8}&nbsp;</p></td> \
                       </tr>'.format(chrom, var_pos, D_base_count.get("A", 0) + D_base_count.get("a", 0),
                                    D_base_count.get("T", 0) + D_base_count.get("t", 0) , D_base_count.get("C", 0) + D_base_count.get("c", 0), 
                                    D_base_count.get("G", 0) + D_base_count.get("g", 0), D_base_count.get("N", 0) + D_base_count.get("n", 0),
                                    D_base_count.get("I", 0) + D_base_count.get("i", 0) , D_base_count.get("_", 0))

    print show_depth_info
    print show_genome_info
    if show_mut:
        show_sorted_seq(mut_base, D_record, max_record_num)
    else:
        for i in show_record_list:
            print i
        
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pos", help="the fusion pos", type=int)
    parser.add_argument("-c", "--chrom", help="the chrom of fusion")
    parser.add_argument("-b", "--bam", help="the input bam file")
    parser.add_argument("-r", "--ref", help="the reference file, default is hg19",
                        default = "/biocluster/data/bioexec/database/genome/"
                        "Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
    parser.add_argument("-w", "--window_size", help="the terminal size default is 300",
                        type=int, default=300)
    parser.add_argument("-s", "--sort" , action="store_true", default=False,
                        help="sort as base like igv, default is False")
    parser.add_argument("-n", "--num_reads", type=int, default=100,
                        help="the max reads number to view default is 100")
    parser.add_argument("-q", "--base_q", type=int, default=20, help="the min base quality")
    parser.add_argument("-Q", "--aln_q", type=int, default=20, help="the min aln quality")
    args = parser.parse_args()
    
    min_bq = chr(args.base_q + 33)
    min_aq = args.aln_q
    user_pos = args.pos
    chrom = args.chrom
    referece_fasta = args.ref
    bam_file = args.bam
    window_size = args.window_size/2
    
    # test data
    #user_pos = 29447670
    #chrom = "2"
    #referece_fasta = "/biocluster/data/bioexec/database/genome/Hom
    #                 "o_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
    #window_size = args.window_size/2
    #bam_file = "/bionfsdate/ctDNA/experiment/wangmingxia/22.tumor_c" \
    #           "fDNA/FC80P-D1/second_aln_str/FC80P-D1/01_aln/FC80P-D1.x.0_realigned.bam"

    window_start = user_pos - window_size
    window_end = user_pos +  window_size

    #print the html first half
    print_html_header("chr{0}:{1}".format(chrom, user_pos))

    # show the reference seq 
    genome_base_list = list(dear_seq.get_genome_seq(chrom, window_start,
                                        window_end, referece_fasta))
    insert_list(genome_base_list, [window_size+1, window_size], "|")

    #genome_base_list_colored = [ dear_seq.color_seq(base) 
    #                           for base in genome_base_list]
    genome_base_seq_colored = dear_seq.generate_color_html_seq(genome_base_list, "s")
    #genome_base_seq_colored += "</br>"
    show_genome_info = '<tr>       \
                           <td><p>{0}</p></td> \
                           <td><p>{1}</p></td> \
                        </tr>'.format(
                                    "Reference",
                                    genome_base_seq_colored
                                    )
    #print show_genome_info

    # show the fusion reads cover user pos
    sam_file = pysam.AlignmentFile(bam_file, "rb")
    sam_iter = sam_file.fetch(chrom, user_pos-150, user_pos+150)
    mut_base = dear_seq.get_genome_seq(chrom, user_pos, user_pos+1, referece_fasta)[0]
    if args.sort:
        show_reads(sam_iter, chrom, user_pos, window_start, window_end, 
                  genome_base_list, mut_base, show_genome_info,
                  min_bq, min_aq, max_record_num = args.num_reads, show_mut = True)
    else:
        show_reads(sam_iter, chrom, user_pos, window_start, window_end, 
                  genome_base_list, mut_base, show_genome_info,
                  min_bq, min_aq, max_record_num = args.num_reads,show_mut = False)

    #print the html second half
    print_html_tail()


if __name__ == "__main__":
    main()
                                           
