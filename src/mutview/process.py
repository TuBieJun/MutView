#-*- coding:utf-8 -*-
# author: liteng
# email: 707704459@qq.com
#!/usr/bin/env python3


import re
import sys
import pysam
import argparse
from . import utils

def out_html_header(title_name, font, handler):
    '''
    this function  we can use jinja2 to instead but that will be slow 
    and a little consumption of memeroy.
    If in the future html code becomes complicated, you can consider using jinja2.
    '''

    handler.write('<style type="text/css">\n')
    handler.write('table {border:1px solid gray;padding:2x;border-collapse:collapse;}\n')
    handler.write('td {border:1px solid gray;padding-left:2px;padding-right:2px;text-align:left}\n')
    handler.write('p {{color:white;font-family: {font}}}\n'.format(font=font))
    handler.write('</style>\n')
    handler.write('<title>{0}</title>\n'.format(title_name))
    handler.write('<body style="background-color:black;font-family: {font}">\n'.format(font=font))
    handler.write('<table>\n')

def out_html_tail(handler):
    handler.write('</table>\n')
    handler.write('</body>\n')

def insert_str(seq, pos_list, char):
    l = list(seq)
    for pos in pos_list:
        l.insert(pos, char)
    return l

def insert_list(l, pos_list, char):
    for pos in pos_list:
        l.insert(pos, char)

def out_sorted_seq(mut_base, D_record, max_n, handler):
    n = 0
    order_rule = ["A", "a", "T", "t", "C", "c", "G", "g", "I", "i", "_", "N", "n"]
    order_rule.remove(mut_base)
    order_rule.remove(mut_base.lower())
    for allele in order_rule:
        for show_seq in D_record.get(allele, []):
            handler.write(show_seq+"\n")
            n += 1
            if n >= max_n:
                return None
    for show_seq in D_record.get(mut_base, []):
        handler.write(show_seq+"\n")
        n += 1
        if n >  max_n:
            return None
    for show_seq in D_record.get(mut_base.lower, []):
        handler.write(show_seq+"\n")
        n += 1
        if n >  max_n:
            return None


def filter_base_qual(seq_l, qual_l, min_bq):
    return [ j if qual_l[i] >= min_bq else j.lower() for i,j in enumerate(seq_l) ]

def out_reads(sam_iter, chrom, var_pos, window_start, 
              window_end, genome_base_list, mut_base, show_genome_info, 
              min_bq, min_aq, handler, max_record_num=100, show_mut=False):
    
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
                #sub_list = [ utils.color_seq(base) for base 
                #           in raw_base_list[index:index+length] ]
                sub_list = filter_base_qual(record.seq[index:index+length],
                                            record.qual[index:index+length],
                                            min_bq)
                sub_list = utils.generate_color_html_seq(sub_list, 
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
                    utils.generate_color_html_seq("|", "s"))

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
        for k in range(len(new_base_list)):
            v = new_base_list[k]
            if ("<" not in v) and ( v.upper() != genome_base_list_cut[k])  \
                and (v.upper() != "N") and ("|" not in v):
                new_base_list[k] = utils.generate_color_html_seq(v, "s")
        if s_aln:
            new_base_list = ["&nbsp;"]*(show_start - window_start) + new_base_list

        # out the result
        show_seq = "".join(new_base_list)
        if record.mapping_quality >= min_aq:
            seq_info_html = '<td><p>{0}</p></td>'.format(show_seq)
        else:
            seq_info_html = '<td bgcolor="#B9B9FF"><p>{0}</p></td>'.format(show_seq)
        show_info = '<tr>\n{0}\n{1}</tr>'.format(read_info_html, seq_info_html)
        if show_mut:
            D_record.setdefault(record_allele, []).append(show_info)
        else:
            #out show_info
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

    handler.write(show_depth_info+"\n")
    handler.write(show_genome_info+"\n")
    if show_mut:
        out_sorted_seq(mut_base, D_record, max_record_num, handler)
    else:
        for i in show_record_list:
            handler.write(i+"\n")
        
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pos", help="the mutation position", type=int, required=True)
    parser.add_argument("-c", "--chrom", help="the chrom of mutation", required=True)
    parser.add_argument("-b", "--bam", help="the input bam file",required=True)
    parser.add_argument("-r", "--ref", help="the reference genome fasta file",required=True)
    parser.add_argument("-o", "--out", help="the output html file",required=True)
    parser.add_argument("-w", "--window_size", help="the window size of reference genome around the mutation, default is 300",
                        type=int, default=300)
    parser.add_argument("-s", "--sort" , action="store_true", default=False,
                        help="sort as base like igv, default is False")
    parser.add_argument("-n", "--num_reads", type=int, default=100,
                        help="the max reads number to view default is 100")
    parser.add_argument("-q", "--base_q", type=int, default=20, 
                        help="the min base quality, reads with bq below the specified threshold will be lowercases. default is 20.")
    parser.add_argument("-Q", "--aln_q", type=int, default=20,
                        help="the min aln quality, bases with mq below the specified threshold will be highlighted. default is 20.")
    parser.add_argument("-f", "--font", type=str, default="consolas", 
                        help="the the font of the out html file, it should be monospaced font, default is consolas for windows, \
                            if you view the out html file in mac, you need set the value to Menlo or other monospaced fonts")
    parser.add_argument("-S", "--samtools", default="samtools",
                        help="the samtools path, default is samtools")

    if len(sys.argv) == 1:
        parser.parse_args(["-h"])
    args = parser.parse_args()
        

    
    min_bq = chr(args.base_q + 33)
    min_aq = args.aln_q
    user_pos = args.pos
    chrom = args.chrom
    referece_fasta = args.ref
    bam_file = args.bam
    out_file = args.out
    font = args.font
    samtools_path = args.samtools
    with open(out_file, "w") as out_handler:
        window_size = int(args.window_size/2)
    
        window_start = user_pos - window_size
        window_end = user_pos +  window_size

        # out the html first half
        out_html_header("{0}:{1}".format(chrom, user_pos), font, out_handler)

        # out the reference seq 
        genome_base_list = list(utils.get_genome_seq(chrom, window_start,
                                            window_end, referece_fasta, samtools_path))
        insert_list(genome_base_list, [window_size+1, window_size], "|")

        genome_base_seq_colored = utils.generate_color_html_seq(genome_base_list, "s")
        show_genome_info = '<tr>       \
                               <td><p>{0}</p></td> \
                               <td><p>{1}</p></td> \
                            </tr>'.format(
                                        "Reference",
                                        genome_base_seq_colored
                                        )

        # out the fusion reads cover user pos
        sam_file = pysam.AlignmentFile(bam_file, "rb")
        sam_iter = sam_file.fetch(chrom, user_pos-150, user_pos+150)
        mut_base = utils.get_genome_seq(chrom, user_pos, user_pos+1, referece_fasta, samtools_path)[0]
        if args.sort:
            out_reads(sam_iter, chrom, user_pos, window_start, window_end, 
                      genome_base_list, mut_base, show_genome_info,
                      min_bq, min_aq, out_handler, max_record_num = args.num_reads, show_mut = True)
        else:
            out_reads(sam_iter, chrom, user_pos, window_start, window_end, 
                      genome_base_list, mut_base, show_genome_info,
                      min_bq, min_aq, out_handler, max_record_num = args.num_reads,show_mut = False)

        out_html_tail(out_handler)

if __name__ == "__main__":
    main()
                                           
