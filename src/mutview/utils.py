#-*- coding:utf-8 -*-
#!/usr/bin/env python

import subprocess

def color_seq(seq):

    # the base of color
    D_base_rule = {
                    "A" : "\033[1;31m", # red
                    "T" : "\033[1;32m", # green
                    "C" : "\033[1;33m", # yellow 
                    "G" : "\033[1;34m", # blue
                    "I" : "\033[1;35m",
                    "-" : "\033[1;36m"
                    }
    
    # the defaul terminal color
    tail_flag = "\033[0m"
    
    color_base_list = []
    for base in seq:
        color_base_list.append("{0}{1}{2}".format(
                                        D_base_rule.get(base, ""),
                                        base,
                                        tail_flag
                               ))
    return "".join(color_base_list)

def generate_color_html_seq(seq, return_type):
    
    D_color_rule = {
            "A" : "#FF0000", #red
            "a" : "#FF0000", #red
            "T" : "#00FF00", #green
            "t" : "#00FF00", #green
            "C" : "#FFFF00", #yellow
            "c" : "#FFFF00", #yellow
            "G" : "#00CCCC", #blue
            "g" : "#00CCCC", #blue
            "I" : "#FF00FF", #purple
            "i" : "#FF00FF",
            "_" : "#00FFFF", #cyan
            "N" : "#C0C0C0",  #gray
            "n" : "#C0C0C0",
            "|" : "#FF9900",
            }
    color_base_list = []
    for base in seq:
        color_base_list.append('<span style="color:{0};font-weight: bold">{1}</span>'.format(
                                                                              D_color_rule.get(base, "white"),
                                                                               base))
        #color_base_list.append('<span style="color:{0};">{1}</span>'.format(
        #                                                                       D_color_rule.get(base, "white"),
        #                                                                       base))

    if return_type == "l":
        return color_base_list
    else:
        return "".join(color_base_list)
    
def get_genome_seq(chrom, start_pos, end_pos, ref, samtools_path):
    samtools_info = subprocess.Popen('{4} faidx {0} {1}:{2}-{3}'.format(
                                    ref, chrom, start_pos, end_pos, samtools_path),
                                    shell = True,
                                    stdout = subprocess.PIPE,
                                    stderr = subprocess.PIPE
                                    )
    stdout_value, stderr_value = samtools_info.communicate()    
    if stderr_value:
        raise(Exception("samtools faidx error!\n" + stderr_value))
    else:
        return  "".join(str(stdout_value, encoding="utf8").split("\n")[1:])
        




