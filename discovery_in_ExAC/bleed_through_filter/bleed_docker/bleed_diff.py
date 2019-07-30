import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

BLEED_CONTEXT=["TG_C_T","TA_C_T","TC_C_T","TT_C_T","TG_C_G",
                    "CG_C_G",
                    "AT_C_T","AG_C_A","AA_C_A","AC_C_A","AT_C_A","AG_C_G",
                    "GT_C_T","GG_C_G","GA_C_G","GC_C_G","GT_C_G",
                    "TG_A_G","TC_A_C","TG_A_T","TA_A_T","TC_A_T","TT_A_T",
                    "CG_A_G","CG_A_C","CA_A_C","CC_A_C","CT_A_C","CT_A_T",
                    "AG_A_G","AC_A_C","AT_A_T",
                    "GG_A_G","GA_A_G","GC_A_G","GT_A_G","GC_A_C","GT_A_T",
                    "TG_G_T","TA_G_T","TC_G_T","TT_G_T","TA_G_A","TC_G_C",
                    "CT_G_T","CA_G_A","CG_G_C","CA_G_C","CC_G_C","CT_G_C",
                    "AT_G_T","AG_G_A","AA_G_A","AC_G_A","AT_G_A","AC_G_C",
                    "GT_G_T","GA_G_A","GC_G_C",
                    "TG_T_G","TC_T_C","TA_T_A",
                    "CG_T_G","CG_T_C","CA_T_C","CC_T_C","CT_T_C","CA_T_A",
                    "AG_T_G","AC_T_C","AG_T_A","AA_T_A","AC_T_A","AT_T_A",
                    "GG_T_G","GA_T_G","GC_T_G","GT_T_G","GC_T_C","GA_T_A"]

def extract_context(raw_cont_str):
    raw_cont_str = raw_cont_str.upper()
    cont_len = len(raw_cont_str)
    cont_len_mid = cont_len // 2

    return raw_cont_str[cont_len_mid-1] + raw_cont_str[cont_len_mid+1]

def main():

    input_filename=str(sys.argv[1])
    output_filename=str(sys.argv[2])
    diff_cutoff=int(sys.argv[3])

    bleed_cs=pd.read_csv(input_filename,sep='\t',comment="#",dtype=object)
    
    bleed_cs['cont_ext']=bleed_cs['context'].apply(extract_context)
    bleed_cs['cont_ext']=bleed_cs['cont_ext']+"_" + bleed_cs["ref_allele"].str.upper() + "_" + bleed_cs["alt_allele"].str.upper()
    
    bleed_cs["position"]=pd.to_numeric(bleed_cs["position"])
    index=~bleed_cs.cont_ext.isin(BLEED_CONTEXT)
    null_model=bleed_cs[index]
    null_model=null_model.loc[null_model['judgement'].str.match("KEEP"),:]
    null_model=null_model.sort_values(by=["cont_ext",'contig',"position"])
    
    null_diff=null_model['position'].diff()
    null_diff_next=null_model['position'].diff(periods=-1)
    null_model['diff_prev']=null_diff
    null_model['diff_next']=null_diff_next    
    
    null_diff=null_diff[null_diff>0]

    bleed_cs=pd.read_csv(input_filename,sep='\t',comment="#",dtype=object)
    bleed_cs['cont_ext']=bleed_cs['context'].apply(extract_context)
    bleed_cs['cont_ext']=bleed_cs['cont_ext']+"_" + bleed_cs["ref_allele"].str.upper() + "_" + bleed_cs["alt_allele"].str.upper()
    
    bleed_cs["position"]=pd.to_numeric(bleed_cs["position"])
    context_group=bleed_cs.loc[bleed_cs['cont_ext'].isin(BLEED_CONTEXT),:]
    context_group=context_group.sort_values(by=["cont_ext",'contig',"position"])
    
    diff_arr=context_group['position'].diff()
    diff_arr_next=context_group['position'].diff(periods=-1)

    diff_arr=diff_arr[diff_arr>0]
    
    x=range(len(diff_arr))
    
    plt.hist(np.log10(diff_arr), edgecolor='black', linewidth=1,bins=100)
    axes = plt.gca()

    plt.savefig(output_filename+"_diff_dist.png")
    plt.gcf().clear()
    
    context_group['diff_prev']=diff_arr
    context_group['diff_next']=diff_arr_next
    context_filtered=context_group[(context_group['diff_prev']>diff_cutoff)|(context_group['diff_next']>diff_cutoff)]
    
    context_filtered=pd.concat([null_model,context_filtered[context_filtered['judgement'].str.contains("KEEP")]])
    context_filtered=context_filtered.sort_values(by=["cont_ext",'contig',"position"])
    
    context_filtered.to_csv(output_filename+"_read_position_filtered.maf",sep='\t',index=False)

    writer=open("num_passed.txt",'w')
    num_pass=len(context_filtered.index)
    writer.write(str(num_pass))
    writer.close()

    writer=open("num_rejected.txt",'w')
    num_rejected=len(context_group.index)-len(context_filtered.index)
    writer.write(str(num_rejected))
    writer.close()
    
    x=range(len(null_diff))
    
    plt.hist(np.log10(null_diff), edgecolor='black', linewidth=1,bins=100)
    axes = plt.gca()
    #axes.set_xlim([0,750000])
    #axes.set_ylim([0,500])
    plt.savefig(output_filename+"_null_dist.png")
    plt.gcf().clear()
    
    
if __name__ == "__main__":
    main()