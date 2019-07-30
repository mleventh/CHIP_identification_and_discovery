import pandas as pd 
import numpy as np
import pickle
import sys
from scipy import stats
import statsmodels.stats.multitest as multi

def remove_exac_sites_from_maf(call_stats_table, exac_file):
    # use ExAC vcf to filter likely germline variants out of candidate sites
    with open(exac_file, 'rb') as handle:
        exac_dict = pickle.load(handle)
        keep = np.ones_like(call_stats_table['Start_position'], dtype=bool)
    for index, row in call_stats_table.iterrows():
        key = str(row['Chromosome']) + '_' + str(row['Start_position'])
        try:
            exac_dict[key]
            #print('removing site '+ key+ ' minor allele fraction = ' + str(exac_dict[key]['AF']))
            keep[index] = False
        except KeyError:
            pass
        except IndexError:
            pass
    return call_stats_table[keep]

def main():
    input_file=sys.argv[1]
    exac_file=sys.argv[2]
    output_filename=sys.argv[3]

    maf=pd.read_csv(input_file,sep='\t')
    afb=np.arange(0,1.01,0.01)

    # pox=pd.to_numeric(maf.pox)>pd.to_numeric(maf.pox_cutoff)
    # GT=maf.Reference_Allele.str.match("G")&maf.Tumor_Seq_Allele2.str.match("T")
    # CA=maf.Reference_Allele.str.match("C")&maf.Tumor_Seq_Allele2.str.match("A")

    # maf=maf[~((GT&pox)|(CA&pox))]

    maf.dropna(subset=["i_t_ALT_F1R2","i_t_ALT_F2R1"],inplace=True)

    maf['orientation_bias_F1']=stats.binom.cdf(pd.to_numeric(maf.i_t_ALT_F1R2),
                              pd.to_numeric(maf.i_t_ALT_F1R2)+pd.to_numeric(maf.i_t_ALT_F2R1),0.96)
    maf['orientation_bias_F2']=stats.binom.cdf(pd.to_numeric(maf.i_t_ALT_F2R1),
                                 pd.to_numeric(maf.i_t_ALT_F1R2)+pd.to_numeric(maf.i_t_ALT_F2R1),0.96)

    maf['orientation_bias_F1']=multi.multipletests(maf['orientation_bias_F1'],method='fdr_bh')[1]
    maf['orientation_bias_F2']=multi.multipletests(maf['orientation_bias_F2'],method='fdr_bh')[1]

    maf=maf.loc[(maf.orientation_bias_F1<0.99)&(maf.orientation_bias_F2<0.99),:]

    # maf=maf[pd.to_numeric(maf.t_lod_fstar)>8.6]
    # print(len(maf.index))
    maf=maf.reset_index(drop=True)

    #maf=maf[pd.to_numeric(maf.keep_max)>60]
    #maf=maf.reset_index(drop=True)

    maf=remove_exac_sites_from_maf(maf,exac_file)
    maf=maf.reset_index(drop=True)

    beta_vals=[]

    for idx,row in maf.iterrows():
        distribution=stats.beta.pdf(afb,row.t_alt_count+1,row.t_ref_count+1)
        beta_vals.append(distribution)

    beta_df=pd.DataFrame(beta_vals,columns=afb)

    maf=maf.loc[:,("Hugo_Symbol","Chromosome","Start_position",
              "Variant_Classification", "Reference_Allele","Tumor_Seq_Allele2",
              "Protein_Change","Tumor_Sample_Barcode")]

    formatted_maf=pd.concat([maf,beta_df],axis=1)

    formatted_maf.to_csv(output_filename,sep='\t',index=False)

if __name__ == '__main__':
    main()