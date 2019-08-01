import pandas as pd 
import numpy as np
import pickle
import sys
from scipy import stats
import statsmodels.stats.multitest as multi

def remove_exac_sites_from_vcf(call_stats_table, exac_file):
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
    sample_name=sys.argv[4]

    vcf=pd.read_csv(input_file,sep='\t',compression='gzip',comment="#")
    afb=np.arange(0,1.01,0.01)
    
    #re-format VCF as a phylogic-ready MAF
    vcf.columns=['Chromosome', 'Start_position',"ID",'Reference_Allele','Tumor_Seq_Allele2','QUAL','FILTER',
                'INFO','FORMAT',sample_name]

    #eliminate chr prefixes
    vcf['Chromosome'].replace(to_replace=r'^chr',value="",regex=True,inplace=True)

    #filter rejected sites
    vcf=vcf[vcf['FILTER'].str.match('PASS')]
    vcf=vcf.reset_index(drop=True)

    #eliminate sites with >1% population allele frequency in ExAC
    vcf=remove_exac_sites_from_vcf(vcf,exac_file)
    vcf=vcf.reset_index(drop=True)

    count_info=vcf[sample_name].str.split(':',expand=True)

    vcf['t_alt_count']=pd.to_numeric(count_info.iloc[:,1].str.split(',',expand=True)[1])
    vcf['t_ref_count']=pd.to_numeric(count_info.iloc[:,1].str.split(',',expand=True)[0])

    vcf['t_ALT_F1R2']=pd.to_numeric(count_info.iloc[:,3].str.split(',',expand=True)[1])
    vcf['t_ALT_F2R1']=pd.to_numeric(count_info.iloc[:,4].str.split(',',expand=True)[1])

    #eliminate orientation bias artifacts and variants with small LOD relative of a mutation rate of .2/Mb
    vcf['orientation_bias_F1']=stats.binom.cdf(pd.to_numeric(vcf.t_ALT_F1R2),
                              pd.to_numeric(vcf.t_ALT_F1R2)+pd.to_numeric(vcf.t_ALT_F2R1),0.96)
    vcf['orientation_bias_F2']=stats.binom.cdf(pd.to_numeric(vcf.t_ALT_F2R1),
                                 pd.to_numeric(vcf.t_ALT_F1R2)+pd.to_numeric(vcf.t_ALT_F2R1),0.96)

    vcf['orientation_bias_F1']=multi.multipletests(vcf['orientation_bias_F1'],method='fdr_bh')[1]
    vcf['orientation_bias_F2']=multi.multipletests(vcf['orientation_bias_F2'],method='fdr_bh')[1]

    vcf=vcf.loc[(vcf.orientation_bias_F1<0.99)&(vcf.orientation_bias_F2<0.99),:]

    vcf=vcf.reset_index(drop=True)

    vcf['t_lod_fstar']=vcf['INFO'].str.extract(r'(TLOD=[0-9]+\.[0-9]+)')
    vcf['t_lod_fstar'].replace(to_replace=r'TLOD=',value="",regex=True,inplace=True)

    vcf=vcf[pd.to_numeric(vcf['t_lod_fstar'])>8]
    vcf=vcf.reset_index(drop=True)

    vcf['Variant_Type']="SNP"

    vcf.loc[(vcf['Tumor_Seq_Allele2'].str.len()==1)&(vcf['Reference_Allele'].str.len()>1),'Variant_Type']="DEL"
    vcf.loc[(vcf['Tumor_Seq_Allele2'].str.len()>1)&(vcf['Reference_Allele'].str.len()==1),'Variant_Type']="INS"
    vcf['Tumor_Seq_Allele1']=vcf['Tumor_Seq_Allele2']
    vcf_ind=vcf[~vcf['Variant_Type'].str.match("SNP")]
    vcf_ind.reset_index(drop=True)
    vcf=vcf[vcf['Variant_Type'].str.match("SNP")]
    vcf=vcf.reset_index(drop=True)

    beta_vals=[]

    #calculate beta pdf for each variant in the defined allele fraction bins
    for idx,row in vcf.iterrows():
        distribution=stats.beta.pdf(afb,row['t_alt_count']+1,row['t_ref_count']+1)
        beta_vals.append(distribution)

    beta_df=pd.DataFrame(beta_vals,columns=afb)

    vcf=vcf.loc[:,("Chromosome","Start_position",
              "Variant_Classification", "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Variant_Type",
              "Tumor_Sample_Barcode")]

    #append beta pdf information
    formatted_maf=pd.concat([vcf,beta_df],axis=1)

    formatted_maf.to_csv(output_filename,sep='\t',index=False)

    vcf_ind.to_csv(output_filename+".indel",sep='\t',index=False)

if __name__ == '__main__':
    main()
