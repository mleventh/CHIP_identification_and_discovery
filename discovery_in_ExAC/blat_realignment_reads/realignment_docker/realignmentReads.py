import subprocess,time,tempfile,sys,os,random,string,socket
import twobitreader
import blat_realign_cython

def chrom2int(x):
   if x == 'X':
      return 23
   elif x == 'Y':
      return 24
   else:
      return int(x)

def get_reference(genome,chrom,start,end):

    contig=genome[str(chrom)]
    return contig[start:end+1]

class Wrapper: #for running in FH
    @staticmethod
    def run_realign(maf,hg19="/opt/hg19.2bit",out_name=None,min_diff=5,cov_fraction=0.05,min_reads_cov_fraction=5):
        Realignment_Filter(maf,hg19,out_name=out_name,min_diff=min_diff,cov_fraction=cov_fraction,min_reads_cov_fraction=min_reads_cov_fraction)


class BLAT_server:
    def __init__(self,ref_2bit,port=random.randint(58142,65530),debug=False):
        self.server=blat_realign_cython.pyRealignmentFilter(ref_2bit,"/opt/hg19.11.ooc")


    def query(self,sequence,chrN,pos):
        return(self.server.blat_seq(sequence,chrN,pos))


class Realignment_Filter:
    @staticmethod
    def get_stats(score_diff_list):
        return map(str, [min(score_diff_list), sum(score_diff_list) / float(len(score_diff_list)), max(score_diff_list)]) if len(score_diff_list) else ["", "", ""]

    def __init__(self,maf,out_name=None,hg19="/opt/hg19.2bit",min_diff=5,cov_fraction=0.05,min_reads_cov_fraction=5):
        self.input_maf = open(maf,"r")
        self.blat_server = BLAT_server(hg19)
        self.min_diff=float(min_diff)
        self.cov_fraction=float(cov_fraction)
        self.min_reads_cov_fraction=int(min_reads_cov_fraction)

        if not out_name:
            out_name=maf.rpartition("/")[-1].rpartition(".maf")[0]
        self.run_filter(out_name)


    def blat_read(self,read,mut):
        score=self.blat_server.query(read,mut[0],mut[1])
        try:

            if score > self.min_diff:
                return [True,score]
            else:
                return [False,score] #If it's a tie, it's inconclusive, but the mate should decide. If both have double then it's questoinable.

        except:
            print "Can't Map! or blat failure. -",
            print read
            return [-1,0]

    def find_reads(self,chrm,pos,alt,ref,ac,depth):
        concordant=0
        discordant=0
        mate_not_discordant=0
        refrence=0
        blat_fail=0
        read_names = set()
        with_minus_without=[]
        reads=[]
        stream=range(10,66,5)

        # initialize 2bitreader genome
        genome = twobitreader.TwoBitFile("/opt/hg19.2bit")
        chr=[genome[c] for c in map(str,range(1,23))]
        chr.append(genome['X'])
        chr.append(genome['Y'])

        for index in stream:

            if alt != "-" and ref != "-":
                read=get_reference(genome,chrm,pos-index,pos+(75-index))
                read=read[0:index]+alt+read[index+1:len(read)]
                reads.append(read)
            elif alt != "-" and ref == "-":
                for point in range(1,len(alt)+1):
                    read_back=get_reference(genome,chrm,pos-index-point,pos+(75-index)-point)

                    if index+point+1+len(alt)<=76:
                        alt_read=read_back[0:index+point+1]+alt+read_back[index+point+1:len(read_back)-len(alt)]
                        reads.append(alt_read)

                    read_forward=get_reference(genome,chrm,pos-index+point,pos+(75-index)+point)

                    if index+1+len(alt)<=76:
                        alt_read=read_forward[0:index+1]+alt+read_forward[index+1:len(read_forward)-len(alt)]                    
                        reads.append(alt_read)

            elif ref != "-" and alt == "-":
                for point in range(1,len(ref)+1):
                    read_back=get_reference(genome,chrm,pos-index-point,pos+(75-index)+len(ref)-point)
                    if index+point<=76:
                        alt_back=read_back[0:index+point]+read_back[index+point+len(ref):(75+index+point+len(ref))]
                        reads.append(alt_back)

                    read_forward=get_reference(genome,chrm,pos-index,pos+(75-index)+len(ref)+point)

                    if point<=index:
                        alt_forward=read_forward[point:index]+read_forward[index+len(ref):(75+index+len(ref)+point)]
                        reads.append(alt_forward)

                    
                    


        for read in reads:

            is_alt=True
            #if alt != "-" and ref != "-": #indels handled below
                #is_alt=True if read[stream] == alt else False

            #if read.query_position or is_alt:
            #read_names.add(read.alignment.qname)

            if is_alt:
                #read_names.add(read.alignment.qname)
                blat_res=self.blat_read(read,[chrm,pos])
                if blat_res[0] is True: #if returns true
                    with_minus_without.append(blat_res[1])

                    if concordant < 3:
                        concordant +=1
                    elif discordant < 1: #If we matched three reads, and didn't map one away, stop checking.
                        break

                elif blat_res[0] == -1:
                    blat_fail+=1
                    continue

                else: #if returns false
                    with_minus_without.append(blat_res[1])
                    discordant +=1
            else:
                refrence+=1

            if concordant+discordant+mate_not_discordant+blat_fail >= ac:
                break
            if discordant > max([5,depth]):
                break

        keep_scores=[]
        tie_scores=[]
        reject_scores=[]
        for score_diff in with_minus_without:
            if score_diff > self.min_diff:
                keep_scores.append(score_diff)
            elif score_diff == self.min_diff:
                tie_scores.append(self.min_diff)
            elif score_diff < self.min_diff:
                reject_scores.append(score_diff)

        return [concordant,discordant,mate_not_discordant,refrence,blat_fail]+list(Realignment_Filter.get_stats(keep_scores))+[str(len(tie_scores))]+list(Realignment_Filter.get_stats(reject_scores))+[",".join(map(str,keep_scores)),",".join(map(str,reject_scores))]
    def run_filter(self,out_name="no_id"):
        self.input_maf.seek(0) #make sure we're at the start
        with open(out_name+".blat.all.maf","w") as new_maf,open(out_name+".blat.maf","w") as pass_maf, open(out_name+".blat.rejected.maf","w") as debug_maf:
            header=self.input_maf.readline()
            while header[0] == "#" or not header.strip():
                #new_maf.write(header)
                header=self.input_maf.readline()

            header=header.strip("\n").split("\t")
            new_maf.write("\t".join(header+["reads_confirmed","reads_mapped_away","mate_not_mapped_away","reference","bad_map","realign_judgment"]+["keep_min","keep_mean","keep_max","n_tie","reject_min","reject_mean","reject_max","keep_scores","reject_scores"])+"\n")
            pass_maf.write("\t".join(header+["reads_confirmed","reads_mapped_away","mate_not_mapped_away","reference","bad_map","realign_judgment"]+["keep_min","keep_mean","keep_max","n_tie","reject_min","reject_mean","reject_max","keep_scores","reject_scores"])+"\n")
            debug_maf.write("\t".join(header+["reads_confirmed","reads_mapped_away","mate_not_mapped_away","reference","bad_map","realign_judgment"]+["keep_min","keep_mean","keep_max","n_tie","reject_min","reject_mean","reject_max","keep_scores","reject_scores"])+"\n")
            try: #maflite
                contig_i=header.index("chr")
                pos_i=header.index("start")
                ref_allele_i=header.index("ref_allele")
                alt_allele_i=header.index("tum_allele2")


            except: #maf
                contig_i=header.index("Chromosome")
                pos_i=header.index("Start_position")
                ref_allele_i=header.index("Reference_Allele")
                alt_allele_i=header.index("Tumor_Seq_Allele2")

            try:
                            t_ref_count=header.index("t_ref_count")
                            t_alt_count=header.index("t_alt_count")
            except:#abs maf
                t_ref_count=header.index("ref_cnt")
                t_alt_count=header.index("alt_cnt")

            start = time.time()/60.0
            last_chk = start
            for idx,line in enumerate(self.input_maf):
                vals=line.strip("\n").split("\t")
                try:
                    ref_count=int(float(vals[t_ref_count]))
                    alt_count=int(float(vals[t_alt_count]))#PCAWG mafs have alt/ref counts as 6.0 for some reason.
                except:
                    print "Detected mutation with no alt count, skipping!"
                    continue
                contig_for_blat = "MT" if vals[contig_i] == "M" else vals[contig_i] ##fix inconsistency 
                out_stats = self.find_reads(contig_for_blat,int(vals[pos_i]),vals[alt_allele_i],vals[ref_allele_i],alt_count,alt_count+ref_count)
                num_rejected=int(out_stats[1])
                num_confirmed=int(out_stats[0])+int(out_stats[2])


                if num_rejected and (((alt_count-num_rejected < 3) or num_confirmed == 0)  or (num_rejected >= max([self.cov_fraction * (alt_count+ref_count),self.min_reads_cov_fraction]))):
                    judj="REJECT"
                    debug_maf.write("\t".join(vals+map(str,out_stats[0:5])+[judj]+out_stats[5:])+"\n")

                else:
                    judj="KEEP"
                    pass_maf.write("\t".join(vals+map(str,out_stats[0:5])+[judj]+out_stats[5:])+"\n")


                new_maf.write("\t".join(vals+map(str,out_stats[0:5])+[judj]+out_stats[5:])+"\n")
                if idx % 100 == 0:
                    curr_time=time.time()/60.0
                    print "read {} sites in {}min. Total sites {}, total time {}hrs".format(100,curr_time-last_chk,idx,(curr_time-start)/60.0)
                    last_chk=curr_time
            end=time.time()/60.0
        print "total time ",end-start," minutes"


#RF=Realignment_Filter(bam,maf,hg19,out_name=out_name,min_diff=min_diff,cov_fraction=cov_fraction,min_reads_cov_fraction=min_reads_cov_fraction)
RF=Realignment_Filter(*sys.argv[1:])
