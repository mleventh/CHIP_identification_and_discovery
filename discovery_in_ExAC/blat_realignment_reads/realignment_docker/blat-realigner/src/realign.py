import subprocess,pysam,time,tempfile,sys,os,random,string,socket


import blat_realign_cython
class Wrapper: #for running in FH
    @staticmethod
    def run_realign(bam,maf,hg19="/opt/hg19.2bit",out_name=None,min_diff=5,cov_fraction=0.05,min_reads_cov_fraction=5):
        Realignment_Filter(bam,maf,hg19,out_name=out_name,min_diff=min_diff,cov_fraction=cov_fraction,min_reads_cov_fraction=min_reads_cov_fraction)


class BLAT_server:
    def __init__(self,ref_2bit,port=random.randint(58142,65530),debug=False):
        self.server=blat_realign_cython.pyRealignmentFilter(ref_2bit,"/opt/hg19.11.ooc")


    def query(self,sequence,chrN,pos):
        return(self.server.blat_seq(sequence,chrN,pos))


class Realignment_Filter:
    @staticmethod
    def get_stats(score_diff_list):
        return map(str, [min(score_diff_list), sum(score_diff_list) / float(len(score_diff_list)), max(score_diff_list)]) if len(score_diff_list) else ["", "", ""]

    def __init__(self,bam,maf,out_name=None,hg19="/opt/hg19.2bit",min_diff=5,cov_fraction=0.05,min_reads_cov_fraction=5):
        self.input_bam = pysam.AlignmentFile(bam, "rb")
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

        for pileupcolumn in self.input_bam.pileup(chrm, pos-1, pos,max_depth=2000000,multiple_iterators=True):

            if pileupcolumn.pos != pos-1:
                continue


            for read in pileupcolumn.pileups:
                base_score_filter = ((ord(read.alignment.qual[read.query_position])-33)<= 5) if read.query_position is not None else False #read.query_position can be none for deletions
                base_score_filter = True if read.query_position is None and alt is not "-" else base_score_filter #If the read has "none" but the mutation isn't a deletion, skip this read since it's not the alt we want.
                if read.alignment.qname in read_names or read.alignment.is_duplicate or read.alignment.mapq <= 5 or read.alignment.is_qcfail or read.alignment.is_secondary or read.alignment.is_unmapped or base_score_filter:
                    #print "read failed quality metrics and therefore not counting",
                    #print  (read.alignment.qname in read_names , read.alignment.is_duplicate , read.alignment.mapq <= 5 , read.alignment.is_qcfail , read.alignment.is_secondary , read.alignment.is_unmapped , base_score_filter)
                    continue


                is_alt=False
                if alt != "-" and ref != "-": #indels handled below
                    is_alt=True if read.alignment.seq[read.query_position] == alt else False

                elif len(read.alignment.blocks)>1:
                    #print read.cigar,read.pos, read.qlen,read.qend,read.blocks
                    #print read.blocks

                    blk_prev=read.alignment.blocks[0][1]
                    for blck in read.alignment.blocks[1:]:
                            blck_size=blck[0]-blk_prev
                            if blck_size==0: #insertion
                                if blk_prev==pos-1: is_alt=True

                            else: #deletion
                                if blk_prev<=pos-1<=blk_prev+blck_size : is_alt=True

                            blk_prev=blck[1]

                #if alt == "-":
                #    is_alt=True if read.is_del == 1 | read.indel else False
                #elif ref == "-":
                #    is_alt=True if self.position_to_ins_del(read.alignment.pos,pos,read.alignment.cigar) == "1" else False



                #if read.query_position or is_alt:
                read_names.add(read.alignment.qname)

                if is_alt:
                    #read_names.add(read.alignment.qname)
                    blat_res=self.blat_read(read.alignment.seq,[chrm,pos])
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
                        try:
                           mate=self.input_bam.mate(read.alignment) #fetch the mate from the bam (might be slow?)
                           #read.reference_start #get the start position to see if it maps away from it
                           mate_res=self.blat_read(mate.query_alignment_sequence,[chrm,mate.reference_start+2])
                           with_minus_without.append(max([blat_res[1],mate_res[1]]))
                           if mate_res[0] > 0: #if the mate maps stays:
                               #print "read discordant, mate not discordant!"
                               mate_not_discordant+=1
                           else: #if mate maps away from it's start
                               discordant +=1
                        except:
                           with_minus_without.append(blat_res[1])
                           print "failed to get mate"
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
