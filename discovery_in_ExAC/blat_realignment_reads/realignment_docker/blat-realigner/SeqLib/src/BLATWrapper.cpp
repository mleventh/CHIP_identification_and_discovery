/****************************************************************
 ******************* LICENSE AND COPYRIGHT **********************
A large portion of this file and BLATWrapper.cpp are copied
from the BLAT source code, which is copyright of Jim Kent. See 
the BLAT license below:

CONTENTS AND COPYRIGHT

This archive contains the entire source tree for BLAT and
associated utilities.  All files are copyrighted, but license 
is hereby granted for personal, academic, and non-profit use.  
A license is also granted for the contents of the top level 
lib and inc directories for commercial users.  Commercial 
users should contact jim_kent@pacbell.net for access to other modules.
**/

#include "SeqLib/BLATWrapper.h"

#include <set>

#define QWARN_SIZE 5000000

static int ssAliCount = 16;	/* Number of alignments returned by ssStitch. */

namespace SeqLib {

  void BLATWrapper::SetHeaderInfo(const BamHeader& h) {
    
    const bam_hdr_t * t = h.get_();
    for (int i = 0; i < t->n_targets; ++i) 
      m_name2id.insert(std::pair<std::string, int>(std::string(t->target_name[i]),i));

  }
  
  int BLATWrapper::__scoreAli(struct ffAli *ali, boolean isProt, 
			      enum ffStringency stringency, 
			      struct dnaSeq *tSeq, struct trans3 *t3List)
  /* Score alignment. */
  {
    int (*scoreFunc)(char *a, char *b, int size);
    struct ffAli *ff, *nextFf;
    int score = 0;
    if (isProt) 
      scoreFunc = aaScoreMatch;
    else
      scoreFunc = dnaScoreMatch;
    for (ff = ali; ff != NULL; ff = nextFf)
      {
	nextFf = ff->right;
	score += scoreFunc(ff->nStart, ff->hStart, ff->nEnd-ff->nStart);
	if (nextFf != NULL)
	  {
	    int nhStart = trans3GenoPos(nextFf->hStart, tSeq, t3List, FALSE);
	    int ohEnd = trans3GenoPos(ff->hEnd, tSeq, t3List, TRUE);
	    int hGap = nhStart - ohEnd;
	    int nGap = nextFf->nStart - ff->nEnd;
	    score -= ffCalcGapPenalty(hGap, nGap, stringency);
	  }
      }
    return score;
  }
  
  struct ssBundle* BLATWrapper::__fastMapClumpsToBundles(struct genoFind *gf, struct gfClump *clumpList, bioSeq *qSeq)
  /* Convert gfClumps ffAlis. */
  {
    struct gfClump *clump;
    struct gfRange *rangeList = NULL, *range;
    bioSeq *targetSeq;
    struct ssBundle *bunList = NULL, *bun;
    
    for (clump = clumpList; clump != NULL; clump = clump->next)
      __clumpToHspRange(clump, qSeq, gf->tileSize, 0, NULL, &rangeList, FALSE, TRUE);
    slReverse(&rangeList);
    slSort(&rangeList, gfRangeCmpTarget);
    rangeList = gfRangesBundle(rangeList, 256);
    for (range = rangeList; range != NULL; range = range->next)
      {
	targetSeq = range->tSeq;
	bun = (ssBundle*)calloc(1, sizeof(ssBundle)); //(ssBundle*)memset(bun, 0, sizeof(ssBundle)); //needMem(sizeof(*bun)); //AllocVar(bun);
	bun->qSeq = qSeq;
	bun->genoSeq = targetSeq;
	bun->ffList = gfRangesToFfItem(range->components, qSeq);
	bun->isProt = FALSE;
	slAddHead(&bunList, bun);
      }
    gfRangeFreeList(&rangeList);
    return bunList;
  }

  void BLATWrapper::__gfLongDnaInMem(struct dnaSeq *query, struct genoFind *gf, 
				     boolean isRc, int minScore, Bits *qMaskBits, 
				     boolean fastMap, boolean band, BamRecordVector& brv)
  /* Chop up query into pieces, align each, and stitch back
   * together again. */
  {
    
    //printf("query: %s  name: %s  fastMap: %d   band: %d\n", query->dna, query->name, fastMap, band);
    
    int hitCount = 0;
    int maxSize = MAXSINGLEPIECESIZE;
    int preferredSize = 4500;
    int overlapSize = 250;
    struct dnaSeq subQuery = *query;
    struct lm *lm = lmInit(0);
    int subOffset, subSize, nextOffset;
    DNA saveEnd, *endPos;
    struct ssBundle *oneBunList = NULL, *bigBunList = NULL, *bun;
    struct hash *bunHash = newHash(8);

    for (subOffset = 0; subOffset<query->size; subOffset = nextOffset)
      {
	struct gfClump *clumpList;
	struct gfRange *rangeList = NULL;
	
	/* Figure out size of this piece.  If query is
	 * maxSize or less do it all.   Otherwise just
	 * do prefered size, and set it up to overlap
	 * with surrounding pieces by overlapSize.  */
	if (subOffset == 0 && query->size <= maxSize) {
	  nextOffset = subSize = query->size;
	} else
	  {
	    subSize = preferredSize;
	    if (subSize + subOffset >= query->size)
	      {
		subSize = query->size - subOffset;
		nextOffset = query->size;
	      }
	    else
	      {
		nextOffset = subOffset + preferredSize - overlapSize;
	      }
	  }
	subQuery.dna = query->dna + subOffset;
	subQuery.size = subSize;
	endPos = &subQuery.dna[subSize];
	saveEnd = *endPos;
	*endPos = 0;
	//gf->maxPat = 1024; // jeremiah

	if (band)
	  {
	    oneBunList = ffSeedExtInMem(gf, &subQuery, qMaskBits, subOffset, lm, minScore, isRc);
	  }
	else
	  {

	    // body of gfFindClumpWithQmask
	    /*
	    struct gfClump *clumpList = NULL;
	    struct gfHit *hitList;
	    int minMatch = gf->minMatch;
	    
	    hitList =  gfFindHitsWithQmask(gf, &subQuery, qMaskBits, subOffset, lm,
	                &hitCount, NULL, 0, 0);
	    cmpQuerySize = subQuery.size; // seq->size
	    clumpList = clumpHits(gf, hitList, minMatch);
	    */
	    clumpList = gfFindClumpsWithQmask(gf, &subQuery, qMaskBits, subOffset, lm, &hitCount);

	    //printf("HIT COUNT: %d   qMaskMits==null: %d   fastMap: %d, subQuery.dna: %s    offset: %d   gf->tileSize : %d  gf->minMatch+1   %d\n", hitCount, qMaskBits == NULL, fastMap, 
	    //       subQuery.dna, subOffset, gf->tileSize, gf->minMatch+1);
	    //printf("maxPat: %d  minMatch: %d   maxGap: %d  tileSize :%d  stepsize %d tileSpaceSize: %d  tileMask %d  sourceCount: %d isPep %d  allowOneMismatch %d   segSize %d\n",
	    //	   gf->maxPat, gf->minMatch, gf->maxGap, gf->tileSize, gf->stepSize, gf->tileSpaceSize, gf->tileMask, gf->sourceCount, gf->isPep, gf->allowOneMismatch, gf->segSize);
	    if (fastMap)
	      {
		oneBunList = __fastMapClumpsToBundles(gf, clumpList, &subQuery);
	      }
	    else
	      {
		oneBunList = __gfClumpsToBundles(clumpList, isRc, &subQuery, minScore, &rangeList);
		gfRangeFreeList(&rangeList);
	      }
	    gfClumpFreeList(&clumpList);
	  }
	__addToBigBundleList(&oneBunList, bunHash, &bigBunList, query);
	*endPos = saveEnd;
      }
#ifdef DEBUG
    dumpBunList(bigBunList);
#endif /* DEBUG */
    
    for (bun = bigBunList; bun != NULL; bun = bun->next)
      {
	ssStitch(bun, ffCdna, minScore, ssAliCount);
	//if (!fastMap && !band)
	//  __refineSmallExonsInBundle(bun);
	
	//jeremiah
	//printf("bigBunList. Name: %s size: %d\n", bun->genoSeq->name, bun->genoSeq->size); //jeremiah
	
	//jeremiah score it
	struct ssFfItem *ffi;
	for (ffi = bun->ffList; ffi != NULL; ffi = ffi->next) 
	  {

	    struct ffAli *ff_loop = ffi->ff;
	    struct trans3 *t3List = NULL;
	    int score;
	    struct dnaSeq *tSeq = bun->genoSeq; 
	    //struct dnaSeq *qSeq = bun->qSeq;
	    score = __scoreAli(ff_loop, bun->isProt, ffCdna, tSeq, NULL);
	    
	    boolean qIsRc = FALSE;
	    //boolean tIsRc = FALSE;
	    
	    //jeremiah
	    //////////////////
	    //////////////////
	    struct ffAli *ali = ff_loop; //jeremiah
	    struct ffAli *ff, *nextFf;
	    struct ffAli *right = ffRightmost(ali);
	    //enum ffStringency stringency = ffCdna; //jeremiah
	    //int minMatch = minScore; // jeremiah
	    DNA *needle = bun->qSeq->dna;
	    DNA *hay = tSeq->dna;
	    int nStart = ali->nStart - needle;
	    int nEnd = right->nEnd - needle;
	    int nInsertBaseCount = 0;
	    int nInsertCount = 0;
	    int hInsertBaseCount = 0;
	    int hInsertCount = 0;
	    int matchCount = 0;
	    int mismatchCount = 0;
	    int repMatch = 0;
	    int countNs = 0;
	    DNA *np, *hp, n, h;
	    int blockSize;
	    int i;
	    Bits *maskBits = NULL;
	    int chromOffset = 0; //jeremiah
	    
	    int hStart = trans3GenoPos(ali->hStart, tSeq, t3List, FALSE) + chromOffset;
	    //	    int hEnd = trans3GenoPos(right->hEnd, tSeq, t3List, TRUE) + chromOffset;	    

	    /* Count up matches, mismatches, inserts, etc. */
	    for (ff = ali; ff != NULL; ff = nextFf)
	      {
		nextFf = ff->right;
		blockSize = ff->nEnd - ff->nStart;
		np = ff->nStart;
		hp = ff->hStart;
		for (i=0; i<blockSize; ++i)
		  {
		    n = np[i];
		    h = hp[i];
		    if (n == 'n' || h == 'n')
		      ++countNs;
		    else
		      {
			if (n == h)
			  {
			    if (maskBits != NULL)
			      {
				int seqOff = hp + i - hay;
				if (bitReadOne(maskBits, seqOff))
				  ++repMatch;
				else
				  ++matchCount;
			      }
			    else
			      ++matchCount;
			  }
			else
			  ++mismatchCount;
		      }
		  }
		if (nextFf != NULL)
		  {
		    int nhStart = trans3GenoPos(nextFf->hStart, tSeq, t3List, FALSE) + chromOffset;
		    int ohEnd = trans3GenoPos(ff->hEnd, tSeq, t3List, TRUE) + chromOffset;
		    int hGap = nhStart - ohEnd;
		    int nGap = nextFf->nStart - ff->nEnd;
		    
		    if (nGap != 0)
		      {
			++nInsertCount;
			nInsertBaseCount += nGap;
		      }
		    if (hGap != 0)
		      {
			++hInsertCount;
			hInsertBaseCount += hGap;
		      }
		  }
	      }

	    // move on if score too low
	    if (score < 30)
	      continue;
	    
	    Cigar cigv;

	    // make the cigar (idea from Heng Li's psl2sam.pl)
	    std::vector<long> x, y, z;
	    for (ff = ali; ff != NULL; ff = ff->right)
	      x.push_back((long)(ff->nEnd - ff->nStart));
	    for (ff = ali; ff != NULL; ff = ff->right)
	      y.push_back((long)(ff->nStart - needle));
	    for (ff = ali; ff != NULL; ff = ff->right)
	      z.push_back(trans3GenoPos(ff->hStart, tSeq, t3List, FALSE) + chromOffset);

	    // 5' clipping
	    if (nStart > 0)
	      cigv.add(CigarField('S', nStart));

	    int gap_open = 0;
	    int gap_ext = 0;
	    int y0 = y[0], z0 = z[0]; 
	    
	    int blockCount = ffAliCount(ali);

	    for (int i = 1; i < blockCount; ++i)
	      {
		int ly = y[i] - y[i-1] - x[i-1];
		int lz = z[i] - z[i-1] - x[i-1];
		if (ly < lz) // del: the reference gap is longer
		  {
		    ++gap_open;
		    gap_ext += lz - ly;
		    cigv.add(CigarField('M', y[i] - y0));
		    cigv.add(CigarField('D', lz - ly));
		    //cig += std::to_string(y[i] - y0) + "M";
		    //cig += std::to_string(lz - ly) + "D";
		    y0 = y[i];
		    z0 = z[i];
		  }
		else if (lz < ly)  // ins: query gap is longer
		  {
		    ++gap_open;
		    gap_ext += ly - lz;
		    cigv.add(CigarField('M', z[i] - z0));
		    cigv.add(CigarField('I', ly - lz));
		    //cig += std::to_string(z[i] - z0) + "M";
		    //cig += std::to_string(ly - lz) + "I";
		    y0 = y[i];
		    z0 = z[i];
		  }
	      }
	    cigv.add(CigarField('M', nEnd - y0));
	    //cig += std::to_string(nEnd - y0) + "M";
	    if (query->size != nEnd) { // 3'-end clipping
	      cigv.add(CigarField('S', query->size - nEnd));
	      //cig += std::to_string(query->size - nEnd) + "S";
	    }

	    SeqHashMap<std::string, int>::const_iterator fff = m_name2id.find(std::string(bun->genoSeq->name));
	    if (fff == m_name2id.end())
	      throw std::out_of_range("No header entry for " + std::string(bun->genoSeq->name));
	    
	    //std::shared_ptr<bam1_t> raw_b;
	    bam1_t * raw_b = bam_init1();

	    raw_b->core.tid = fff->second; 
	    raw_b->core.pos = hStart;
	    raw_b->core.qual = 0; //a.mapq; // TODO
	    raw_b->core.flag = 0; //a.flag;
	    raw_b->core.n_cigar = cigv.size(); //a.n_cigar;
	    
	    // set dumy mate
	    raw_b->core.mtid = -1;
	    raw_b->core.mpos = -1;
	    raw_b->core.isize = 0;

	    // if alignment is reverse, set it
	    if (qIsRc) 
	      raw_b->core.flag |= BAM_FREVERSE;

	    // allocate all the data
	    raw_b->core.l_qname = strlen(query->name) + 1; //name.length() + 1;
	    raw_b->core.l_qseq = query->size; //seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
	    raw_b->l_data = raw_b->core.l_qname + (raw_b->core.n_cigar<<2) + ((raw_b->core.l_qseq+1)>>1) + (raw_b->core.l_qseq);
	    raw_b->data = (uint8_t*)malloc(raw_b->l_data);

	    // allocate the qname
	    memcpy(raw_b->data, query->name, raw_b->core.l_qname);

	    // allocate the sequence
	    uint8_t* m_bases = raw_b->data + raw_b->core.l_qname + (raw_b->core.n_cigar<<2);
	    
	    // TODO move this out of bigger loop
	    int slen = query->size; //seq.length();
	    int j = 0;
	    if (qIsRc/* && false*/) {
	      for (int i = slen-1; i >= 0; --i) {
		
		// bad idea but works for now
		// this is REV COMP things
		uint8_t base = 15;
		if (query->dna[i] == 't') //seq.at(i) == 'T')
		  base = 1;
		else if (query->dna[i] == 'g') //seq.at(i) == 'G')
		  base = 2;
		else if (query->dna[i] == 'c') //seq.at(i) == 'C')
		  base = 4;
		else if (query->dna[i] == 'a') //seq.at(i) == 'A')
		  base = 8;
		
		m_bases[j >> 1] &= ~(0xF << ((~j & 1) << 2));   ///< zero out previous 4-bit base encoding
		m_bases[j >> 1] |= base << ((~j & 1) << 2);  ///< insert new 4-bit base encoding
		++j;
	      }
	    } else {
	      for (int i = 0; i < slen; ++i) {
		// bad idea but works for now
		uint8_t base = 15;
		if (query->dna[i] == 'a') //seq.at(i) == 'a')
		  base = 1;
		else if (query->dna[i] == 'c') //seq.at(i) == 'C')
		  base = 2;
		else if (query->dna[i] == 'g') //seq.at(i) == 'G')
		  base = 4;
		else if (query->dna[i] == 't') //seq.at(i) == 'T')
		  base = 8;
		
		m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
		m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
		
	      }
	    }
	    
	    // allocate the quality to NULL
	    uint8_t* s = bam_get_qual(raw_b);
	    s[0] = 0xff;
	    
	    // allocate the cigar. 32 bits per elem (4 type, 28 length)

	    uint32_t * cigr = bam_get_cigar(raw_b);
	    for (size_t i = 0; i < cigv.size(); ++i)
	      cigr[i] = cigv[i].raw(); //Length << BAM_CIGAR_SHIFT | BAM_CMATCH;

	    // assign the read by making a copy of this bam1_t
	    BamRecord b;
	    b.assign(raw_b);
	    //free(raw_b);

	    // add the alignment score, as calculate by psl2sam
	    int a_param = 1;
	    int b_param = 3;
	    int q_param = 5;
	    int r_param = 2;

	    // add the score
	    int as_score = a_param * score - b_param * mismatchCount - q_param * gap_open - r_param * gap_ext;
	    b.AddIntTag("AS", as_score);
	    b.AddIntTag("BS", score);
	    b.AddIntTag("SR", brv.size());

	    brv.push_back(b);

	    //printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%s\t%d\t%d\n", score, mismatchCount, repMatch, countNs,
	    //nInsertCount, nInsertBaseCount, hInsertCount, hInsertBaseCount, (qIsRc ? '-' : '+'),
	    //   qSeq->name, qSeq->size, nStart, nEnd, bun->genoSeq->name, hStart, hEnd);
	  }
	///////////////
	///////////////
	
	//saveAlignments(bun->genoSeq->name, bun->genoSeq->size, 0, bun, NULL, isRc, FALSE, ffCdna, minScore, out);
      }

    ssBundleFreeList(&bigBunList);
    freeHash(&bunHash);
    lmCleanup(&lm);
  }
  
  
  void BLATWrapper::searchOneStrand(struct dnaSeq *seq, struct genoFind *gf, 
				    boolean isRc, Bits *qMaskBits, BamRecordVector& brv)
  /* Search for seq in index, align it, and write results to psl. */
  {
    if (fastMap && (seq->size > MAXSINGLEPIECESIZE))
      errAbort("Maximum single piece size (%d) exceeded by query %s of size (%d). "
	       "Larger pieces will have to be split up until no larger than this limit "
	       "when the -fastMap option is used."	
	       , MAXSINGLEPIECESIZE, seq->name, seq->size);
    __gfLongDnaInMem(seq, gf, isRc, minScore, qMaskBits, fastMap, false /* fine? */, brv);
  }
 
  Bits* BLATWrapper::maskQuerySeq(struct dnaSeq *seq, boolean isProt, 
		     boolean maskQuery, boolean lcMask)
  /* Massage query sequence a bit, converting it to correct
   * case (upper for protein/lower for DNA) and optionally
   * returning upper/lower case info , and trimming poly A. */
  {
    Bits *qMaskBits = NULL;
    verbose(2, "%s\n", seq->name);
    if (isProt)
      faToProtein(seq->dna, seq->size);
    else
      {
	if (maskQuery)
	  {
	    if (lcMask)
	      toggleCase(seq->dna, seq->size);
	    qMaskBits = maskFromUpperCaseSeq(seq);
	  }
	faToDna(seq->dna, seq->size);
      }
    if (seq->size > QWARN_SIZE)
      {
	warn("Query sequence %s has size %d, it might take a while.",
	     seq->name, seq->size);
      }
    return qMaskBits;
  }
  
  
  void BLATWrapper::searchOne(bioSeq *seq, struct genoFind *gf, struct hash *maskHash, Bits *qMaskBits, BamRecordVector& brv)
  /* Search for seq on either strand in index. */
  {
    //if (isProt)
    //  {
    //	searchOneProt(seq, gf, f);
    // }
    //else
      {
	gvo->maskHash = maskHash; //debug
	searchOneStrand(seq, gf, FALSE, qMaskBits, brv);
	reverseComplement(seq->dna, seq->size);
	searchOneStrand(seq, gf, TRUE, qMaskBits, brv);
	reverseComplement(seq->dna, seq->size);
      }

      FILE * f = NULL;
    gfOutputQuery(gvo, f);
  }
  

  void BLATWrapper::trimSeq(struct dnaSeq *seq, struct dnaSeq *trimmed)
  /* Copy seq to trimmed (shallow copy) and optionally trim
   * off polyA tail or polyT head. */
  {
    DNA *dna = seq->dna;
    int size = seq->size;
    *trimmed = *seq;
    if (trimT)
      maskHeadPolyT(dna, size);
    if (trimA || trimHardA)
      {
	int trimSize = maskTailPolyA(dna, size);
	if (trimHardA)
	  {
	    trimmed->size -= trimSize;
	    dna[size-trimSize] = 0;
	  }
      }
  }
  

  void BLATWrapper::searchOneMaskTrim(struct dnaSeq *seq, boolean isProt,
				      struct genoFind *gf,
		       struct hash *maskHash,
				      long long *retTotalSize, int *retCount, BamRecordVector& brv)

  /* Search a single sequence against a single genoFind index. */
  {
    boolean maskQuery = false ; //(qMask != NULL);
    boolean lcMask = false ; //(qMask != NULL && sameWord(qMask, "lower"));
    Bits *qMaskBits = maskQuerySeq(seq, isProt, maskQuery, lcMask);
    struct dnaSeq trimmedSeq;
    ZeroVar(&trimmedSeq);
    trimSeq(seq, &trimmedSeq);
    if (qType == gftRna || qType == gftRnaX)
      memSwapChar(trimmedSeq.dna, trimmedSeq.size, 'u', 't');

    searchOne(&trimmedSeq, gf,maskHash, qMaskBits, brv);
    *retTotalSize += seq->size;
    *retCount += 1;
    bitFree(&qMaskBits);
  }
  
  
  void BLATWrapper::__searchOneIndex(int fileCount, char *files[], struct genoFind *gf, char *outName, 
				     struct hash *maskHash, FILE *outFile, boolean showStatus)
  /* Search all sequences in all files against single genoFind index. */
  {
    int i;
    char *fileName;
    int count = 0; 
    long long totalSize = 0;
    
    std::cerr << "file count " << fileCount << std::endl;

    //gfOutputHead(gvo, outFile);
    
    for (i=0; i<fileCount; ++i)
      {
	fileName = files[i];
	if (nibIsFile(fileName))
	  {
	    /*	    struct dnaSeq *seq;
	    
	    if (isProt)
	      errAbort("%s: Can't use .nib files with -prot or d=prot option\n", fileName);
	    seq = nibLoadAllMasked(NIB_MASK_MIXED, fileName);
	    freez(&seq->name);
	    seq->name = cloneString(fileName);
	    searchOneMaskTrim(seq, isProt, gf, outFile,
			      maskHash, &totalSize, &count);
	    freeDnaSeq(&seq);
	    */
	  }
	else if (twoBitIsSpec(fileName))
	  {
	    /*
	      struct twoBitSpec *tbs = twoBitSpecNew(fileName);
	      struct twoBitFile *tbf = twoBitOpen(tbs->fileName);
	      if (isProt)
	      errAbort("%s is a two bit file, which doesn't work for proteins.", 
		       fileName);
	    if (tbs->seqs != NULL)
	      {
		struct twoBitSeqSpec *ss = NULL;
		for (ss = tbs->seqs;  ss != NULL;  ss = ss->next)
		  {
		    struct dnaSeq *seq = twoBitReadSeqFrag(tbf, ss->name,
							   ss->start, ss->end);
		    searchOneMaskTrim(seq, isProt, gf, outFile,
				      maskHash, &totalSize, &count);
		    dnaSeqFree(&seq);
		  }
	      }
	    else
	      {
		  struct twoBitIndex *index = NULL;
		  for (index = tbf->indexList; index != NULL; index = index->next)
		  {
		  struct dnaSeq *seq = twoBitReadSeqFrag(tbf, index->name, 0, 0);
		  searchOneMaskTrim(seq, isProt, gf, outFile,
		  maskHash, &totalSize, &count);
		  dnaSeqFree(&seq);
		  }
	      }
	    twoBitClose(&tbf);
	    */
	  }
	else
	  {
	    static struct dnaSeq seq;
	    struct lineFile *lf = lineFileOpen(fileName, TRUE);
	    while (faMixedSpeedReadNext(lf, &seq.dna, &seq.size, &seq.name))
	      {
		//searchOneMaskTrim(&seq, isProt, gf, outFile,
		//		  maskHash, &totalSize, &count);
	      }
	    lineFileClose(&lf);
	  }
      }
    carefulClose(&outFile);
    if (showStatus)
      printf("Searched %lld bases in %d sequences\n", totalSize, count);
  }
  
  
  void BLATWrapper::LoadIndex(const std::string& file, const std::string& oocfile) {
    
    int dbCount;
    bool tIsProt = false;
    
    // copy over from const char to char
    char *dbFile = (char*)malloc(file.length() + 1);
    strcpy(dbFile, file.c_str());
    ooc = (char*)malloc(oocfile.length() + 1);
    strcpy(ooc, oocfile.c_str());
    
    bool showStatus = true;
    
    gfClientFileArray(dbFile, &dbFiles, &dbCount);

    dbSeqList = gfClientSeqList(dbCount, dbFiles, tIsProt, false, repeats,
				minRepDivergence, showStatus);

    int databaseSeqCount = slCount(dbSeqList);
    bool qIsProt = false;
    char *databaseName = dbFile;

    unsigned long databaseLetters = 0;	/* Number of bases in database. */

    struct dnaSeq *seq;
    for (seq = dbSeqList; seq != NULL; seq = seq->next)
      databaseLetters += seq->size;

    FILE * f = NULL;

    gvo = gfOutputAny(outputFormat, minIdentity*10, qIsProt, tIsProt, noHead, 
		      databaseName, databaseSeqCount, databaseLetters, minIdentity, f);

    gf = gfIndexSeq(dbSeqList, minMatch, maxGap, tileSize, repMatch, ooc, 
		    tIsProt, oneOff, false, stepSize);

  }
  
  void BLATWrapper::AlignSequence(const std::string& name, const std::string& sequence, BamRecordVector& brv) {
    
    long long totalSize = 0;
    int count = 0;
    struct hash *maskHash = NULL;

    // convert to lowercase
    std::string lower = sequence;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

        // make the dna seq
    struct dnaSeq seq;
    seq.name = (char*)malloc(name.length() + 1);
    memcpy(seq.name, name.c_str(), name.length() + 1);
    seq.dna = (char*)malloc(sequence.length() + 1);
    memcpy(seq.dna, lower.c_str(), sequence.length());  //cloneString();
    seq.dna[sequence.length()] = 0;
    seq.size = sequence.length();
    seq.next = NULL;
    seq.mask = NULL;

    //setup the output

    // search the sequence
    bool isProt = false;
    FILE * outFile = NULL;
    BamRecordVector bavr;
    searchOneMaskTrim(&seq, isProt, gf,maskHash, &totalSize, &count, bavr);

    free(seq.name);
    free(seq.dna);

    if (!bavr.size())
      return;

    std::vector<BamRecord> best_hit;  // which contig is best at loc
    std::vector<int> best_score;      // score of that contig

    // get covered positions
    /*    for (int i = 0; i < sequence.length(); ++i) {
      best_hit.push_back(bavr[0]);
      int bs = bavr[0].GetIntTag("BS");
      for (auto& r : bavr) {
	if (r.coveredBase(i) && r.GetIntTag("BS") > bs) {
	  best_hit[i] = r;
	  bs = r.GetIntTag("BS");
	}
      }
    }

    std::set<int> sset;
    for (auto& i : best_hit) {		       
      int sr = i.GetIntTag("SR");
      if (!sset.count(sr))
	{
	  brv.push_back(i);
	  sset.insert(sr);
	}
	}*/
    brv = bavr;
  }

  void BLATWrapper::QueryFile(const std::string& file, const std::string& ofile) {
    
    int queryCount;
    bool showStatus = true;
    struct hash *maskHash = NULL;
    
    char** queryFiles;
    
    // copy over from const char to char
    char *queryFile = (char*)malloc(file.length() + 1);
    strcpy(queryFile, file.c_str());

    char *outName = (char*)malloc(ofile.length() + 1);
    strcpy(outName, ofile.c_str());

    FILE *f = mustOpen(outName, "w");
    
    gfClientFileArray(queryFile, &queryFiles, &queryCount);

    __searchOneIndex(queryCount, queryFiles, gf, outName, maskHash, f, showStatus);

    free(queryFile);
    free(outName);

  }

  void BLATWrapper::__addToBigBundleList(struct ssBundle **pOneList, struct hash *bunHash, 
					 struct ssBundle **pBigList, struct dnaSeq *query)
  /* Add bundles in one list to bigList, consolidating bundles that refer
   * to the same target sequence.  This will destroy oneList in the process. */
  {
    struct ssBundle *oneBun, *bigBun;
    
    for (oneBun = *pOneList; oneBun != NULL; oneBun = oneBun->next)
      {
	char *name = oneBun->genoSeq->name;
	if ((bigBun = (ssBundle*)hashFindVal(bunHash, name)) == NULL)
	  {
	    //AllocVar(bigBun);
	    bigBun = (ssBundle*)calloc(1, sizeof(ssBundle)); //(ssBundle*)memset(bigBun, 0, sizeof(ssBundle));//jeremiah
	    slAddHead(pBigList, bigBun);
	    hashAdd(bunHash, name, bigBun);
	    bigBun->qSeq = query;
	    bigBun->genoSeq = oneBun->genoSeq;
	    bigBun->isProt = oneBun->isProt;
	    bigBun->avoidFuzzyFindKludge = oneBun->avoidFuzzyFindKludge;
	  }
	bigBun->ffList = (ssFfItem*)slCat(bigBun->ffList, oneBun->ffList); //jeremiah
	oneBun->ffList = NULL;
      }
    ssBundleFreeList(pOneList);
  }

  struct ssBundle* BLATWrapper::__gfClumpsToBundles(struct gfClump *clumpList, 
						    boolean isRc, struct dnaSeq *seq, int minScore,  
						    struct gfRange **retRangeList)
  /* Convert gfClumps to an actual alignments (ssBundles) */ 
  {
    struct ssBundle *bun, *bunList = NULL;
    struct gfRange *rangeList = NULL, *range;
    struct dnaSeq *targetSeq;
    
    int usualExpansion = 500; //jeremiah

    rangeList = __seqClumpToRangeList(clumpList, 0);
    slSort(&rangeList, gfRangeCmpTarget);
    rangeList = gfRangesBundle(rangeList, 2000);
    for (range = rangeList; range != NULL; range = range->next)
      {
    targetSeq = range->tSeq;
    gfiExpandRange(range, seq->size, targetSeq->size, FALSE, isRc, 
		   usualExpansion);
    range->tStart = 0;
    range->tEnd = targetSeq->size;
    bun = (ssBundle*)calloc(1, sizeof(ssBundle)); //memset(bun, 0, sizeof(ssBundle)); //AllocVar(bun);
    bun->qSeq = seq;
    bun->genoSeq = targetSeq;
    __alignComponents(range, bun, ffCdna);
    ssStitch(bun, ffCdna, minScore, ssAliCount);
    slAddHead(&bunList, bun);
      }
    slReverse(&bunList);
    *retRangeList = rangeList;
    return bunList;
  }
  
  
  struct gfRange* BLATWrapper::__seqClumpToRangeList(struct gfClump *clumpList, int frame)
  /* Convert from clump list to range list. */
  {
    struct gfRange *rangeList = NULL, *range;
    struct gfClump *clump;
    char *name;
    int tOff;
    
    for (clump = clumpList; clump != NULL; clump = clump->next)
      {
	tOff = clump->target->start;
	//AllocVar(range);
	range = (gfRange*)calloc(1,sizeof(gfRange)); //(gfRange*)memset(range, 0, sizeof(gfRange));
	range->qStart = clump->qStart;
	range->qEnd = clump->qEnd;
	name = __clumpTargetName(clump);
	range->tName = cloneString(name);
	range->tStart = clump->tStart - tOff;
	range->tEnd = clump->tEnd - tOff;
	range->tSeq = clump->target->seq;
	range->frame = frame;
	slAddHead(&rangeList, range);
      }
    slReverse(&rangeList);
    return rangeList;
  }
  
  boolean BLATWrapper::__alignComponents(struct gfRange *combined, struct ssBundle *bun, 
					 enum ffStringency stringency)
  /* Align each piece of combined->components and put result in
   * bun->ffList. */
  {
    struct gfRange *range;
    struct dnaSeq *qSeq = bun->qSeq, *tSeq = bun->genoSeq;
    struct ssFfItem *ffi;
    struct ffAli *ali;
    int qStart, qEnd, tStart, tEnd;
    int extra = 250;
    boolean gotAny = FALSE;
    
    for (range = combined->components; range != NULL; range = range->next)
      {
	/* Expand to include some extra sequence around range. */
	qStart = range->qStart - extra;
	tStart = range->tStart - extra;
	qEnd = range->qEnd + extra;
	tEnd = range->tEnd + extra;
	if (range == combined->components)
	  {
	    qStart -= extra;
	    tStart -= extra;
	  }
	if (range->next == NULL)
	  {
	    qEnd += extra;
	    tEnd += extra;
	  }
	if (qStart < combined->qStart) qStart = combined->qStart;
	if (tStart < combined->tStart) tStart = combined->tStart;
	if (qEnd > combined->qEnd) qEnd = combined->qEnd;
	if (tEnd > combined->tEnd) tEnd = combined->tEnd;
	ali = ffFind(qSeq->dna + qStart,
		     qSeq->dna + qEnd,
		     tSeq->dna + tStart - combined->tStart,
		     tSeq->dna + tEnd - combined->tStart,
		     stringency);
	if (ali != NULL)
	  {
	    //AllocVar(ffi);
	    ffi = (ssFfItem*)calloc(1, sizeof(ssFfItem));//(ssFfItem*)memset(ffi, 0, sizeof(ssFfItem));
	    ffi->ff = ali;
	    slAddHead(&bun->ffList, ffi);
	    gotAny = TRUE;
	  }
      }
    return gotAny;
  }
  
  char* BLATWrapper::__clumpTargetName(struct gfClump *clump)
  /* Return target name of clump - whether it is in memory or on disk. */
  {
    struct gfSeqSource *target = clump->target;
    char *name;
    if (target->seq != NULL)
      name = target->seq->name;
    else
      name = target->fileName;
    if (name == NULL)
      internalErr();
    return name;
  }

void BLATWrapper::__clumpToHspRange(struct gfClump *clump, bioSeq *qSeq, int tileSize,
				    int frame, struct trans3 *t3, struct gfRange **pRangeList, 
				    boolean isProt, boolean fastMap)
/* Covert clump->hitList to HSPs (high scoring local sequence pair,
 * that is longest alignment without gaps) and add resulting HSPs to
 * rangeList. */
{
  struct gfSeqSource *target = clump->target;
  aaSeq *tSeq = target->seq;
  BIOPOL *qs, *ts, *qe, *te;
  struct gfHit *hit;
  int qStart = 0, tStart = 0, qEnd = 0, tEnd = 0, newQ = 0, newT = 0;
  boolean outOfIt = TRUE;		/* Logically outside of a clump. */
  struct gfRange *range;
  BIOPOL *lastQs = NULL, *lastQe = NULL, *lastTs = NULL; //, *lastTe = NULL;
  int (*scoreMatch)(char a, char b) = (isProt ? aaScore2 : dnaScore2);
  int maxDown, minSpan;
  
  if (fastMap)
    {
      maxDown = 1;
      minSpan = 50;
    }
  else
    {
      maxDown = 10;
      minSpan = 0;
    }
  
  
  if (tSeq == NULL)
    internalErr();
  
  /* The termination condition of this loop is a little complicated.
   * We want to output something either when the next hit can't be
   * merged into the previous, or at the end of the list.  To avoid
   * duplicating the output code we're forced to complicate the loop
   * termination logic.  Hence the check for hit == NULL to break
   * the loop is not until near the end of the loop. */
  for (hit = clump->hitList; ; hit = hit->next)
    {
      if (hit != NULL)
        {
	  newQ = hit->qStart;
	  newT = hit->tStart - target->start;
	}
      
      /* See if it's time to output merged (diagonally adjacent) hits. */
      if (!outOfIt)	/* Not first time through. */
        {
	/* As a micro-optimization handle strings of adjacent hits
	 * specially.  Don't do the extensions until we've merged
	 * all adjacent hits. */
	  if (hit == NULL || newQ != qEnd || newT != tEnd)
	    {
	      qs = qSeq->dna + qStart;
	      ts = tSeq->dna + tStart;
	      qe = qSeq->dna + qEnd;
	      te = tSeq->dna + tEnd;
	      __extendHitRight(qSeq->size - qEnd, tSeq->size - tEnd,
			     &qe, &te, scoreMatch, maxDown);
	      __extendHitLeft(qStart, tStart, &qs, &ts, scoreMatch, maxDown);
	      if (qs != lastQs || ts != lastTs || qe != lastQe || qs !=  lastQs)
		{
		  lastQs = qs;
		  lastTs = ts;
		  lastQe = qe;
		  //lastTe = te;
		  if (qe - qs >= minSpan)
		    {
		      //AllocVar(range);
		      range = (gfRange*)calloc(1, sizeof(gfRange)); //(gfRange*)memset(range, 0, sizeof(gfRange));
		      range->qStart = qs - qSeq->dna;
		      range->qEnd = qe - qSeq->dna;
		      range->tName = cloneString(tSeq->name);
		      range->tSeq = tSeq;
		      range->tStart = ts - tSeq->dna;
		      range->tEnd = te - tSeq->dna;
		      range->hitCount = qe - qs;
		      range->frame = frame;
		      range->t3 = t3;
		      assert(range->tEnd <= tSeq->size);
		      slAddHead(pRangeList, range);
		    }
		}
	      outOfIt = TRUE;
	    }
	}
      if (hit == NULL)
        break;
      
      if (outOfIt)
        {
	  qStart = newQ;
	  qEnd = qStart + tileSize;
	  tStart = newT;
	  tEnd = tStart + tileSize;
	  outOfIt = FALSE;
	}
      else
        {
	  qEnd = newQ + tileSize;
	  tEnd = newT + tileSize;
	}
    }
}

void BLATWrapper::__extendHitRight(int qMax, int tMax,
				   char **pEndQ, char **pEndT, int (*scoreMatch)(char a, char b), 
				   int maxDown)
/* Extend endQ/endT as much to the right as possible. */
{
  int maxScore = 0;
  int score = 0;
  int maxPos = -1;
  int last = Blatmin(qMax, tMax);
  int i;
  char *q = *pEndQ, *t = *pEndT;
  
  for (i=0; i<last; ++i)
    {
      score += scoreMatch(q[i], t[i]);
      if (score > maxScore)
	{
	  maxScore = score;
	  maxPos = i;
	 }
      else if (i > maxPos + maxDown)
	{
	  break;
	}
    }
  *pEndQ = q+maxPos+1;
  *pEndT = t+maxPos+1;
}

/*
struct gfHit* BlatWrapper::__gfFindHitsWithQmask(struct genoFind *gf, bioSeq *seq,
						 Bits *qMaskBits, int qMaskOffset, struct lm *lm, int *retHitCount, 
						 struct gfSeqSource *target, int tMin, int tMax)
// Find hits associated with one sequence soft-masking seq according to qMaskBits.
// The hits will be in genome rather than chromosome coordinates. 
{
  
  struct gfHit *hitList = NULL;
  if (gf->segSize == 0 && !gf->isPep && !gf->allowOneMismatch)
    {
      hitList = gfFastFindDnaHits(gf, seq, qMaskBits, qMaskOffset, lm, retHitCount,
				  target, tMin, tMax);
    }
  else
    {
    if (gf->segSize == 0)
      {
	if (gf->allowOneMismatch)
	  {
	    hitList = gfStraightFindNearHits(gf, seq, qMaskBits, qMaskOffset, lm, 
					     retHitCount, target, tMin, tMax);
	  }
	else
	  {
	    hitList = gfStraightFindHits(gf, seq, qMaskBits, qMaskOffset, lm, retHitCount, 
					 target, tMin, tMax);
	  }
      }
    else
      {
	if (gf->allowOneMismatch)
	  {
	    hitList = gfSegmentedFindNearHits(gf, seq, qMaskBits, qMaskOffset, lm, retHitCount,
					      target, tMin, tMax);
	  }
	else
	  {
	    hitList = gfSegmentedFindHits(gf, seq, qMaskBits, qMaskOffset, lm, retHitCount,
					  target, tMin, tMax);
	  }
      }
    }
  return hitList;
}
*/

void BLATWrapper::__extendHitLeft(int qMax, int tMax,
				  char **pStartQ, char **pStartT, int (*scoreMatch)(char a, char b),
				  int maxDown)
// Extend startQ/startT as much to the left as possible. 
{
  int maxScore = 0;
  int score = 0;
  int maxPos = 0;
  int last = -Blatmin(qMax, tMax);
  int i;
  char *q = *pStartQ, *t = *pStartT;
  
  for (i=-1; i>=last; --i)
    {
      score += scoreMatch(q[i], t[i]);
      if (score > maxScore)
	{
	  maxScore = score;
	  maxPos = i;
	}
      else if (i < maxPos - maxDown)
	{
	  break;
	}
    }
  *pStartQ = q+maxPos;
  *pStartT = t+maxPos;
}


/*
struct gfHit* BLATWrapper::__gfFastFindDnaHits(struct genoFind *gf, struct dnaSeq *seq, 
					       Bits *qMaskBits,  int qMaskOffset, struct lm *lm, int *retHitCount,
					       struct gfSeqSource *target, int tMin, int tMax)
// Find hits associated with one sequence. This is is special fast
// case for DNA that is in an unsegmented index. 
{
  struct gfHit *hitList = NULL, *hit;
  int size = seq->size;
  int tileSizeMinusOne = gf->tileSize - 1;
  int mask = gf->tileMask;
  DNA *dna = seq->dna;
  int i, j;
  bits32 bits = 0;
  bits32 bVal;
  int listSize;
  bits32 qStart, *tList;
  int hitCount = 0;
  
  for (i=0; i<tileSizeMinusOne; ++i)
    {
      bVal = ntValNoN[(int)dna[i]];
      bits <<= 2; 
      bits += bVal;
    }
  
  for (i=tileSizeMinusOne; i<size; ++i)
    {
      bVal = ntValNoN[(int)dna[i]];
      bits <<= 2;
      bits += bVal;
      bits &= mask;
      listSize = gf->listSizes[bits];
      
      if (listSize != 0)
	{
	  qStart = i-tileSizeMinusOne;
	  if (qMaskBits == NULL || bitCountRange(qMaskBits, qStart+qMaskOffset, gf->tileSize) == 0)
	    {
	      tList = gf->lists[bits];
	      for (j=0; j<listSize; ++j)
		{
		  int tStart = tList[j];
		  if (target == NULL || 
		      (target == __findSource(gf, tStart) && tStart >= tMin && tStart < tMax) ) 
		    {
		      
		      lmAllocVar(lm, hit);
		      hit->qStart = qStart;
		      hit->tStart = tStart;
		      hit->diagonal = tStart + size - qStart;
		      slAddHead(&hitList, hit);
		      ++hitCount;
		    }
		}
	    }
	}
    }
  *retHitCount = hitCount;
  return hitList;
}

struct gfSeqSource* BLATWrapper::__findSource(struct genoFind *gf, bits32 targetPos)
// Find source given target position. 
{
  //printf("running find source on pos: %d\n");
  struct gfSeqSource *ss =  bsearch(&targetPos, gf->sources, gf->sourceCount, 
				    sizeof(gf->sources[0]), bCmpSeqSource);
  if (ss == NULL)
    errAbort("Couldn't find source for %d", targetPos);
  return ss;
}
*/

}
