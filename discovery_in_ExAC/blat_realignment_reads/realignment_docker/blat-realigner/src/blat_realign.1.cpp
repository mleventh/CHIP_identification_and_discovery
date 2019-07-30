#include "SeqLib/BLATWrapper.h"
class RealignmentFilter
{
    private:
		SeqLib::BLATWrapper* bw;
		SeqLib::BLATWrapper* setup_blat_server(const std::string& genome,const std::string& ooc){
			SeqLib::BLATWrapper * b = new SeqLib::BLATWrapper();
			SeqLib::BamHeader h("@HD\tVN:1.5\tSO:coordinate\n@SQ\tSN:1\tLN:249250621\n@SQ\tSN:2\tLN:243199373\n@SQ\tSN:3\tLN:198022430\n@SQ\tSN:4\tLN:191154276\n@SQ\tSN:5\tLN:180915260\n@SQ\tSN:6\tLN:171115067\n@SQ\tSN:7\tLN:159138663\n@SQ\tSN:8\tLN:146364022\n@SQ\tSN:9\tLN:141213431\n@SQ\tSN:10\tLN:135534747\n@SQ\tSN:11\tLN:135006516\n@SQ\tSN:12\tLN:133851895\n@SQ\tSN:13\tLN:115169878\n@SQ\tSN:14\tLN:107349540\n@SQ\tSN:15\tLN:102531392\n@SQ\tSN:16\tLN:90354753\n@SQ\tSN:17\tLN:81195210\n@SQ\tSN:18\tLN:78077248\n@SQ\tSN:19\tLN:59128983\n@SQ\tSN:20\tLN:63025520\n@SQ\tSN:21\tLN:48129895\n@SQ\tSN:22\tLN:51304566\n@SQ\tSN:X\tLN:155270560\n@SQ\tSN:Y\tLN:59373566\n@SQ\tSN:MT\tLN:16569\n@SQ\tSN:GL000207.1\tLN:4262\n@SQ\tSN:GL000226.1\tLN:15008\n@SQ\tSN:GL000229.1\tLN:19913\n@SQ\tSN:GL000231.1\tLN:27386\n@SQ\tSN:GL000210.1\tLN:27682\n@SQ\tSN:GL000239.1\tLN:33824\n@SQ\tSN:GL000235.1\tLN:34474\n@SQ\tSN:GL000201.1\tLN:36148\n@SQ\tSN:GL000247.1\tLN:36422\n@SQ\tSN:GL000245.1\tLN:36651\n@SQ\tSN:GL000197.1\tLN:37175\n@SQ\tSN:GL000203.1\tLN:37498\n@SQ\tSN:GL000246.1\tLN:38154\n@SQ\tSN:GL000249.1\tLN:38502\n@SQ\tSN:GL000196.1\tLN:38914\n@SQ\tSN:GL000248.1\tLN:39786\n@SQ\tSN:GL000244.1\tLN:39929\n@SQ\tSN:GL000238.1\tLN:39939\n@SQ\tSN:GL000202.1\tLN:40103\n@SQ\tSN:GL000234.1\tLN:40531\n@SQ\tSN:GL000232.1\tLN:40652\n@SQ\tSN:GL000206.1\tLN:41001\n@SQ\tSN:GL000240.1\tLN:41933\n@SQ\tSN:GL000236.1\tLN:41934\n@SQ\tSN:GL000241.1\tLN:42152\n@SQ\tSN:GL000243.1\tLN:43341\n@SQ\tSN:GL000242.1\tLN:43523\n@SQ\tSN:GL000230.1\tLN:43691\n@SQ\tSN:GL000237.1\tLN:45867\n@SQ\tSN:GL000233.1\tLN:45941\n@SQ\tSN:GL000204.1\tLN:81310\n@SQ\tSN:GL000198.1\tLN:90085\n@SQ\tSN:GL000208.1\tLN:92689\n@SQ\tSN:GL000191.1\tLN:106433\n@SQ\tSN:GL000227.1\tLN:128374\n@SQ\tSN:GL000228.1\tLN:129120\n@SQ\tSN:GL000214.1\tLN:137718\n@SQ\tSN:GL000221.1\tLN:155397\n@SQ\tSN:GL000209.1\tLN:159169\n@SQ\tSN:GL000218.1\tLN:161147\n@SQ\tSN:GL000220.1\tLN:161802\n@SQ\tSN:GL000213.1\tLN:164239\n@SQ\tSN:GL000211.1\tLN:166566\n@SQ\tSN:GL000199.1\tLN:169874\n@SQ\tSN:GL000217.1\tLN:172149\n@SQ\tSN:GL000216.1\tLN:172294\n@SQ\tSN:GL000215.1\tLN:172545\n@SQ\tSN:GL000205.1\tLN:174588\n@SQ\tSN:GL000219.1\tLN:179198\n@SQ\tSN:GL000224.1\tLN:179693\n@SQ\tSN:GL000223.1\tLN:180455\n@SQ\tSN:GL000195.1\tLN:182896\n@SQ\tSN:GL000212.1\tLN:186858\n@SQ\tSN:GL000222.1\tLN:186861\n@SQ\tSN:GL000200.1\tLN:187035\n@SQ\tSN:GL000193.1\tLN:189789\n@SQ\tSN:GL000194.1\tLN:191469\n@SQ\tSN:GL000225.1\tLN:211173\n@SQ\tSN:GL000192.1\tLN:547496\n@SQ\tSN:NC_007605\tLN:171823\n@SQ\tSN:hs37d5\tLN:35477943");
			
			b->SetHeaderInfo(h);
			b->LoadIndex(genome, ooc);
			return b;
		}
	
	public:

		RealignmentFilter(){	
		}
		RealignmentFilter(const std::string& genome,const std::string& ooc){
						bw = setup_blat_server(genome,ooc);
		}

		int blat_seq(const std::string& seq,const std::string& contig,const int mut_position){

			SeqLib::BamRecordVector blat_res;
			bw->AlignSequence("blat",seq,blat_res);

		        auto contains_mut=0;
		        auto does_not_contain_mut=0;

			for (auto& i : blat_res){

				if (i.Position() < mut_position && i.PositionEnd() >= mut_position )
					contains_mut = std::max(contains_mut,i.GetIntTag("AS"));
				else
					does_not_contain_mut = std::max(does_not_contain_mut,i.GetIntTag("AS"));
			}
			return contains_mut-does_not_contain_mut;
		}

};

