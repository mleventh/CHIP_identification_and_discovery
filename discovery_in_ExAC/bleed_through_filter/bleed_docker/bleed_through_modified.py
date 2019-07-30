from argparse import ArgumentParser

def filter_bleed_through(IN_FILE, PASSED_FILE, REJECTED_FILE):
    count_rejected = 0
    count_passed = 0

    with open(IN_FILE,'rb') as reader, open(PASSED_FILE, 'w') as passed_writer, open(REJECTED_FILE, 'w') as rejected_writer:
        for line in reader:
            if not line.startswith("#"):
                values = line.strip("\n").split("\t")
                if 'contig' in values:                
                    T_ALT_COUNT_IDX = values.index('t_alt_count')
                    T_REF_COUNT_IDX = values.index('t_ref_count')
                    T_ALT_SUM_IDX = values.index('t_alt_sum')
                    T_REF_SUM_IDX = values.index('t_ref_sum')
                    ALT_IDX = values.index('alt_allele')
                    REF_IDX = values.index('ref_allele')
                    #VARIANT_TYPE_IDX = values.index('Variant_Type')
                    passed_writer.write('\t'.join(values + ['bleed_through_stat']) + '\n')
                    rejected_writer.write('\t'.join(values + ['bleed_through_stat']) + '\n')

                elif not line.startswith('#'): 

                    if values[ALT_IDX] != '-' and values[ALT_IDX] != '-':
                           
                        if float(values[T_ALT_COUNT_IDX]) != 0 and float(values[T_REF_COUNT_IDX]) != 0:                    
                            alt_data = float(values[T_ALT_SUM_IDX])/float(values[T_ALT_COUNT_IDX])
                            ref_data = float(values[T_REF_SUM_IDX])/float(values[T_REF_COUNT_IDX])
                            if ref_data - alt_data > 7:
                                rejected_writer.write('\t'.join(values + [str(ref_data-alt_data)]) + '\n')
                                count_rejected += 1
                            else:
                                passed_writer.write('\t'.join(values + [str(ref_data-alt_data)]) + '\n')
                                count_passed += 1
                    else:
                        passed_writer.write('\t'.join(values + ['']) + '\n')
    count_rejected_file = 'count_rejected.txt'
    count_passed_file = 'count_passed.txt'
    with open(count_passed_file, 'w') as writer:
        writer.write(str(count_passed))
    with open(count_rejected_file, 'w') as writer:
        writer.write(str(count_rejected))


if __name__ == '__main__':
    parser = ArgumentParser(description='Filter Bleed through mutations')
    parser.add_argument('IN_FILE',help='input File (Mutation Validated MAF)')
    parser.add_argument('PASSED_FILE',help='passed mutations File (MAF)')
    parser.add_argument('REJECTED_FILE',help='rejected mutations File (MAF)')
    args = parser.parse_args()
    if(args):
        print 'Picked up   IN_FILE  : ' + str(args.IN_FILE)
        print 'Picked up   PASSED_FILE  : ' + str(args.PASSED_FILE)
        print 'Picked up   REJECTED_FILE  : ' + str(args.REJECTED_FILE)
        filter_bleed_through(args.IN_FILE, args.PASSED_FILE, args.REJECTED_FILE)
    else:
        parser.print_help() 