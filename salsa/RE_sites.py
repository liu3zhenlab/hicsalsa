
import re 
import argparse

def parse_fasta(fh):
    fa = {}
    current_short_name = None
    # Part 1: compile list of lines per sequence
    for ln in fh:
        if ln[0] == '>':
            # new name line; remember current sequence's short name
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            # append nucleotides to current sequence
            fa[current_short_name].append(ln.rstrip())
    # Part 2: join lists into strings
    for short_name, nuc_list in fa.iteritems():
        # join this sequence's lines into one long string
        fa[short_name] = ''.join(nuc_list)
    return fa

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--assembly",help="assembled contigs")
    parser.add_argument("-e","--enzyme",help="restriction enzyme")
    #parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format")
    #parser.add_argument("-d","--dir",help="output directory for results",default='out')
    args = parser.parse_args()
    f = parse_fasta(open(args.assembly,'r'))
    enzymes_input = args.enzyme.replace(' ','').split(',')
    final_enzymes = []
    for each in enzymes_input:
        if 'N' in each:
            final_enzymes.append(each.replace('N','G'))
            final_enzymes.append(each.replace('N','A'))
            final_enzymes.append(each.replace('N','T'))
            final_enzymes.append(each.replace('N','C'))
        else:
            final_enzymes.append(each)


    for key in f:
        
        id,seq = key, f[key]
        left_count = 0
        rigt_count = 0
        for enzyme in final_enzymes:
            pos  = [m.start(0) for m in re.finditer(enzyme,seq)]
         
            length = len(seq)    
            for each in pos:
            	if each < length/2:
            		left_count += 1
            	else:
            		rigt_count += 1

        # pos  = [m.start(0) for m in re.finditer(r"GA[ACTG]TC",seq)]
     
        # length = len(seq)

        # left_count = 0
        # rigt_count = 0
        # for each in pos:
        #     if each < length/2:
        #         left_count += 1
        #     else:
        #         rigt_count += 1

        print id, left_count, rigt_count

    

if __name__ == '__main__':
    main()
