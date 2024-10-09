import argparse


def get_contigs(fai, target_contigs, n_splits) :
    
    if not n_splits : 
        splits = 1
    else : 
        splits = n_splits

    lines = ''
    with open(fai, 'r') as f : 
        for line in f : 
            info = line.strip().split()
            chrom = info[0]
            if chrom in target_contigs : 
                size = int(info[1])
                batch_size = round(size/splits)

                curr = 1
                for i in range(splits) : 
                    if i == splits-1 : 
                        lines += f"{chrom}:{curr}-{size}\n"
                        curr = curr+batch_size+1
                    else : 
                        lines += f"{chrom}:{curr}-{curr+batch_size}\n"
                        curr = curr+batch_size+1
    f.close()
    op = open('contigs.txt', 'w')
    op.write(lines)
    op.close()

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-fai", type = str, required = True)
    parser.add_argument("-target_contigs", type = str, required = True)
    parser.add_argument("-splits", type = int, required = False)

    return parser.parse_args()

def main() : 

    args = get_args() 
    get_contigs(args.fai, args.target_contigs.split(" "), args.splits)

if __name__ == "__main__" : 

    main()