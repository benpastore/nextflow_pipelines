import argparse


def get_contigs(fai) :
    
    splits = 10
    lines = ''
    with open(fai, 'r') as f : 
        for line in f : 
            info = line.strip().split()
            chrom = info[0]
            size = int(info[1])
            batch_size = round(size/splits)

            curr = 1
            for i in range(splits) : 
                if i == splits-1 : 
                    lines += f"{chrom}:{curr}-{size}\n"
                    curr = curr+batch_size
                else : 
                    lines += f"{chrom}:{curr}-{curr+batch_size}\n"
                    curr = curr+batch_size
    f.close()
    op = open('contigs.txt', 'w')
    op.write(lines)
    op.close()

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-fai", type = str, required = True)

    return parser.parse_args()

def main() : 

    args = get_args() 
    get_contigs(args.fai)

if __name__ == "__main__" : 

    main()