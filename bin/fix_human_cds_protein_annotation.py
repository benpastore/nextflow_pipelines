import re
import argparse

def extract_ENST_id(string):
    pattern = r'ENST\d+\.\d+'
    match = re.search(pattern, string)
    if match:
        return match.group(0)
    else:
        return None

def reformat(fasta, outname) : 

    lines = ''
    header = False
    seq = ''
    with open(fasta, 'r') as f : 
        for line in  f : 
            if line.startswith(">") : 
                if header : 
                    enst_id = extract_ENST_id(header)
                    lines += f'>{enst_id}\n{seq}\n'
                    seq = ''
                    header = line.strip()
                else : 
                    header = line.strip()
            else : 
                seq += line.strip()
        else :
            enst_id = extract_ENST_id(header)
            lines += f'>{enst_id}\n{seq}\n'
    f.close()

    op = open(outname, 'w') 
    op.write(lines)
    op.close()


def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--fasta", type = str, required = True)
    required.add_argument("-o", "--outname", type = str, required = False)
    return parser.parse_args()


def main() : 

    args = get_args()

    reformat(args.fasta, args.outname)

if __name__ == "__main__" : 

    main()


