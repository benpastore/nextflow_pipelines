
import os
import sys
import argparse
import pandas as pd


class create_table() : 

    def __init__(self, files, annotation, outprefix, id) :
        self._files = files
        self._outprefix = outprefix
        self._annotation = annotation
        self._id = id if id is not None else 'gene_id'

    def cbind_tables(self) :
        
        for i,f in enumerate(self._files) : 

            sample_name = os.path.basename(f).replace(".counts.tsv", "")
            df = pd.read_table(f, sep = "\t", names = [f'{self._id}', f'{sample_name}'])

            if i == 0 : 
                df_out = df
            else : 
                df_out = df_out.merge(df, on = [f'{self._id}'], how = 'outer')
                df_out = df_out.fillna(0)

            self._df = df_out

    def merge_to_annotation(self) : 
        
        """
        annotation can be gene_id, blah, blah blah. 
        make sure is tab separated and has header, make sure one column is named gene_id
        """
        if self._annotation is not None : 
            annotation = pd.read_csv(self._annotation, sep = "\t")
            annotation = annotation.drop_duplicates()
            df_out = self._df.merge(annotation, on = [f'{self._id}'], how = 'left')
            df_out.to_csv(f"{self._outprefix}.tsv", header = True, index = False, sep = "\t")
        else : 
            self._df.to_csv(f"{self._outprefix}.tsv", header = True, index = False, sep = "\t")


def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)

    ### required 
    required = parser.add_argument_group('required arguments')
    required.add_argument("-f", "--files", type = str, required = True)
    required.add_argument("-a", "--annotation", type = str, required = False)
    required.add_argument("-o", "--outprefix", type = str, required = True)
    required.add_argument("-i", "--identifier", type = str, required = False)

    return parser.parse_args()


def main() : 

    args = get_args()

    files = args.files.replace("[","").replace("]","").split(", ")
    
    f = [ i for i in files if ".tsv" in i ]
    
    x = create_table(f, args.annotation, args.outprefix, args.identifier)
    x.cbind_tables()
    x.merge_to_annotation()

if __name__ == "__main__" : 

    main()

        

