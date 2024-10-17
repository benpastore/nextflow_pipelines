#!/usr/bin/env python3 

import glob
import os

def get_subdirectories(dir) : 

    return [ i[0] for i in os.walk(dir) ]

def combine_lanes(dir, outdir) : 

    if not os.path.exists(outdir) : 
        os.makedirs(outdir)

    subdirs = get_subdirectories(dir)
    for s in subdirs : 
        basename = os.path.basename(s)

        R1 = glob.glob(f"{s}/{basename}*_1*fq*gz")
        R2 = glob.glob(f"{s}/{basename}*_2*fq*gz")
        
        if len(R1) >= 1 : 
            print(R1)
            R1_sort = sorted(R1)
            R2_sort = sorted(R2)
            os.system(f"cat {' '.join(R1_sort)} > {outdir}/{basename}_1.fq.gz")
            os.system(f"cat {' '.join(R2_sort)} > {outdir}/{basename}_2.fq.gz")
            print(os.path.basename(s))
            print("++++++++++++++++++++++++++++")

        #if len(l1) > 0 : 
        #    for sample in l1 : 
        #        basename = os.path.basename(sample).split("_L1_")[0]
        #        R1 = glob.glob(f"{s}/{basename}*1.fq.gz")
        #        R2 = glob.glob(f"{s}/{basename}*2.fq.gz")
        #        
        #        R1_sort = sorted(R1)
        #        R2_sort = sorted(R2)

                #os.system(f"cat {' '.join(R1_sort)} > {outdir}/{basename}_combn_1.fq.gz")
                #os.system(f"cat {' '.join(R2_sort)} > {outdir}/{basename}_combn_2.fq.gz")

#combine_lanes("/fs/ess/PAS0631/kholman/deep_sequencing/new_sib_pairs_fq/RawData", "/fs/scratch/PAS0631/combined_fastq")

#combine_lanes("/fs/scratch/PAS0631/sibs_analysis/SibsRawData", "/fs/scratch/PAS0631/combined_fastq")

#combine_lanes("/fs/scratch/PAS0631/sibs_analysis/SibsRawData2ndRound", "/fs/scratch/PAS0631/combined_fastq")

#combine_lanes("/fs/scratch/PAS0631/192sibconcordant", "/fs/scratch/PAS0631/combined_fastq")

#combine_lanes("/fs/ess/PAS1473/deep_sequencing/2024_04_28_sequencing", "/fs/ess/PAS1473/deep_sequencing/2024_04_28_sequencing/combined_fastq")

#combine_lanes("/fs/scratch/PAS0631/recessive_modifier_fastq", "/fs/scratch/PAS0631/recessive_modifier_fastq/recessive_modifier_combined_fastq")
combine_lanes("/fs/ess/PAS1473/deep_sequencing/2024_08_25_Novogene/01.RawData", "/fs/ess/PAS1473/deep_sequencing/2024_08_25_Novogene/lanes_combined")