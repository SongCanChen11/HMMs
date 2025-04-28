import pandas as pd
import os
import re
import sys
import argparse


def get_best_hit_hmmscan(hmmscan_file):
    best_hit_hmmscan={}
    for l in open(hmmscan_file):
        if not l.startswith('#'):
            a=re.split(r'\s+',l.strip())
            gene, seqid, score = a[0],a[2],a[5]
            if not seqid in best_hit_hmmscan:
                best_hit_hmmscan[seqid]=re.sub(r'\.aln',r'',gene)
    return best_hit_hmmscan




parser = argparse.ArgumentParser(description='Parse the hmmscan results')
parser.add_argument('-i','--hmmscan_file', dest='hmmscan_file',required=True,
                    help='hmmscan results')
parser.add_argument('-f','--hmm_info_file', dest='hmm_info_file',required=True,
                    help='HMM info')
parser.add_argument('-o','--anno_outfile', dest='anno_outfile',required=True,
                    help='Output')

args = parser.parse_args()
hmmscan_file=args.hmmscan_file
hmm_info_file=args.hmm_info_file
anno_outfile=args.anno_outfile




# Read the information about sulfur-cycling marker genes
sulfur_cycling_marker_genes= ['ttrA','shyB','phsA','otr','npsr','mccA','dsrA','aprA','dsrB','asrA','sudA','fccB','tsdA','sqr','soxB','sorT','sorA','soeA','sdo','doxD','doxDA','yihQ','tmm','ssuD','sfnG','mtoX','msmA','mddA','dsyB','dsoB','dmsA','dmoA','dddY','dddX','dddQ','dddP','dddL','dddW','dddD','mtsA','sHdrB1','mmtN'] #'alma1','sor','tetH'# note that dddK is missing
hmm_info=pd.read_excel(hmm_info_file,index_col=[0,1]).reset_index()
marker_gene_hmm_info=hmm_info[hmm_info['Gene'].isin(sulfur_cycling_marker_genes)]


# Get the best hit from the hmmscan
best_hit_hmmscan=get_best_hit_hmmscan(hmmscan_file)
best_hit_hmmscan_df=pd.DataFrame([[key,value] for key, value in best_hit_hmmscan.items()],columns=['Accession','monophyletic_clade'])


# Export the annotation for sulfur-cycling marker genes
anno=best_hit_hmmscan_df.merge(marker_gene_hmm_info,left_on='monophyletic_clade',right_on='monophyletic_clade',how='inner')
export_df=anno[['Accession','monophyletic_clade','Gene','Product_Name','Category']].sort_values(by=['Category','Gene'])
export_df.to_csv(anno_outfile, header=True, sep = '\t', index=False)


