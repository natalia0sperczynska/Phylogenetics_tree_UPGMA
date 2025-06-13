
from __future__ import annotations
from Algorithm.Node.Node import *
import numpy as np
import pandas as pd
from pandas.core.interchange.dataframe_protocol import DataFrame
from Sequences.Sequences import SequenceUser
#steps:
# align and names
# compare sequence - pairwise seq alignment
# count the mismatches and record them in the table
# complete table by comparing all seq
# look for the lowest value to find the first group
# arithmetic mean step by step, construct three

def algorithm_implementation(seq1 : SequenceUser, seq2 :SequenceUser, match=1, gap = -2, mismatch=-1)->DataFrame:
    """Implement the Needleman-Wunsch global alignment algorithm."""
    score = 0
    seq1 = "-" + seq1.seq()
    seq2 = "-" + seq2.seq()
    row_seq=[label for label in seq1]
    col_seq =[label for label in seq2]
    df = pd.DataFrame(np.zeros((len(seq1), len(seq2))), index=row_seq, columns=col_seq)
    for i in range(len(seq2)):
        df.iloc[0, i] = i * gap
    for j in range(len(seq1)):
        df.iloc[j, 0] = j * gap

    for i in range(1,len(seq1)):
        for j in range(1, len(seq2)):
            df.iloc[i,j]=max(
                df.iloc[i-1,j] + gap,
                df.iloc[i, j - 1] + gap,
                df.iloc[i - 1, j - 1] + match if df.index[i] == df.columns[j] else df.iloc[i - 1, j - 1] + mismatch)
    return df

def traceback(df : pd.DataFrame, match=1, gap=-2, mismatch=-1):
    """Perform traceback through scoring matrix to find optimal alignments."""
    accumulator = []
    tracebackr(len(df.index)-1, len(df.columns) -1, '', '', df, accumulator, gap=gap, mismatch=mismatch, match=match)
    return accumulator

def tracebackr(i, j, align1, align2, df : pd.DataFrame,  accumulator : list, match=1, gap=-2, mismatch=-1):
    """Recursive helper function for traceback operation."""
    if i == 0 and  j == 0:
        accumulator.append((align1, align2))
    if i > 0 and df.iloc[i,j] == df.iloc[i-1,j]+gap:
        tracebackr(i-1, j, '-' + align1, df.index[i]+align2, df, accumulator, match, gap, mismatch)
    if j>0 and df.iloc[i,j] == df.iloc[i,j-1]+gap:
        tracebackr(i,j-1, df.columns[j]+align1, "-"+align2, df, accumulator,match, gap, mismatch)
    if i>0 and j>0 and (df.iloc[i,j] == df.iloc[i-1,j-1]+mismatch or df.iloc[i,j] == df.iloc[i-1,j-1]+match):
        tracebackr(i-1,j-1, df.columns[j]+align1, df.index[i]+align2, df, accumulator, match, gap, mismatch)

def get_score(df:pd.DataFrame)->int:
    """Get the final alignment score from the scoring matrix."""
    return df.iloc[-1,-1]

def msa_alghoritm(sequences, match=1, gap=-2, mismatch=-1):
    """Calculate similarity matrix between all sequence pairs."""
    labels = [seq.seq() if hasattr(seq, 'seq') else str(seq) for seq in sequences]
    similarity_matrix = pd.DataFrame(0, index=labels, columns=labels)

    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            if i >= j:
                continue
            df = algorithm_implementation(seq1, seq2, match, gap, mismatch)
            score = int(get_score(df))
            similarity_matrix.loc[labels[i], labels[j]] = score
            similarity_matrix.loc[labels[j], labels[i]] = score

    similarity_matrix['sum'] = similarity_matrix.sum(axis=1)
    return similarity_matrix

# def upgma_distance_matrix(sequences, match=1, gap=-2, mismatch=-1):
#     labels = [chr(65+i) for i in range(len(sequences))]
#     upgma_matrix_mismatch_matrix = pd.DataFrame(0, index=labels, columns=labels)
#
#     for i, seq1 in enumerate(sequences):
#         for j, seq2 in enumerate(sequences):
#             if i >= j:
#                 continue
#             df = algorithm_implementation(seq1, seq2, match, gap, mismatch)
#             aligment = traceback(df, match, gap, mismatch)
#             get_distance()
#             upgma_matrix.iloc[i,j] =
#             upgma_matrix.iloc[j,i] =
#
#     return upgma_matrix_mismatch_matrix

# def get_distance():
#     mismatch_score = 0
#
#     distance =
#
#     return distance


def center_star(similarity_matrix,sequences, match=1, gap=-2, mismatch=-1):
    """Align all sequences to the center sequence."""
    center_label = similarity_matrix['sum'].idxmax()
    center_seq = next(seq for seq in sequences if str(seq) == center_label)
    final_alignments = []
    for seq in sequences:
        if seq==center_seq:
            continue
        df=algorithm_implementation(seq, center_seq, match, gap, mismatch)
        alignment=traceback(df, match, gap, mismatch)[0]
        final_alignments.append(alignment)

    return final_alignments

def get_main_seq(similarity_matrix, sequences):
    """Get the center sequence with the highest similarity score."""
    center_label = similarity_matrix['sum'].idxmax()
    return next(seq for seq in sequences if str(seq) == center_label)

def merge_alignments(final_alignments, main_seq):
    """Merge all alignments into a multiple sequence alignment."""
    main_seq = list(main_seq.seq())
    msa = [main_seq]
    final_alignments = [([x for x in guide], [x for x in aligned]) for (guide, aligned) in final_alignments]

    for algn_index, (guide_seq, aligned_seq) in enumerate(final_alignments):
        i = 0
        while i < len(guide_seq):
            if guide_seq[i] == '-' and (i >= len(msa[0]) or msa[0][i] != '-'):
                add_gap(msa, i)
                add_gap_in_remaining_alignments(final_alignments[algn_index + 1:], i)
            i += 1
        msa.append(aligned_seq)
    return msa


def add_gap(msa: list[list[str]], i: int):
    for seq in msa:
        seq.insert(i, '-')

def add_gap_in_remaining_alignments(final_alignments: list[tuple[list[str], list[str]]], i: int):
    for _, aligned_seq in final_alignments:
        aligned_seq.insert(i, '-')


def percentage_for_all_matches(alignments: list[tuple[str, str]]):
    return [match_percentage(*alm) for alm in alignments]

def match_percentage(seq1:str,seq2:str):
    identical =  int(sum(np.array([ord(s) for s in seq1]) == np.array([ord(s) for s in seq2])))
    gaps = seq1.count("-")+seq2.count("-")
    return identical/len(seq1), gaps/len(seq1)

def print_msa_table(msa: list[list[str]]):
    for row in msa:
        print("".join(row))

def score_msa(msa: list[list[str]], match=1, mismatch=-1, gap=-2) -> int:
    """Score the MSA by summing all pairwise scores for each column."""
    total_score = 0
    num_seqs = len(msa)
    alignment_length = len(msa[0])

    for col in range(alignment_length):
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                a = msa[i][col]
                b = msa[j][col]
                if a == '-' or b == '-':
                    total_score += gap
                elif a == b:
                    total_score += match
                else:
                    total_score += mismatch
    return total_score

def msa_to_table(msa_results: list[list[str]]):
    table = pd.DataFrame(msa_results)
    return table

def remove_gaps(df:DataFrame) -> DataFrame:
    df=df.drop(columns=df.columns[(df == '-').any()])
    return df

def get_distances(df:DataFrame)->DataFrame:
    label = [ chr(65+i) for i in range(df.shape[0]) ]
    distance_matrix = pd.DataFrame(np.zeros((df.shape[0], df.shape[0])), index=label, columns=label)
    data=df.to_numpy()
    for i in range(data.shape[0]):
        diff = (data != data[i]).sum(axis=1)
        distance_matrix.iloc[i, :] = diff
    return distance_matrix

def get_smallest_distance(df:DataFrame):
    df_positive=df.replace(0,np.nan)
    return df_positive.min(skipna=True).min(skipna=True)

def get_position_of_smallest_distance(min:float,df:DataFrame):
    min_positions = np.where(df == min)
    for row, col in zip(*min_positions):
        return (df.index[row], df.columns[col])

def new_distance_matrix(df:DataFrame,position):
    first_seq,second_seq=position
    new_cluster_label = f"({first_seq}{second_seq})"
    remaining_seq=[label for label in df.index if label not in position]
    new_dist=[]
    for seq in remaining_seq:
        dist=(df.loc[first_seq,seq]+df.loc[second_seq,seq])/2
        new_dist.append(dist)

    new_labels=[new_cluster_label]+remaining_seq
    size=len(new_labels)

    new_matrix=pd.DataFrame(np.zeros((size,size)),index=new_labels,columns=new_labels)

    for i, label_i in enumerate(new_labels):
        for j, label_j in enumerate(new_labels):
            if i==j:
               new_matrix.iloc[i,j]=0
            elif label_i == new_cluster_label:
                    new_matrix.iloc[i,j]=new_dist[j-1]
            elif label_j == new_cluster_label:
                    new_matrix.iloc[i,j]=new_dist[i-1]
            else:
                new_matrix.iloc[i,j]=df.loc[label_i,label_j]
    return new_matrix

def calcualte_internal_node_depth(node : Node, initial_distance_matrix : DataFrame):
    depth = 0
    right_leaves = node.right.get_leaves()
    left_leaves = node.left.get_leaves()
    for right_leaf in right_leaves:
        for left_leaf in left_leaves:
            depth+=initial_distance_matrix.loc[left_leaf.label,right_leaf.label]
    depth/=len(right_leaves)*len(left_leaves)
    return depth/2

def create_final_matrix(distance_matrix:DataFrame):
    tree={k:Node(k) for k in  distance_matrix.columns}
    print(tree)
    current_matrix=distance_matrix.copy()
    print(current_matrix)
    while len(current_matrix)>1:
        min_distance=get_smallest_distance(current_matrix)
        label1, label2=get_position_of_smallest_distance(min_distance,current_matrix)
        new_label = f"({label1}{label2})"
        node = Node(new_label)
        node.left = tree[label1]
        node.right=tree[label2]
        tree[label1].parent = node
        tree[label2].parent = node
        node.weight = calcualte_internal_node_depth(node,distance_matrix)
        tree[new_label] = node
        current_matrix=new_distance_matrix(current_matrix,(label1, label2))
        #print(current_matrix)
    return node

if __name__ == '__main__':
    pass




