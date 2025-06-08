
from Algorithm.Alghorithm import *
from Visulisation.Tree_visualisation import *

if __name__ == '__main__':
    # seg = [convert_to_sequence("ATTGCCATT"), convert_to_sequence("ATGGCCATT"), convert_to_sequence("ATCCATTTTT"), convert_to_sequence("ATCTTCTT"),
    #        convert_to_sequence("ACTGACC")]
    # matrix_similarit = msa_alghoritm(seg)
    # final_al = center_star(matrix_similarit, seg)
    # centerSeq = get_main_seq(matrix_similarit, seg)
    # msa_result = merge_alignments(final_al, centerSeq)
    # table = msa_to_table(msa_result)
    #
    # #upgma_matrix = upgma_algorithm(seg)
    #
    # print("Similarity Matrix:")
    # print(matrix_similarit)
    #
    # #print("UPGMA matrix :")
    # #print(upgma_matrix)
    #
    # print("\nFinal Alignments:")
    # print(final_al)
    # print(merge_alignments(final_al, centerSeq))
    #
    # print("\nMSA Result:")
    # print_msa_table(msa_result)
    # score = score_msa(msa_result)
    # print(f"Total MSA Score: {score}")
    # table1 = msa_to_table(msa_result)
    # table_without_gaps=remove_gaps(table1)
    # print(table_without_gaps)
    # distance_table=get_distances(table_without_gaps)
    #
    # print(distance_table)
    labels = ['A', 'B', 'C', 'D', 'E','F','G']
    distance_table=pd.DataFrame(np.zeros((7,7)), index=labels, columns=labels)
    distance_table.loc['A', 'B':'G'] = [19, 27, 8, 33, 18, 13]
    distance_table.loc['B', 'C':'G'] = [31, 18, 36, 1, 13]
    distance_table.loc['C', 'D':'G'] = [26, 41, 32, 29]
    distance_table.loc['D', 'E':'G'] = [31, 17, 14]
    distance_table.loc['E', 'F':'G'] = [35, 28]
    distance_table.loc['F', 'G'] = 12
    distance_table = distance_table + distance_table.T
    np.fill_diagonal(distance_table.values, 0)
    print(distance_table)

    min=get_smallest_distance(distance_table)
    position=get_position_of_smallest_distance(min,distance_table)
    depth=min/2

    print(min)
    print(position)

    new_table1 = new_distance_matrix(distance_table,position)
    print(new_table1)
    root_node = create_final_matrix(distance_table)

    print("TREE:")
    root_node.print()
    visualize_tree(root_node)
