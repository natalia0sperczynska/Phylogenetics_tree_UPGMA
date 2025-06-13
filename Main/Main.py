from Sequences.Sequences import convert_to_sequence
from Algorithm.Alghorithm import *
from Visulisation.Tree_visualisation import *

def get_user_input():
    print("\nUPGMA Tree Creator")
    print("1.Enter sequences manually")
    print("2.Enter distance matrix manually")
    print("3.Use example data")
    try:
        choice = int(input("Your choice (1-3): ").strip())
        if choice==1:
            sequences_manually()
        elif choice==2:
            distance_matrix_manually()
        elif choice==3:
            example_data()
        else:
            print("Invalid choice")
            return get_user_input()
    except ValueError:
        print("Please enter a number")
        return get_user_input()

def sequences_manually():
    sequences = []
    print("\nEnter sequences (one per line, blank line to finish):")
    while True:
        seq = input().strip().upper()
        if not seq:
            break
        if not all(c in 'ACGT' for c in seq):
            print("Error: Only A,C,G,T are allowed")
            continue
        sequences.append(convert_to_sequence(seq))
    if len(sequences) < 2:
        print("Error: Need at least 2 sequences")
        return sequences_manually()
    try:
        matrix_similarit = msa_alghoritm(sequences)
        final_al = center_star(matrix_similarit, sequences)
        centerSeq = get_main_seq(matrix_similarit, sequences)
        msa_result = merge_alignments(final_al, centerSeq)
        alignments = merge_alignments(final_al, centerSeq)
        msa = print_msa_table(msa_result)
        score = score_msa(msa_result)
        table1 = msa_to_table(msa_result)
        table_without_gaps = remove_gaps(table1)
        distance_table=get_distances(table_without_gaps)
        root_node = create_final_matrix(distance_table)
        final_result = results_sequences_text(sequences, matrix_similarit, final_al,
                                              alignments, score, distance_table,
                                              root_node)
        print(final_result)
        visualize_tree(root_node)
        save_results(final_result)
    except Exception as e:
        print(f"Error processing sequences: {str(e)}")

def distance_matrix_manually():
    print("\nEnter number of sequences: ")
    try:
        n = int(input())
        if n < 2:
            print("Error: Need at least 2 sequences")
            return distance_matrix_manually()
    except ValueError:
        print("Error: Please enter a number")
        return distance_matrix_manually()

    print("Enter labels (one per line):")
    labels = [input(f"Label {i+1}: ").strip() for i in range(n)]
    matrix = pd.DataFrame(np.zeros((n, n)), index=labels, columns=labels)
    print("\nEnter pairwise distances:")
    for i in range(n):
        for j in range(i + 1, n):
            while True:
                dist = input(f"Distance between {labels[i]} and {labels[j]}: ").strip()
                if not dist:
                    if i != j:
                        matrix.iloc[i, j] = matrix.iloc[j, i]
                    break
                try:
                    dist = float(dist)
                    if dist < 0:
                        print("Distance cannot be negative")
                        continue
                    matrix.iloc[i, j] = dist
                    matrix.iloc[j, i] = dist
                    break
                except ValueError:
                    print("Please enter a number")
    np.fill_diagonal(matrix.values, 0)
    try:
        root_node = create_final_matrix(matrix)
        final_result = results_sequences_distance_matrix(matrix, root_node)
        print(final_result)
        visualize_tree(root_node)
        save_results(final_result)
    except Exception as e:
        print(f"Error processing sequences: {str(e)}")
        return None


def show_results_distance_matrix(distance_table:DataFrame):
    root_node = create_final_matrix(distance_table)
    print("TREE:")
    root_node.print()
    visualize_tree(root_node)

def validate_input(matrix:DataFrame):
    if not (matrix.values.diagonal() == 0).all():
        print("Error: Diagonal must be 0")
        return None
    if not np.allclose(matrix, matrix.T):
        print("Error: Matrix must be symmetric")
        return None
    else:
        return matrix

def example_data():
    labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    distance_table = pd.DataFrame(np.zeros((7, 7)), index=labels, columns=labels)
    distance_table.loc['A', 'B':'G'] = [19, 27, 8, 33, 18, 13]
    distance_table.loc['B', 'C':'G'] = [31, 18, 36, 1, 13]
    distance_table.loc['C', 'D':'G'] = [26, 41, 32, 29]
    distance_table.loc['D', 'E':'G'] = [31, 17, 14]
    distance_table.loc['E', 'F':'G'] = [35, 28]
    distance_table.loc['F', 'G'] = 12
    distance_table = distance_table + distance_table.T
    np.fill_diagonal(distance_table.values, 0)
    print(distance_table)

    min = get_smallest_distance(distance_table)
    position = get_position_of_smallest_distance(min, distance_table)
    root_node = create_final_matrix(distance_table)

    print("TREE:")
    root_node.print()

    visualize_tree(root_node)
def results_sequences_text(sequences: list, matrix_similarit: pd.DataFrame,
                           final_al: pd.DataFrame, alignments, score,
                           table_without_gaps, root_node):
    results = [
        "*===* UPGMA Phylogenetic Tree Analysis Results *===*",
        "\n=== INPUT SEQUENCES ==="
    ]

    for i, seq in enumerate(sequences, 1):
        results.append(f"Sequence {i}: {str(seq)}")

    sections = [
        ("\n=== SIMILARITY MATRIX ===", matrix_similarit),
        ("\n=== FINAL ALIGNMENTS ===", final_al),
        ("\n=== MULTIPLE ALIGNMENTS ===", alignments),
        ("\n=== MSA SCORE ===", f"Total Score: {score}"),
        ("\n=== DISTANCE MATRIX ===", table_without_gaps),
        ("\n=== PHYLOGENETIC TREE ===", get_tree_text(root_node))
    ]

    for header, content in sections:
        results.append(header)
        results.append(str(content))

    return "\n".join(results)


def results_sequences_distance_matrix(distance_matrix,root_node):
    results = [
        "*===* UPGMA Phylogenetic Tree Analysis Results *===*",
        "\n=== INPUT DISTANCE MATRIX ==="
    ]

    sections = [
        ("\n=== DISTANCE MATRIX ===", distance_matrix),
        ("\n=== PHYLOGENETIC TREE ===", get_tree_text(root_node))
    ]
    for header, content in sections:
        results.append(header)
        results.append(str(content))
    return "\n".join(results)


def save_results(results_text: str, filename: str = "upgma_results.txt"):
    try:
        full_path = os.path.abspath(filename)
        with open(filename, 'w',encoding='utf-8') as f:
            f.write(results_text)
        print(f"\nResults successfully saved to {full_path}")
    except Exception as e:
        print(f"Error saving results: {str(e)}")

def get_tree_text(root_node):
    from io import StringIO
    import sys

    old_stdout = sys.stdout
    sys.stdout = buffer = StringIO()

    print("TREE:")
    root_node.print()

    sys.stdout = old_stdout
    return buffer.getvalue()

if __name__ == '__main__':
    import os
    print(f"\nCurrent working directory: {os.getcwd()}")
    get_user_input()
    # seg = [convert_to_sequence("ATTGCCATT"), convert_to_sequence("ATGGCCATT"), convert_to_sequence("ATCCATTTTT"), convert_to_sequence("ATCTTCTT"),
    #        convert_to_sequence("ACTGACC")]