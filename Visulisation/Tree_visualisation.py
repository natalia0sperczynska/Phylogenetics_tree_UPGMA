
from Algorithm.Alghorithm import *
def visualize_tree(root_node: Node):
    """Convert your Node-based tree to ete3 Tree and visualize it with correct distances."""
    ete_tree = Tree()
    build_tree(ete_tree, root_node, root_node.weight if not root_node.is_leaf() else 0)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.title.add_face(TextFace("UPGMA Tree", fsize=20), column=0)
    ete_tree.show(tree_style=ts)


def build_tree(ete_node, your_node: Node, parent_height: float = 0):
    """Recursively build ete3 tree from your Node structure with correct distances."""
    if your_node.is_leaf():
        ete_node.name = your_node.label
        ete_node.dist = parent_height
    else:
        current_height = your_node.weigh
        left_child = ete_node.add_child()
        left_child.dist = current_height - (your_node.left.weight if not your_node.left.is_leaf() else 0)
        build_tree(left_child, your_node.left, current_height)

        right_child = ete_node.add_child()
        right_child.dist = current_height - (your_node.right.weight if not your_node.right.is_leaf() else 0)
        build_tree(right_child, your_node.right, current_height)