from ete3 import add_face_to_node
from numpy.random import randint

from Algorithm.Alghorithm import *
def visualize_tree(root_node: Node):
    ete_tree = Tree()
    build_tree(ete_tree, root_node, root_node.weight if not root_node.is_leaf() else 0)

    ts = TreeStyle()
    ts.show_leaf_name =True
    ts.show_branch_length = True
    ts.branch_vertical_margin = 10
    nstyle = NodeStyle()
    ts.title.add_face(TextFace("UPGMA Tree", fsize=15), column=0)
    nstyle["shape"] = "sphere"
    nstyle["size"] = 10
    nstyle["fgcolor"] = "darkred"
    nstyle["hz_line_type"] = 1
    nstyle["hz_line_color"] = "#cccccc"
    for n in ete_tree.traverse():
        n.set_style(nstyle)

    ete_tree.show(tree_style=ts)

def build_tree(ete_node, node: Node, parent_height: float = 0):
    if node.is_leaf():
        ete_node.name = node.label
        ete_node.dist = parent_height
    else:
        current_height = node.weight
        left_child = ete_node.add_child()
        left_child.dist = current_height - (node.left.weight if not node.left.is_leaf() else 0)
        build_tree(left_child, node.left, current_height)

        right_child = ete_node.add_child()
        right_child.dist = current_height - (node.right.weight if not node.right.is_leaf() else 0)
        build_tree(right_child, node.right, current_height)