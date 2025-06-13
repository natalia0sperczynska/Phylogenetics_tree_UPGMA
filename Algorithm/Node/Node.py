from __future__ import annotations
from collections import deque
class Node:
    def __init__(self,label):
        self.parent : Node = None
        self.left :Node = None
        self.right :Node = None
        self.weight = 0.0
        self.label = label

    def __str__(self):
        if self.is_leaf():
            return f'<{self.label}>'
        else:
            return f'<{self.label} [[{self.weight}]]: ↙:{self.left.label}, ↘:{self.right.label}>'

    def __repr__(self):
        return self.__str__()

    def is_leaf(self)->bool:
        return self.left is None and self.right is None

    def get_leaves(self, leaves : list[Node] = None) -> list[Node]:
        if not leaves:
            leaves = []
        if self.is_leaf():
            leaves.append(self)
            return leaves
        leaves = self.right.get_leaves(leaves)
        leaves = self.left.get_leaves(leaves)
        return leaves

    def print(self):
        lvl=0
        q=deque()
        q.append((self,lvl))
        while len(q)>0:
            node, node_lvl = q.popleft()
            if node_lvl > lvl:
                lvl=node_lvl
                print()
            print(node, end="\t")
            if not node.is_leaf():
                q.append((node.left,node_lvl+1))
                q.append((node.right,node_lvl+1))

