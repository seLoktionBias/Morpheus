"""Minimal Newick reader/writer with tip pruning.

Only what the pipeline needs: read one tree, keep a subset of tips, collapse the
resulting single-child internal nodes, and write it back. Written here so tree
pruning does not depend on ete3/dendropy being installed.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Sequence, Set


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name: str = "", length: Optional[float] = None):
        self.name = name
        self.length = length
        self.children: List["Node"] = []
        self.parent: Optional["Node"] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def add(self, child: "Node") -> None:
        child.parent = self
        self.children.append(child)

    def leaves(self) -> List["Node"]:
        if self.is_leaf:
            return [self]
        out: List[Node] = []
        for c in self.children:
            out.extend(c.leaves())
        return out

    def leaf_names(self) -> List[str]:
        return [n.name for n in self.leaves()]


_TOKEN = re.compile(r"\s*([(),;])\s*|\s*([^(),;]+)\s*")


def parse(text: str) -> Node:
    """Parse a Newick string into a tree. Supports names, branch lengths, comments."""
    text = re.sub(r"\[[^\]]*\]", "", text).strip()
    if not text:
        raise ValueError("empty Newick string")

    pos = 0

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if pos < len(text) and text[pos] == "(":
            pos += 1
            while True:
                node.add(parse_node())
                if pos < len(text) and text[pos] == ",":
                    pos += 1
                    continue
                if pos < len(text) and text[pos] == ")":
                    pos += 1
                break
        start = pos
        while pos < len(text) and text[pos] not in "(),;":
            pos += 1
        label = text[start:pos].strip()
        if ":" in label:
            name, _, length = label.rpartition(":")
            node.name = name.strip().strip("'\"")
            try:
                node.length = float(length)
            except ValueError:
                node.length = None
        else:
            node.name = label.strip().strip("'\"")
        return node

    root = parse_node()
    return root


def write(node: Node, with_lengths: bool = True) -> str:
    def render(n: Node) -> str:
        if n.is_leaf:
            body = n.name
        else:
            body = "(" + ",".join(render(c) for c in n.children) + ")" + n.name
        if with_lengths and n.length is not None and n.parent is not None:
            body += f":{n.length:g}"
        return body

    return render(node) + ";"


def read_file(path) -> Node:
    with open(path) as fh:
        return parse(fh.read())


def prune(root: Node, keep: Iterable[str]) -> Optional[Node]:
    """Return a copy of the tree containing only `keep` tips.

    Internal nodes left with a single child are collapsed, summing branch
    lengths so distances between surviving tips are preserved.
    """
    keep = set(keep)

    def copy_keep(n: Node) -> Optional[Node]:
        if n.is_leaf:
            if n.name in keep:
                out = Node(n.name, n.length)
                return out
            return None
        kids = [c for c in (copy_keep(c) for c in n.children) if c is not None]
        if not kids:
            return None
        if len(kids) == 1:
            child = kids[0]
            if n.length is not None and child.length is not None:
                child.length += n.length
            elif n.length is not None:
                child.length = n.length
            return child
        out = Node(n.name, n.length)
        for k in kids:
            out.add(k)
        return out

    pruned = copy_keep(root)
    if pruned is not None:
        pruned.length = None
        pruned.parent = None
    return pruned


def prune_to_file(tree_path, keep: Iterable[str], out_path) -> Optional[List[str]]:
    """Prune a Newick file to `keep` and write it. Returns the surviving tips."""
    root = read_file(tree_path)
    pruned = prune(root, keep)
    if pruned is None or len(pruned.leaf_names()) < 2:
        return None
    from pathlib import Path
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(write(pruned) + "\n")
    return pruned.leaf_names()
