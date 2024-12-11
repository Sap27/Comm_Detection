import networkx as nx
from itertools import combinations
from queue import Queue
import threading
from collections import deque, defaultdict

class Pick:
    def __init__(self, score, child, parent):
        self.score = score
        self.child = child
        self.parent = parent

class Picker(threading.Thread):
    def __init__(self, from_idx, to_idx, array, unprocessed, processed, out, lock):
        threading.Thread.__init__(self)
        self.from_idx = from_idx
        self.to_idx = to_idx
        self.array = array
        self.unprocessed = unprocessed
        self.processed = processed
        self.out = out
        self.lock = lock

    def run(self):
        for i in range(self.from_idx, self.to_idx):
            max_score = -1.0
            childe = -1
            parent = -1
            val = self.processed[i]
            with self.lock:
                for child in self.unprocessed:
                    if self.array[val][child] > max_score:
                        max_score = self.array[val][child]
                        childe = child
                        parent = val
            self.out.put(Pick(max_score, childe, parent))

class TreeNode:
    def __init__(self, id, min_module_size, parent=None, d=float('inf'), size=1):
        self.id = id
        self.parent = parent
        self.children = set()
        self.size = size
        self.d = d
        self.fragment_size = 0.0
        self.min_module_size = min_module_size
        self.ext_id = id

    def collect_member_ids(self):
        result = []
        q = deque([self])
        while q:
            n = q.popleft()
            result.append(n.ext_id)
            q.extend(n.children)
        return result

    def grow(self, new_node_id, d):
        self.size += 1
        if self.parent:
            self.parent.update_size(1.0)
        child = TreeNode(new_node_id, self.min_module_size, parent=self, d=d, size=1)
        self.children.add(child)
        return child

    def update_size(self, change):
        self.size += change
        if self.parent:
            return self.parent.update_size(change)
        else:
            self.fragment_size += change
            return self.size

    def find_cut_point(self):
        children = []
        scores = []
        q = deque([self])
        while q:
            c = q.popleft()
            q.extend(c.children)
            p = c.parent
            if not p:
                continue
            above_size = self.fragment_size - c.size
            if above_size >= 3.0 or c.size >= 3.0:
                children.append(c)
                if above_size!=0 and c.size!=0:
                    scores.append(c.d / min(above_size, c.size))
        if not scores:
            return None
        min_score_idx = scores.index(min(scores))
        return children[min_score_idx]

    def cut(self):
        self.fragment_size = self.size
        if self.parent:
            self.parent.update_size(-self.size)
            self.parent.children.remove(self)
            self.parent = None
        self.d = float('inf')

def read_input(file):
    graph = nx.Graph()
    with open(file, 'r') as f:
        for line in f:
            node1, node2, weight = line.strip().split('\t')
            weight = float(weight)
            graph.add_edge(node1, node2, weight=weight)
    return graph

def get_jaccard(graph, node1, node2):
    neighbors1 = set(graph.neighbors(node1))
    neighbors2 = set(graph.neighbors(node2))
    intersection = len(neighbors1 & neighbors2)
    union = len(neighbors1 | neighbors2)
    return intersection / union if union != 0 else 0

def run_dcut(threads, max_module_size, min_module_size,graph=None,directed=False):
    
    nodes = list(graph.nodes())
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    idx_to_node = {idx: node for idx, node in enumerate(nodes)}
    N = len(nodes)

    adjacency_matrix = [[0] * N for _ in range(N)]
    for node1, node2 in combinations(nodes, 2):
            i, j = node_to_idx[node1], node_to_idx[node2]
            weight = get_jaccard(graph, node1, node2)
            adjacency_matrix[i][j] = adjacency_matrix[j][i] = weight

    unprocessed = list(range(1, N))
    processed = [0]
    tree = {0: TreeNode(0, min_module_size)}
    out = Queue()
    lock = threading.Lock()

    while unprocessed:
        to_update = list(range(len(processed)))
        N = len(to_update)
        L = N // threads if N % threads == 0 else N // (threads - 1)
        L = max(L, 1)
        threads_list = []
        for i in range(0, N, L):
            t = Picker(i, min(N, i + L), adjacency_matrix, unprocessed, processed, out, lock)
            t.start()
            threads_list.append(t)
        for t in threads_list:
            t.join()

        best_score = -1.0
        best_child = -1
        best_parent = -1
        while not out.empty():
            pick = out.get()
            if pick.score > best_score:
                best_score = pick.score
                best_child = pick.child
                best_parent = pick.parent

        if best_child != -1 and best_parent != -1:
            with lock:
                unprocessed.remove(best_child)
                processed.append(best_child)
            parent_node = tree[best_parent]
            new_node = parent_node.grow(best_child, best_score)
            tree[best_child] = new_node

    clusters = []
    queue = deque([tree[0]])
    while queue:
        current = queue.popleft()
        if current.size <= max_module_size:
            clusters.append([idx_to_node[i] for i in current.collect_member_ids()])
        else:
            cut_point = current.find_cut_point()
            if cut_point:
                cut_point.cut()
                queue.append(current)
                queue.append(cut_point)

    return [[int(j) for j in i] for i in clusters]

# Example usage

"""import networkx as nx
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from collections import deque, defaultdict
from queue import PriorityQueue

class Pick:
    def __init__(self, score, child, parent):
        self.score = score
        self.child = child
        self.parent = parent

    def __lt__(self, other):
        return self.score > other.score  # Higher scores have higher priority


class TreeNode:
    def __init__(self, id, min_module_size, parent=None, d=float('inf'), size=1):
        self.id = id
        self.parent = parent
        self.children = set()
        self.size = size
        self.d = d
        self.fragment_size = 0.0
        self.min_module_size = min_module_size
        self.ext_id = id

    def collect_member_ids(self):
        result = []
        q = deque([self])
        while q:
            n = q.popleft()
            result.append(n.ext_id)
            q.extend(n.children)
        return result

    def grow(self, new_node_id, d):
        self.size += 1
        if self.parent:
            self.parent.update_size(1.0)
        child = TreeNode(new_node_id, self.min_module_size, parent=self, d=d, size=1)
        self.children.add(child)
        return child

    def update_size(self, change):
        self.size += change
        if self.parent:
            return self.parent.update_size(change)
        else:
            self.fragment_size += change
            return self.size

    def find_cut_point(self):
        children = []
        scores = []
        q = deque([self])
        while q:
            c = q.popleft()
            q.extend(c.children)
            p = c.parent
            if not p:
                continue
            above_size = self.fragment_size - c.size
            if above_size >= 3.0 or c.size >= 3.0:
                children.append(c)
                if above_size != 0 and c.size != 0:
                    scores.append(c.d / min(above_size, c.size))
        if not scores:
            return None
        min_score_idx = scores.index(min(scores))
        return children[min_score_idx]

    def cut(self):
        self.fragment_size = self.size
        if self.parent:
            self.parent.update_size(-self.size)
            self.parent.children.remove(self)
            self.parent = None
        self.d = float('inf')

def get_jaccard(graph, node1, node2):
    neighbors1 = set(graph.neighbors(node1))
    neighbors2 = set(graph.neighbors(node2))
    intersection = len(neighbors1 & neighbors2)
    union = len(neighbors1 | neighbors2)
    return intersection / union if union != 0 else 0

def build_adjacency_matrix(graph, nodes, node_to_idx):
    N = len(nodes)
    adjacency_matrix = np.zeros((N, N), dtype=np.float32)

    for node1, node2 in combinations(nodes, 2):
        i, j = node_to_idx[node1], node_to_idx[node2]
        weight = get_jaccard(graph, node1, node2)
        adjacency_matrix[i][j] = adjacency_matrix[j][i] = weight

    return adjacency_matrix

def process_chunk(start_idx, end_idx, adjacency_matrix, unprocessed, processed, out):
    max_score = -1.0
    childe = -1
    parent = -1

    for i in range(start_idx, end_idx):
        val = processed[i]
        for child in unprocessed:
            if adjacency_matrix[val][child] > max_score:
                max_score = adjacency_matrix[val][child]
                childe = child
                parent = val

    if childe != -1 and parent != -1:
        out.put(Pick(max_score, childe, parent))

def run_dcut(threads, max_module_size, min_module_size, graph=None, directed=False):
    nodes = list(graph.nodes())
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    idx_to_node = {idx: node for idx, node in enumerate(nodes)}
    N = len(nodes)

    adjacency_matrix = build_adjacency_matrix(graph, nodes, node_to_idx)

    unprocessed = list(range(1, N))
    processed = [0]
    tree = {0: TreeNode(0, min_module_size)}
    out = PriorityQueue()

    while unprocessed:
        to_update = list(range(len(processed)))
        chunk_size = max(1, len(to_update) // threads)
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for i in range(0, len(to_update), chunk_size):
                futures.append(executor.submit(
                    process_chunk,
                    i, min(len(to_update), i + chunk_size),
                    adjacency_matrix, unprocessed, processed, out
                ))

            for future in futures:
                future.result()

        best_score, best_child, best_parent = -1.0, -1, -1
        while not out.empty():
            pick = out.get()
            if pick.score > best_score:
                best_score = pick.score
                best_child = pick.child
                best_parent = pick.parent

        if best_child != -1 and best_parent != -1:
            unprocessed.remove(best_child)
            processed.append(best_child)
            parent_node = tree[best_parent]
            new_node = parent_node.grow(best_child, best_score)
            tree[best_child] = new_node

    clusters = []
    queue = deque([tree[0]])
    while queue:
        current = queue.popleft()
        if current.size <= max_module_size:
            clusters.append([idx_to_node[i] for i in current.collect_member_ids()])
        else:
            cut_point = current.find_cut_point()
            if cut_point:
                cut_point.cut()
                queue.append(current)
                queue.append(cut_point)

    return [[int(j) for j in i] for i in clusters]
"""