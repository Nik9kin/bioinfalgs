import itertools as it

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from tqdm.auto import tqdm


class Assembler:
    def __init__(self, k):
        """
        Genome assembler based on De Bruijn graph

        :param k: k-mer length
        """
        self.G = nx.MultiDiGraph()
        self.k = k

    def build_graph(self, reads, verbose=False):
        self.G = nx.MultiDiGraph()

        if verbose:
            reads = tqdm(reads)

        for read in reads:
            for i in range(len(read) - self.k + 1):
                kmer = str(read[i: i + self.k])
                prefix, suffix = kmer[:-1], kmer[1:]
                if not self.G.has_edge(prefix, suffix):
                    self.G.add_edge(
                        prefix, suffix,
                        key=kmer,
                        coverage=1.0
                    )
                else:
                    self.G[prefix][suffix][kmer]['coverage'] += 1.0

    def compactify(self, verbose=False):
        nodes_list = list(self.G.nodes)
        np.random.shuffle(nodes_list)

        if verbose:
            nodes_list = tqdm(nodes_list)

        for node in nodes_list:
            if self.G.in_degree(node) == self.G.out_degree(node) == 1:
                pred = next(self.G.predecessors(node))
                succ = next(self.G.successors(node))
                if pred != node != succ:
                    seq1 = next(iter(self.G[pred][node]))
                    seq2 = next(iter(self.G[node][succ]))
                    seq_new = seq1 + seq2[self.k - 1:]

                    cov1 = self.G[pred][node][seq1]['coverage']
                    cov2 = self.G[node][succ][seq2]['coverage']
                    len1, len2 = len(seq1), len(seq2)
                    cov_new = (cov1 * len1 + cov2 * len2) / (len1 + len2)

                    self.G.add_edge(pred, succ, key=seq_new, coverage=cov_new)
                    self.G.remove_node(node)

    def _del_isolates(self):
        self.G.remove_nodes_from(list(nx.isolates(self.G)))

    def cut_tails(self, factor=0.3):
        tails = []

        for pred, succ, seq, cov in self.G.edges(keys=True, data='coverage'):
            if (self.G.degree(pred) == 1) != (self.G.degree(succ) == 1):
                tails.append((len(seq) * cov, pred, succ, seq))

        if not tails:
            return

        tails = sorted(tails)
        max_len_cov = tails[-1][0]

        for len_cov, pred, succ, seq in tails:
            if len_cov < factor * max_len_cov:
                self.G.remove_edge(pred, succ, seq)

        self._del_isolates()

    def burst_bubbles(self, cov_thres=1.001):
        bubble_pairs = set()

        for pred, succ in self.G.edges(keys=False):
            if len(self.G[pred][succ]) > 1:
                bubble_pairs.add((pred, succ))

        for pred, succ in bubble_pairs:
            bubble_edges = []
            for seq, attrs in self.G[pred][succ].items():
                cov = attrs['coverage']
                if len(seq) == 2 * self.k - 1:
                    bubble_edges.append((cov, seq))

            if not bubble_edges:
                continue

            for cov, seq in bubble_edges:
                if cov <= cov_thres:
                    self.G.remove_edge(pred, succ, seq)

    def drop_low_covered_edges(self, cov_thres=1.0):
        low_covered_edges = []

        for pred, succ, seq, cov in self.G.edges(keys=True, data='coverage'):
            if cov <= cov_thres:
                low_covered_edges.append((pred, succ, seq))

        self.G.remove_edges_from(low_covered_edges)
        self._del_isolates()

    def local_resolve(self, avg_cov, conf=0.4):
        nodes_list = list(self.G.nodes)
        np.random.shuffle(nodes_list)

        nodes_with_selfloops = set(nx.nodes_with_selfloops(self.G))

        for node in nodes_list:
            if node not in nodes_with_selfloops:
                in_edges = list(self.G.in_edges(node, keys=True, data='coverage'))
                out_edges = list(self.G.out_edges(node, keys=True, data='coverage'))

                cov_sum_in = sum(cov for *_, cov in in_edges)
                cov_sum_out = sum(cov for *_, cov in out_edges)

                if ((len(in_edges) == 1 or len(out_edges) == 1) and
                        abs(cov_sum_in - cov_sum_out) < conf * avg_cov):
                    if len(in_edges) == 1:
                        in_edge = in_edges[0]
                        pred = in_edge[0]
                        seq_in = in_edge[2]
                        l_in = len(seq_in)

                        for _, succ, seq_out, cov_out in out_edges:
                            seq_new = seq_in + seq_out[self.k - 1:]

                            cov_in = in_edge[3] * (cov_out / cov_sum_out)

                            l_out = len(seq_out)
                            cov_new = (cov_in * l_in + cov_out * l_out) / (l_in + l_out)

                            self.G.add_edge(pred, succ, seq_new, coverage=cov_new)
                            if pred == succ:
                                nodes_with_selfloops.add(node)
                    else:
                        out_edge = out_edges[0]
                        succ = out_edge[1]
                        seq_out = out_edge[2]
                        l_out = len(seq_out)

                        for pred, _, seq_in, cov_in in in_edges:
                            seq_new = seq_in + seq_out[self.k - 1:]

                            cov_out = out_edge[3] * (cov_in / cov_sum_in)

                            l_in = len(seq_in)
                            cov_new = (cov_in * l_in + cov_out * l_out) / (l_in + l_out)

                            self.G.add_edge(pred, succ, seq_new, coverage=cov_new)
                            if pred == succ:
                                nodes_with_selfloops.add(node)

                    self.G.remove_node(node)

    def run(self,
            verbose=False):
        """
        Run standard pipeline with default params: `compactify`, `burst_bubbles`,
        `cut_tails`, `drop_low_covered_edges` and `compactify` again

        :param verbose: (default *False*)
        """
        self.compactify(verbose=verbose)
        self.burst_bubbles()
        self.cut_tails()
        self.drop_low_covered_edges()
        self.compactify(verbose=verbose)

    def get_node_degrees(self):
        return list(self.G.degree)

    def get_tails(self):
        tails = []

        for pred, succ, seq, cov in self.G.edges(keys=True, data='coverage'):
            if (self.G.degree(pred) == 1) != (self.G.degree(succ) == 1):
                tails.append((len(seq), cov, seq))

        return tails

    def get_edges(self):
        return [(cov, seq)
                for *_, seq, cov in self.G.edges(keys=True, data='coverage')]

    def get_contigs(self):
        # todo: should be modified with using euler paths in graph
        return sorted(self.get_edges(), reverse=True)

    def print_graph_size(self):
        print(f'Graph size: {self.G.number_of_nodes()} nodes '
              f'and {self.G.number_of_edges()} edges')

    def plot_graph(self,
                   subgraph=None,
                   ax=None,
                   edge_labels='auto',
                   show_node_labels='auto',
                   layout=nx.kamada_kawai_layout,
                   font_size=10,
                   figsize=(12, 12)):
        if subgraph is None:
            subgraph = self.G

        if ax is None:
            plt.figure(figsize=figsize)

        connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.30] * 4)]
        pos = layout(subgraph)

        nx.draw_networkx_nodes(subgraph, pos, ax=ax)

        if show_node_labels == 'auto':
            if self.k <= 9:
                nx.draw_networkx_labels(subgraph, pos, ax=ax)
        elif show_node_labels:
            nx.draw_networkx_labels(subgraph, pos, ax=ax)

        nx.draw_networkx_edges(
            subgraph, pos, edge_color="grey", connectionstyle=connectionstyle, ax=ax
        )

        if edge_labels == 'auto':
            labels = {
                tuple(edge): f"{self._short_view(edge[-1])}\ncov={cov:.2f}"
                for *edge, cov in subgraph.edges(keys=True, data='coverage')
            }
        elif edge_labels == 'coverage':
            labels = {
                tuple(edge): f"{cov:.2f}"
                for *edge, cov in subgraph.edges(keys=True, data='coverage')
            }
        else:
            raise ValueError(f"Unexpected value for edge_labels: {edge_labels}")

        nx.draw_networkx_edge_labels(
            subgraph,
            pos,
            labels,
            connectionstyle=connectionstyle,
            label_pos=0.3,
            font_color="blue",
            font_size=font_size,
            bbox={"alpha": 0},
            ax=ax,
        )

    def plot_graph_componentwise(self, *args, **kwargs):
        for nodes in nx.weakly_connected_components(self.G):
            subgraph = self.G.subgraph(nodes)
            self.plot_graph(subgraph, *args, **kwargs)
            plt.show()

    @staticmethod
    def _short_view(seq: str) -> str:
        if len(seq) <= 9:
            return str(seq)
        else:
            return f"<{len(seq)}bp>"


if __name__ == '__main__':
    from Bio.SeqIO import parse

    for filename, k in [
        # ('test/reads_0.fasta', 3),
        # ('test/reads_1.fasta', 9),
        # ('test/reads_2.fasta', 21),
        # ('test/reads_3.fasta', 41),
        ('test/reads_4.fasta', 91),
    ]:
        assembler = Assembler(k=k)

        with open(filename, 'r') as f_in:
            reads = (record.seq for record in parse(f_in, 'fasta'))
            assembler.build_graph(reads, verbose=True)

        assembler.print_graph_size()
        assembler.compactify(verbose=True)
        assembler.cut_tails()
        assembler.burst_bubbles()
        assembler.compactify()
        assembler.print_graph_size()

        edge_coverage = [cov for cov, seq in assembler.get_edges()]
        plt.hist(edge_coverage, bins=100)
        plt.show()

        assembler.drop_low_covered_edges(cov_thres=3.0)

        edge_coverage = [cov for cov, seq in assembler.get_edges()]
        plt.hist(edge_coverage, bins=100)
        plt.show()

        assembler.compactify()
        assembler.print_graph_size()
        assembler.plot_graph(figsize=(14, 5))
        plt.show()

        edge_coverage = [cov for cov, seq in assembler.get_edges()]
        plt.hist(edge_coverage, bins=100)
        plt.show()

        assembler.local_resolve(avg_cov=20.0)
        assembler.compactify()
        assembler.print_graph_size()
        assembler.plot_graph_componentwise(figsize=(14, 12))
