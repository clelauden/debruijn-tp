#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
import textwrap
matplotlib.use("Agg")

__author__ = "Clemence Lauden"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Clemence Lauden"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Clemence Lauden"
__email__ = "clemence.lauden@hotmail.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()

def read_fastq(path):
    with open(path, "r") as file:
        for line in file:
            if "@" not in line and "+" not in line and "J" not in line:
                line = line.replace("\n", "")
                yield line

def cut_kmer(seq,k_length):
    for i in range(0, len(seq)-k_length+1, 1):
        k_mer = seq[i: i + k_length]
        yield k_mer

def build_kmer_dict(path, k_lenght):
    dict_kmer = {}
    for i in read_fastq(path):
        for k_mer in cut_kmer(i,k_lenght):
            if k_mer in dict_kmer:
                dict_kmer[k_mer]+=1
            else:
                dict_kmer[k_mer] = 1
    return dict_kmer

def build_graph(dict_kmer):
    digraph = nx.DiGraph()
    for u, v in dict_kmer.items():
        digraph.add_edge(u[:-1], u[1:], weight=v)
    return digraph


def remove_paths(digraph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if not delete_entry_node :
            path = path[1:]
        if not delete_sink_node :
            path = path[:-1]
        digraph.remove_nodes_from(path)
    return digraph

def std(data):
    std = statistics.stdev(data)
    return std

def select_best_path(digraph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):

    if std(weight_avg_list) > 0:
        del path_list[weight_avg_list.index(max(weight_avg_list))]
    elif std(path_length) > 0:
        del path_list[path_length.index(max(path_length))]
    else :
        del path_list[randint(0, len(path_list))]
    
    digraph = remove_paths(digraph, path_list, delete_entry_node, delete_sink_node)
    return digraph
        
def path_average_weight(digraph, path):
    mean_w = statistics.mean([d["weight"] for (u, v, d) in digraph.subgraph(path).edges(data=True)])
    return mean_w

def solve_bubble(digraph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths (digraph, ancestor_node, descendant_node))
    path_length = [len(path) for path in path_list]
    weight_avg_list = [path_average_weight(digraph, path) for path in path_list]
    digraph = select_best_path(digraph, path_list, path_length, weight_avg_list)
    return digraph

def simplify_bubbles(digraph):
    bubble = False
    for node in digraph.nodes():
        predecessors_list = [n for n in digraph.predecessors(node)]
        if len(predecessors_list) > 1:
            for i in range(0,len(predecessors_list)-1):
                ancestor_node = nx.lowest_common_ancestor(digraph,predecessors_list[i],predecessors_list[i+1])
                if ancestor_node != None:
                    bubble = True
                    break
            if bubble == True:
                break
    if bubble == True:
        digraph = simplify_bubbles(solve_bubble(digraph,ancestor_node,node))
    return digraph

def solve_entry_tips(digraph, starting_nodes):
    path_list = []
    path_length = []
    weight_avg = []
    for node in digraph.nodes:
        predecessors = list(digraph.predecessors(node))
        if len(predecessors) > 1:
            for start in starting_nodes:
                path_list += list(nx.all_simple_paths(digraph, start, node))
    if len(path_list) != 0:
        for path in path_list:
            path_length.append(len(path))
            weight_avg.append(path_average_weight(digraph, path))
        digraph = select_best_path(digraph, path_list, path_length, weight_avg, delete_entry_node=True, delete_sink_node=False)
    return digraph

def solve_out_tips(digraph, list_node_out):
    path_list = []
    path_length = []
    weight_avg = []
    for node in digraph.nodes:
        successors_list = list(digraph.successors(node))
        if len(successors_list) > 1:
            for node_out in list_node_out:
                path_list += list(nx.all_simple_paths(digraph, node, node_out))
    if len(path_list) != 0:
        for path in path_list:
            path_length.append(len(path))
            weight_avg.append(path_average_weight(digraph, path))
        digraph = select_best_path(digraph, path_list, path_length, weight_avg, delete_entry_node=False, delete_sink_node=True)
    return digraph

def get_starting_nodes(digraph):
    start_node = [n for n,d in digraph.in_degree() if d==0]
    return start_node

def get_sink_nodes(digraph):
    stop_node = [n for n,d in digraph.out_degree() if d==0]
    return stop_node

def get_contigs(digraph, start_node, stop_node):
    contig = []
    for u in start_node:
        for v in stop_node:
            if nx.has_path(digraph, u, v) == True:
                for path in nx.all_simple_paths(digraph, u, v):
                    seq = path[0]
                    for node in path[1:]:
                        seq += node[-1]
                    contig.append((seq, len(seq)))
    return contig

def save_contigs(contig, output_file):
    with open(output_file, "w") as file:
        for i in range(len(contig)):
            file.write(">contig_" + str(i) + " len=" + str(contig[i][1]) + "\n" + textwrap.fill((contig[i][0]), width=80) + "\n")

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    dict_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dict_kmer)
    start_n = get_starting_nodes(graph)
    stop_n = get_sink_nodes(graph)
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, start_n)
    graph = solve_out_tips(graph, stop_n)
    contigs = get_contigs(graph,start_n, stop_n)
    save_contigs(contigs, args.output_file)
    draw_graph(graph,'mygraph')
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()









