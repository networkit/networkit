#!/usr/bin/python2

import pickle as pickle
import os, sys, glob
from collections import defaultdict

# Depends on the python library networkx for
# exporting in GML and GraphML.
import networkx as nx

# Depends on scipy to parse matlab files.
import scipy.io

"""
A parser for the Facebook100 dataset.

Each of the school .mat files has an A matrix (sparse) and a
"local_info" variable, one row per node: ID, a student/faculty status
flag, gender, major, second major/minor (if applicable), dorm/house,
year, and high school. Missing data is coded 0.

See a bug? Contact me at conradlee@gmail.com
"""

node_attributes = ['second_major', 'gender', 'major_index', 'year', 'dorm', 'high_school', 'student_fac']

attribute_dict = {
    "student_fac" : 0,
    "gender" : 1,
    "major_index" : 2,
    "second_major" : 3,
    "dorm" : 4,
    "year" : 5,
    "high_school" : 6,
    }

def export_graph(G, write_filename):
    write_dir = "./output/" + write_filename + "/"
    if not os.path.isdir(write_dir):
        os.mkdir(write_dir)


    # Remove pho edge weights
    for n1 in G.edge:
        for n2 in G.edge[n1]:
            G.edge[n1][n2]={}
    print("\twriting gml")
    for node in G.nodes_iter():
        for key, val in list(G.node[node].items()):
            G.node[node][key]=int(val)
    nx.write_gml(G, write_dir + write_filename + ".gml")
    print("\twriting graphml")
    nx.write_graphml(G, write_dir + write_filename + ".graphml")
    print("\twriting edgelist")
    f = open(write_dir + write_filename + ".edgelist","w")
    for edge in G.edges_iter():
        f.write("\t".join([str(end) for end in list(edge)[:2]])+"\n")
    f.close()
    f = open(write_dir + write_filename + ".nodelist","w")
    print("\twriting nodelist")
    f.write("\t".join(["node_id"] + node_attributes) + "\n")
    for node in G.nodes_iter():
        f.write("\t".join([str(node)] + [str(G.node[node][attribute]) for attribute in node_attributes]) + "\n")

def get_attribute_partition(matlab_object, attribute):
    attribute_rows = matlab_object["local_info"]
    # Indicies as defined in comment above
    try:
        index = attribute_dict[attribute]
    except KeyError:
        raise KeyError("Given attribute " + attribute + " is not a valid choice.\nValid choices include\n" + str(list(attribute_dict.keys(\
))))
    current_id = 0
    partition = defaultdict(set)
    for row in attribute_rows:
        if not(len(row) == 7):
            raise ValueError("Row " + str(current_id) + " has " + str(len(row)) + " rather than the expected 8 rows!")
        value = row[index]
        partition[value].add(current_id)
        current_id += 1
    return dict(partition)

if __name__ == '__main__':
    if not os.path.isdir("./output/"):
        os.mkdir("./output/")

    matlab_files = glob.glob("./*.mat")
    matlab_files = [fn for fn in matlab_files if not("schools.mat" in fn)]
    for matlab_filename in matlab_files:
        network_name = matlab_filename.strip(".").strip("/").split("/")[-1].split(".")[0]
        print("Now parsing " + network_name)
        matlab_object = scipy.io.loadmat(matlab_filename)
        scipy_sparse_graph = matlab_object["A"]
        G = nx.from_scipy_sparse_matrix (scipy_sparse_graph)
        print("\tsize of network %d" % G.size())
        print("\torder of network %d" % G.order())
        print("\tmax degree of network %d" % max(G.degree().values()))
        for attribute in attribute_dict:
            partition = get_attribute_partition(matlab_object, attribute)
            for value in partition:
                for node in partition[value]:
                    G.node[node][attribute] = value
        export_graph(G, network_name)
