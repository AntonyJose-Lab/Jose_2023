#! /usr/bin/env python
"""This program aims to represent entity-sensor-property systems with up/down regulation and to compute all possible heritable regulatory architectures, at least for n < 5. This can potentially be extended to compare topologies of networks reported in the literature to add missing nodes/links that are necessary for heredity. This could be useful for building regulatory networks from the ground up. Missing regulatory links and/or entities can be predicted through comparison with the non-isomoprphic set of graphs with heritable regulatory architectures."""

import matplotlib.pyplot as plt
import networkx as nx
import itertools
import math
from itertools import chain
from networkx.algorithms import isomorphism

non_iso_regulatory_architectures_2_and_3 = []
heritable_non_iso_regulatory_architectures_2_and_3 = []
headers = ["Entities","Heritable","Regulated","Heritably-regulated", "Bits"]
print(headers)

for entities in range (1,4):
    non_iso_regulatory_architectures = []
    heritable_regulatory_architectures = []
    heritable_non_iso_regulatory_architectures = []
    print_line = []
    all_nodes = list( range( 1, entities+1 ) )
    all_edges = list( itertools.permutations( all_nodes, r=2 ) )
    all_possible_edge_sets = []
    for x in range (0 , len(all_edges) + 1):
        all_possible_edge_sets.extend (list ( itertools.combinations( all_edges, r=x )))
    print_line.append(entities)
    ## Get the list of all weakly connected non-isomoprhic Digraphs with every node in-degree != 0.
    weakly_connected_heritable_non_iso_digraphs = []
    for edge_set in all_possible_edge_sets:
        g = nx.DiGraph()
        g.add_nodes_from( all_nodes )
        g.add_edges_from( edge_set )
        if nx.is_weakly_connected( g ):
            heritable = True
            for i in range ( 1, entities+1 ):
                if g.in_degree(i) == 0:
                    heritable = False
            if heritable == True:
                found = False
                for existing in weakly_connected_heritable_non_iso_digraphs:
                    DiGM = isomorphism.DiGraphMatcher(g, existing)
                    if DiGM.is_isomorphic():
                        found = True
                        break
                if not found:
                    weakly_connected_heritable_non_iso_digraphs.append( g )
    print_line.append(len(weakly_connected_heritable_non_iso_digraphs))
    ## Add regulation and then keep only those that satisfy heritability condition.
    for g in weakly_connected_heritable_non_iso_digraphs:
        edge_set = g.edges()
        ## Generate all possible regulated graphs for a given DiGraph architecture.
        regulated_edges = list(itertools.combinations_with_replacement(['gray','black'], len(edge_set))) ## "up" = 'gray' and "down" = 'black' to conveniently call colors for plotting.
        permuted_regulated_edges_list = []
        for m in range ( 0, len(regulated_edges)): ## Note index of a list is zero.
            permuted_regulated_edges_list.append(list(set(list(itertools.permutations(regulated_edges[m], len(edge_set))))))
        permuted_regulated_edges = list(itertools.chain(*permuted_regulated_edges_list))
        for n in range ( 0, len(permuted_regulated_edges) ):
            graph_count = 1
            regulated_g = nx.DiGraph()
            regulated_g.add_nodes_from( all_nodes )
            now_regulated_edges = list(permuted_regulated_edges[n])
            regulated_g.add_edges_from(edge_set)
            ## The topology of each regulated_g is the same as that of each Digraph g, but each regulated_g needs to use a different set of regulatory assignments.
            number_regulated = ['regulation']*len(edge_set)
            ## Need to create nested dictionary to properly assign attributes to edges.
            regulation_to_edges = {}
            count = 0
            for edge in edge_set:
                    regulation_to_edges[edge] = {}
                    regulation_to_edges[edge]['regulation'] = now_regulated_edges[count]
                    count = count + 1
            nx.set_edge_attributes(regulated_g, regulation_to_edges)

            # Do not add if isomorphic to previous regulatory architectures.
            found = False
            for iso in non_iso_regulatory_architectures:
                DiGM = isomorphism.DiGraphMatcher(regulated_g, iso, edge_match= lambda e1,e2: e1['regulation'] == e2['regulation'])
                if DiGM.is_isomorphic():
                    found = True
                    break
            if not found:
                non_iso_regulatory_architectures.append( regulated_g )

            ## Check for regulation being heritable, given the architecture. Rule: If any entity is only downregulated, such an architecture is not heritable.
            for o in range( 1, entities+1 ):
                entity_count = 0
                for p in range( 1, entities+1 ):
                    if regulated_g.has_edge(p,o):
                        if regulated_g.get_edge_data(p,o) == {'regulation': 'gray'}:
                            entity_count = entity_count + 1
                            #print(o, p, entity_count)
                if entity_count == 0:
                    graph_count = 0
            if graph_count == 1:
                heritable_regulatory_architectures.append ( regulated_g )

        # Remove isomorphic heritable_regulatory_architectures, if any.
        heritable_non_iso_regulatory_architectures = []
        for non_iso in heritable_regulatory_architectures:
            found = False
            for iso in heritable_non_iso_regulatory_architectures:
                DiGM = isomorphism.DiGraphMatcher(non_iso, iso, edge_match= lambda e1,e2: e1['regulation'] == e2['regulation'])
                if DiGM.is_isomorphic():
                    found = True
                    break
            if not found:
                heritable_non_iso_regulatory_architectures.append( non_iso )

    print_line.append(len(non_iso_regulatory_architectures))
    print_line.append(len(heritable_non_iso_regulatory_architectures))
    if len(heritable_non_iso_regulatory_architectures) == 0:
        print_line.append("na")
    if len(heritable_non_iso_regulatory_architectures) != 0:
        print_line.append(round(math.log2(len(heritable_non_iso_regulatory_architectures)), 2))
    print(print_line)
    ## Concatenate heritable architectures to plot.
    if entities == 2 or entities == 3: ## Since 1 has no architectures and 4 has too many to plot (19559 non-iso regulatory and 5604 heritable non-iso regulatory).
        non_iso_regulatory_architectures_2_and_3 = non_iso_regulatory_architectures_2_and_3 + non_iso_regulatory_architectures
        heritable_non_iso_regulatory_architectures_2_and_3 =  heritable_non_iso_regulatory_architectures_2_and_3 + heritable_non_iso_regulatory_architectures

# Draw all non-isomeric regulatory architectures.
columns = min( len( non_iso_regulatory_architectures_2_and_3 ), 10 )
rows = (len( non_iso_regulatory_architectures_2_and_3 ) + columns - 1) // columns

for i, g in enumerate( non_iso_regulatory_architectures_2_and_3 ):
    ax = plt.subplot( rows, columns, i+1 )
    ax.axis( 'off' )
    pos = nx.drawing.layout.circular_layout( g )
    node_color_map = []
    nodes_now = len(g.nodes())
    for k in range ( 1, nodes_now+1 ):
        if g.out_degree(k) == 0:
            node_color_map.append('blue')
        else:
            node_color_map.append('red')
    edge_color_map = []
    edges = g.edges()
    edge_color_map = [g[u][v]['regulation'] for u,v in edges]
    nx.drawing.nx_pylab.draw_networkx_nodes( g, node_size=30, node_color=node_color_map, pos=pos )
    nx.drawing.nx_pylab.draw_networkx_edges( g, node_size=30, connectionstyle="arc3,rad=0.2", edge_color=edge_color_map, pos=pos )
plt.show()

# Draw all non-isomeric heritable regulatory architectures.
columns = min( len( heritable_non_iso_regulatory_architectures_2_and_3 ), 5 )
rows = (len( heritable_non_iso_regulatory_architectures_2_and_3 ) + columns - 1) // columns

for i, g in enumerate( heritable_non_iso_regulatory_architectures_2_and_3 ):
    ax = plt.subplot( rows, columns, i+1 )
    ax.axis( 'off' )
    pos = nx.drawing.layout.circular_layout( g )
    node_color_map = []
    nodes_now = len(g.nodes())
    for k in range ( 1, nodes_now+1 ):
        if g.out_degree(k) == 0:
            node_color_map.append('blue')
        else:
            node_color_map.append('red')
    edge_color_map = []
    edges = g.edges()
    edge_color_map = [g[u][v]['regulation'] for u,v in edges]
    nx.drawing.nx_pylab.draw_networkx_nodes( g, node_size=30, node_color=node_color_map, pos=pos )
    nx.drawing.nx_pylab.draw_networkx_edges( g, node_size=30, connectionstyle="arc3,rad=0.2", edge_color=edge_color_map, pos=pos )
plt.show()
