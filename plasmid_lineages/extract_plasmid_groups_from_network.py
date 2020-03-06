# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:13:51 2017

@author: mcastro
"""
import click
import networkx as nx

def get_columns( infileName ):
    """Parses a tab-delimited file and stores it in a list of lists. 
    This allows to easily access the columns"""
    
    data= open(infileName,"r")
    lines= data.readlines()
    precolumns= [l.strip("\n") for l in lines]
    data.close()
    columns= [l.split("\t") for l in precolumns]
    return columns


@click.command()
@click.option("--replicon_sizes", help="Complete path to a table with replicon sizes. The first column must be plasmid identifiers")
@click.option("--input_network", help="Complete path to input network. The first two columns must be plasmid identifiers")
@click.option("--output_dir", help="Complete path to writable output directory")
@click.option("--output_prefix", help="A prefix with wich all output files will begin")
def extract_groups(replicon_sizes, input_network, output_dir, output_prefix):
    """Extract plasmid groups from network
    """
        # Check if path is ending with slash; otherwise, add it
    if not output_dir.endswith("/"):
        output_dir = output_dir+"/"
    
    # Load replicon sizes
    replicon_size_dict = get_columns(replicon_sizes)
    replicon_size_dict = {c[0].strip(): int(c[1]) for c in replicon_size_dict if '#' not in c[0]} 
    
    # Load network and make neccesary name changes to avoid conflicts
    input_columns = get_columns(input_network)
    network = nx.Graph()
    for i in range(len(input_columns)):
        nodeA = input_columns[i][0]
        if 'lcl' in nodeA:
            nodeA = nodeA.split('lcl|')[1]
        if '.' in nodeA:
            nodeA = nodeA.split('.')[0]
        
        nodeB = input_columns[i][1]
        if 'lcl' in nodeB:
            nodeB = nodeB.split('lcl|')[1]
        if '.' in nodeB:
            nodeB = nodeB.split('.')[0]
        
        if nodeA == nodeB: #Ignore self-associations
            continue
        
        network.add_edge(nodeA, nodeB)
    
    
    # Calculate complete network degree
    net_degree = network.degree()
    net_degree = {item[0]: item[1] for item in net_degree}
    #net_degree_list=[net_degree.values().count(x) for x in sorted(set(net_degree.values()))]
    max_degree = max(net_degree.values()) 
    
    # Extract connected components
    conn_components = nx.connected_components(network)
    conn_components = sorted(conn_components,key=len, reverse=True) #sort by len; first the largest connected component
    
    # Open outfiles handles and write headers when needed
    f=open(output_dir+output_prefix+"_info_by_group.txt","w")
    g=open(output_dir+output_prefix+"_info_by_plasmid.txt","w")
    neigh=open(output_dir+output_prefix+"_reference_neighbors.txt","w")
    f.write("#group\tplasmids_in_group\n")
    g.write("#group\tplasmid\tplasmid_len\n")
    
    # Initialize variables used in loop
    c=0
    reference={}
    plasmids_in_conn_components=[]
    # Iterate through each connected component and create a subgraph with it
    for component in conn_components:
        c+=1
        graph=network.subgraph(component)
        # Get & sort node degrees
        # Choose reference (Among the most connected plasmids, the largest)
        degrees = graph.degree()
        degrees = {item[0]: item[1] for item in degrees}
        max_degree = max(degrees.values())
        hubs = [n for n in degrees if degrees[n] == max_degree]
        hubs = {node: replicon_size_dict[node] for node in hubs}
        sortedNodes = sorted(hubs, key=hubs.get, reverse=True) #sort by max value first
        selected_ref = sortedNodes[0]
        reference.update({c: [selected_ref, hubs[selected_ref], max_degree]})
        
        # Get reference neighbors and write result to file
        neighs=graph.neighbors(selected_ref)
        neigh.write(str(c)+"\t")
        for item in neighs:
            neigh.write(item+" ")
        neigh.write("\n")
        
        # Update list of plasmids in connected components
        plasmids_in_conn_components.extend(graph.nodes())
        
        # Write list of members of current component to file
        comp_string=str(component)
        comp_string=comp_string.replace('[', '').replace(']','')            
        comp_string=comp_string.replace('set','')
        comp_string=comp_string.replace('(', '').replace(')','')
        comp_string=comp_string.replace('\'', '')
        comp_string=comp_string.replace('{','')
        comp_string=comp_string.replace('}','')
        f.write(str(c)+"\t"+comp_string+"\n")
        
        # Write features of each plasmid to file
        plasmids_in_curr_comp=comp_string.split(', ')
        plasmids_in_curr_comp=[pl.strip() for pl in plasmids_in_curr_comp]
        for curr_plasmid in plasmids_in_curr_comp:
            g.write(str(c)+"\t"+curr_plasmid+"\t"+str(replicon_size_dict[curr_plasmid])+"\n")
        
    
    # Write plasmid info to file for special cases
    for plasmid in replicon_size_dict:
        if plasmid not in network.nodes():
            # Plasmids not in blast results
            g.write("NA"+"\t"+plasmid+"\t"+str(replicon_size_dict[plasmid])+"\n")
        elif plasmid not in plasmids_in_conn_components:
            # Plasmids in blast results, but outside connected components
            g.write("OUT"+"\t"+plasmid+"\t"+str(replicon_size_dict[plasmid])+"\n")
    
    # Close outfiles
    g.close()
    f.close()
    neigh.close()
    
    # Open outfile and write header
    f=open(output_dir+output_prefix+"_references.txt","w")
    f.write("#group\treference\tref_len\tdegree\n")
    for group in reference:
        f.write(str(group)+"\t"+reference[group][0]+"\t"+str(reference[group][1])+"\t"+str(reference[group][2])+"\n")
    
    f.close()
    return


# Run function
extract_groups()

