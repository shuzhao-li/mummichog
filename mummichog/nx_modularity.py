import networkx as nx

def find_communities(G, method='louvain'):
    '''
    Wrapper of nx.community functions
    
    '''
    if method == 'louvain':
        return nx.community.louvain_communities(G)
    elif method == 'Clauset-Newman-Moore':
        return nx.community.greedy_modularity_communities(G)
    elif method == 'girvan_newman':
        # no optiaml, returning all levels
        coms = nx.community.girvan_newman(G)
        return [list(x) for x in coms]
    else:
        print("Warning, method not implemented")