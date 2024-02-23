"""Get the personal connectome of neuron(s) that are inputed by the user. """

# import the necessary libraries here
import pandas as pd
import numpy as np
import datetime
import sys

def log_msg(*args, out=sys.stdout, **kwargs):
    """Print message m with a timestamp if out is not None."""
    if out:
        print(datetime.datetime.now().strftime("%Y %m %d %H:%M:%S "), *args, **kwargs, file=out)

def get_connectome(main_neurons, exclude_main_neurons=False, connectome_scope='full', weight_threshold=1, connectome_by_type=False):
    """Get the personal connectome of neuron or neurons that are inputed by the user. 
    This function returns a connectome dataframe that contains the weighted connections between bodyIds. The synaptic weights 
    are collapsed across ROIs. This dataframe can be used to create a graph of the connectome in NetworkX using
    from_pandas_edgelist. However, the dataframe will need to be reformatted in order to run the clustering algorithms.

    Usage Notes: 
        - Please makes sure to have your Neuprint client already created before running this function.
        - Remember to set the outputs to different well-named variable for future use.
    
    main_neurons: can be a single bodyId, a list of bodyIds, or NeuronCriteria

    Options:
        - include the main neurons or not
        - input, output, or full connectome
        - weight threshold for the connection strengths to include in the connectome
    """
    
    from neuprint import fetch_adjacencies
    from neuprint import NeuronCriteria as NC
    from neuprint import fetch_neurons

    # the 1st df returns all the neurons involved in making the specified connections
    pre, pre_conns = fetch_adjacencies(None, main_neurons)
    post, post_conns = fetch_adjacencies(main_neurons, None)
    
    # it will now be necessary for main_neurons to be a list of bodyIds
    if not (isinstance(main_neurons, int) or isinstance(main_neurons, list)):
        main_neurons_df, roi_counts_df = fetch_neurons(main_neurons)
        main_neurons = main_neurons_df['bodyId'].tolist()

    if connectome_scope == 'input':
        log_msg(f'Getting input connectome')

        if exclude_main_neurons:
            log_msg(f'Excluding main neurons from the input connectome')
            # remove the main neurons from the pre
            pre = pre[~pre.bodyId.isin(main_neurons)]

        # get connections among neurons using the bodyIds from pre
        log_msg(f'Fetching connections among neurons using the bodyIds from pre')
        partners_, connectome = fetch_adjacencies(pre['bodyId'], pre['bodyId'])

    elif connectome_scope == 'output':
        log_msg(f'Getting output connectome')

        if exclude_main_neurons:
            log_msg(f'Excluding main neurons from the output connectome')
            # remove the main neurons from the post
            post = post[~post.bodyId.isin(main_neurons)]

        # get connections among neurons using the bodyIds from post
        log_msg(f'Fetching connections among neurons using the bodyIds from post')
        partners_, connectome = fetch_adjacencies(post['bodyId'], post['bodyId'])

    elif connectome_scope == 'full':
        log_msg(f'Getting full connectome')

        # combine unique pre and post bodyIds
        partners = pd.concat([pre['bodyId'], post['bodyId']]).unique()
        # turn it back into a series
        partners = pd.Series(partners)

        if exclude_main_neurons:
            log_msg(f'Excluding main neurons from the full connectome')
            # remove the main neurons from the partners
            partners = partners[~partners.isin(main_neurons)]

        # get connections among neurons using the bodyIds from partners
        log_msg(f'Fetching connections among neurons using the bodyIds from partners')
        partners_, connectome = fetch_adjacencies(partners, partners)

    # get rid of the ROI column and group bodyId_pre and bodyId_post by summing weights across ROIs
    connectome = connectome.groupby(['bodyId_pre', 'bodyId_post'], as_index=False)['weight'].sum()

    # if weight_threshold is specified, remove connections with weights less than the threshold
    log_msg(f'Checking weight threshold')
    if weight_threshold > 1:
        log_msg(f'Removing connections with weights less than {weight_threshold}')
        connectome = connectome[connectome['weight'] >= weight_threshold]

     # if connectome_by_type is specified, merge the type information into the connectome and grouby type
    if connectome_by_type:
        # merge type_pre information from partners_ into connectome and rename type column to type_pre
        connectome = connectome.merge(partners_[['bodyId','type']], left_on='bodyId_pre', right_on='bodyId').rename(columns={'type':'type_pre'})

        # merge type_post information from partners_ into connectome and rename type column to type_post
        connectome = connectome.merge(partners_[['bodyId','type']], left_on='bodyId_post', right_on='bodyId').rename(columns={'type':'type_post'})

        # group by type_pre and type_post and sum the weights
        connectome = connectome[['type_pre','type_post','weight']].groupby(['type_pre','type_post'], as_index=False).sum()

# This function creates an undirected list and accounts for bidirectionality
def create_undirected(df):
    """Get the undirected list of synapses from a connectome df. 
    This function returns a undirected list that contains the weighted connections between bodyIds checked for 
    bidirectionality. This dataframe can be used for modularity analysis but will need to be formatted to run 
    through the generalized modularity analysis using the format_edgelist written by A. Kunin.
    """
    undirected_edges = {}  # Dictionary to store the undirected edges and their weights
    for index, row in df.iterrows():
        source = row['bodyId_pre']
        target = row['bodyId_post']
        weight = row['weight']

        # Check if the edge already exists in the reverse
        if (target, source) in undirected_edges:
            # Update the weight of the existing edge
            undirected_edges[(target, source)] += weight
        else:
            # Add a new edge to dict
            undirected_edges[(source, target)] = weight

    # Create a DataFrame from the undirected edges dictionary
    undirected_edgelist = pd.DataFrame(list(undirected_edges.keys()), columns=['source', 'target'])
    undirected_edgelist['weight'] = list(undirected_edges.values())
    
    return undirected_edgelist