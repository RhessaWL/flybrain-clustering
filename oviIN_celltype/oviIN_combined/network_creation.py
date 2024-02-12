# Import necessary libraries
import neuprint as neu
from neuprint import Client, fetch_adjacencies, NeuronCriteria as NC, SynapseCriteria as SC
import pandas as pd
import argparse
import logging
import time

# Define  query parameters
parser = argparse.ArgumentParser()
parser.add_argument("neuron_list", help="List of neurons to query")
parser.add_argument('--auth_file', default='flybrain.auth', help='File with auth token for neuprint.janelia.org')
parser.add_argument('--auth_token', default='', help='Auth token for neuprint.janelia.org; overrides --auth_file')
parser.add_argument('--version', default='1.1', help="Hemibrain version (format as version number, with or without leading 'v').")

args = parser.parse_args()

# Variables needed
hemibrain_version = "v" + args.version.replace('v', '')
logging.info("Processing hemibrain data %s", hemibrain_version)


# Connect to Neuprint
client = neu.Client('neuprint.janelia.org')

if args.auth_token:
    auth_token = args.auth_token.strip()
else:
    auth_token_file = open(args.auth_file, 'r')
    auth_token = next(auth_token_file).strip()
np_client = Client('neuprint.janelia.org', dataset='hemibrain:' + hemibrain_version, token=auth_token)
logging.info("Connected to Neuprint")

# 1a. Importing upstream and downstream neurons by celltype
crit = NC(bodyId = args.neuron_list)
n_df_in, conn_df_in = fetch_adjacencies(crit, None) # weirdly this was named as oviIN being the "input"
n_df_o, conn_df_o = fetch_adjacencies(None, crit) # This was also named with this system (sorry for the overcomplicated confusion)

# Sorting neuron lists by weight
conn_df_in = conn_df_in.sort_values(by='weight', ascending= False)
conn_df_o = conn_df_o.sort_values(by='weight', ascending= False)
logging.info("Fetched data from Neuprint")

# b. Seperating into a list and appending from upstream and downstream dataframes
neuron_ids = conn_df_o.bodyId_pre._append(conn_df_in.bodyId_post).unique().tolist()

# Fetching celltype adjacencies from Neuprint
start_time = time.time()
conn_df, conn_df_ = fetch_adjacencies(neuron_ids, neuron_ids)
conn_in, conn_in_ = fetch_adjacencies(conn_df_o.bodyId_pre, conn_df_o.bodyId_pre)
conn_out, conn_out_ = fetch_adjacencies(conn_df_in.bodyId_post, conn_df_in.bodyId_post)
elapsed_time = time.time() - start_time
logging.info("Elapsed time: %.2f seconds", elapsed_time)
logging.info("Fetched data for each celltype")

#2b. This function creates an undirected list and accounts for bidirectionality
def create_undirected(df):
    undirected_edges = {}  # Dictionary to store the undirected edges and their weights

    for index, row in df.iterrows():
        source = row['type_pre']
        target = row['type_post']
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

#3. Create undirected edges list for all datasets
# Dataframe manipulation for easy merging of celltype specific data
conn_df['bodyId_pre'] = conn_df.bodyId.astype(int)
conn_df.rename(columns={'bodyId': 'bodyId_post'}, inplace=True)

# merging data frames
full_pre = conn_df_.merge(conn_df, on = 'bodyId_pre', how = 'right')
full = full_pre.merge(conn_df, on = 'bodyId_post', how = 'right')

# Data type conversion so bodyIDs are correct (probably excessive but I am paranoid)
full.bodyId_post = full.bodyId_post.astype('Int64')
full.bodyId_pre = full.bodyId_pre.astype('Int64')

# Create undirected edges list for full dataset
df=full[['type_pre','type_post', 'weight']]
undirected_full = create_undirected(df)

# Dataframe manipulation for easy merging of celltype specific data
conn_in['bodyId_pre'] = conn_in.bodyId.astype(int)
conn_in.rename(columns={'bodyId': 'bodyId_post'}, inplace=True)

# merging data frames
full_pre_in = conn_in_.merge(conn_in, on = 'bodyId_pre', how = 'right')
full_in = full_pre_in.merge(conn_in, on = 'bodyId_post', how = 'right')

# Data type conversion so bodyIDs are correct (probably excessive but I am paranoid)
full_in.bodyId_post = full_in.bodyId_post.astype('Int64')
full_in.bodyId_pre = full_in.bodyId_pre.astype('Int64')

# Create undirected edges list for input dataset
df_in=full_in[['type_pre','type_post', 'weight']]
undirected_in = create_undirected(df_in)

# Dataframe manipulation for easy merging of celltype specific data
conn_out['bodyId_pre'] = conn_out.bodyId.astype(int)
conn_out.rename(columns={'bodyId': 'bodyId_post'}, inplace=True)

# merging data frames
full_pre_out = conn_out_.merge(conn_out, on = 'bodyId_pre', how = 'right')
full_out = full_pre_out.merge(conn_out, on = 'bodyId_post', how = 'right')

# Data type conversion so bodyIDs are correct (probably excessive but I am paranoid)
full_out.bodyId_post = full_out.bodyId_post.astype('Int64')
full_out.bodyId_pre = full_out.bodyId_pre.astype('Int64')

# Create undirected edges list for output dataset
df_out=full_out[['type_pre','type_post', 'weight']]
undirected_out = create_undirected(df_out)

# Close the connection to Neuprint
client.close()