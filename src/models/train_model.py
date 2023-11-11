import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import pickle
import argparse
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, classification_report, f1_score,  log_loss, mean_squared_error, balanced_accuracy_score
from math import sqrt
import os
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import networkx as nx
from networkx.algorithms import community
from src.visualization.visualize import plot_network
import random
# Set the random seed for reproducibility
random.seed(42)

def matthews(y_true, y_pred):
    """
    Calculate the Matthews Correlation Coefficient and other metrics.
    """
    if type(y_true) == pd.Series:
        y_true = y_true.values

    P = len([x for x in y_true if x == 1])
    N = len([x for x in y_true if x == 0])

    Tp, Fp = 0, 0
    for i in range(len(y_true)):
        if y_true[i] == 1 and y_pred[i] == 1: Tp += 1
        elif y_true[i] == 0 and y_pred[i] == 1: Fp += 1

    Tn = N - Fp
    Fn = P - Tp

    try:
        mcc = (Tp * Tn - Fp * Fn) / sqrt(
            (Tn + Fn) * (Tn + Fp) * (Tp + Fn) * (Tp + Fp))
    except ZeroDivisionError:
        mcc = 0

    return (mcc, f" \n \
    P: {P:_} \n \
    Tp: {Tp:_} \n \
    Fp: {Fp:_} \n \
    N: {N:_} \n \
    Tn: {Tn:_} \n \
    Fn: {Fn:_}")


#TODO: change output name of reports and models to reflect current day, instead of being fixed names
parser = argparse.ArgumentParser(description="Train the model using processed data.")
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
# optional.add_argument("-d","--data", dest="data", action='store', default='data/processed/model_data_pharokka.csv', help="Path to the data to use in the model (default: data/processed/model_data_pharokka.csv).") #TODO: remove test data
optional.add_argument("-d","--data", dest="data", action='store', default='data/processed/model_data_test.csv', help="Path to the data to use in the model (default: data/processed/model_data_test.csv).")
optional.add_argument("-i","--distance", dest="distance", action='store', default=None, help="Path to distances of all phages procued by MASH (default: None).")
optional.add_argument("-o","--output_directory", dest="output", action='store', default='models/random_forest.pkl', help="Output path for the model (default: models/random_forest.pkl).")
optional.add_argument("-r","--report", dest="report", action='store', default='reports/training_random_forest.txt', help="Output path for the model (default: reports/training_random_forest.txt).")
optional.add_argument("-n","--network", dest="network", action="store_true", help="Create a cluster map with MASH distance. Blue clusters are positive and red means negative gram staining.")

args = parser.parse_args()

print("--- Initiating training of the model ---")

df = pd.read_csv(args.data)


#TODO: Clustering training data with MSH - output to a csv
print ("-- Clustering phages by DNA sequence with MASH...")
# -----------
# Obtaining sequences
genomes = {}

# Iterate over the DataFrame rows
for index, row in df.iterrows():
    genomes[row['id']] = row['sequence_inphared']

def unix_call(command):
    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.decode('utf-8')  # Decode stdout to string

def mash_dist_pair(file_1, file_2):
    command = f"mash dist {file_1} {file_2}"
    try:
        dist_output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        result = dist_output.decode('utf-8').strip()

        # Parse the output to extract shared hashes and total hashes
        parts = result.split('\t')  # Split the mash output into parts
        shared_hashes_info = parts[4]  # Get the shared hashes information, e.g., "1/1000"
        shared, total = map(int, shared_hashes_info.split('/'))  # Split into shared and total hashes
        identity = (shared / total) * 100 if total > 0 else 0  # Calculate identity percentage

        # Return the result with identity percentage added
        return (os.path.splitext(os.path.basename(file_1))[0],  # Genome1 ID
                os.path.splitext(os.path.basename(file_2))[0],  # Genome2 ID
                parts[2],  # Distance
                parts[3],  # P-value
                shared_hashes_info,  # Shared-hashes
                identity)  # Identity percentage

    except subprocess.CalledProcessError as e:
        error_message = f"Error: {e.output.decode()}"
        # Provide a formatted error tuple consistent with the successful output
        return (os.path.splitext(os.path.basename(file_1))[0],
                os.path.splitext(os.path.basename(file_2))[0],
                "error",  # Placeholder for Distance
                "error",  # Placeholder for P-value
                "error",  # Placeholder for Shared-hashes
                "error")  # Placeholder for Identity percentage
    
def cluster_genomes(genomes, create_sketch=True, batch_size = 500):
    temp_dir = '/home/alvaroguleon/temp_mash'
    sketch_files = []
    
    # Check if temp_dir already exists and has files
    if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
        # Check if the directory has sketch files or fasta files depending on create_sketch flag
        existing_files = os.listdir(temp_dir)
        for genome_id in genomes.keys():
            expected_file = f"{genome_id}.msh" if create_sketch else f"{genome_id}.fna"
            if expected_file in existing_files:
                sketch_files.append(os.path.join(temp_dir, expected_file))
            else:
                # If expected files are not found, they need to be created
                print(f"Expected file {expected_file} not found in {temp_dir}.")
                # Set flag to indicate file creation is required
                create_files = True
                break
        else:
            # Set flag to indicate no file creation is required if all files were found
            create_files = False
    else:
        # If the directory doesn't exist, create it and set flag to create files
        os.makedirs(temp_dir, exist_ok=True)
        create_files = True

    if create_files:
        print("Producing temporary files for MASH...")
        # print("FASTA files are being created")

        # Create FASTA or sketch files if needed
        for genome_id, sequence in genomes.items():
            file_path = os.path.join(temp_dir, f"{genome_id}.fna")
            with open(file_path, 'w') as file:
                file.write(f">{genome_id}\n{sequence}\n")
            
            if create_sketch:
                sketch_prefix = file_path.replace('.fna', '')  # Remove the '.fna' for mash output prefix
                # print(f"Sketch files are being created for {sketch_prefix} {file_path}")
                unix_call(f"mash sketch -o {sketch_prefix} {file_path}")
                sketch_files.append(sketch_prefix + '.msh')  # Append the .msh path
            else:
                sketch_files.append(file_path)  # Append the .fna path

    pairs = [(sketch_files[i], sketch_file) for i in range(len(sketch_files)) for sketch_file in sketch_files[i + 1:]]
    

    with ProcessPoolExecutor() as executor, open('data/interim/clustering/mash_distances.csv', 'w') as output_file:
        output_file.write('Genome1,Genome2,Distance,P-value,Shared-hashes,Identity\n')  # Update the header with the new column
        for batch_start in range(0, len(sketch_files), batch_size):
            batch_end = min(batch_start + batch_size, len(sketch_files))
            batch_pairs = [(sketch_files[i], sketch_file) for i in range(batch_start, batch_end) for sketch_file in sketch_files[i+1:batch_end]]
            
            futures = {executor.submit(mash_dist_pair, pair[0], pair[1]): pair for pair in batch_pairs}
            for future in as_completed(futures):
                file_1, file_2, distance, p_value, shared_hashes, identity = future.result()  # Update to receive identity
                if not distance.startswith("Error"):  # Only write if no error occurred
                    # Write including the new identity percentage
                    output_file.write(f"{file_1},{file_2},{distance},{p_value},{shared_hashes},{identity:.2f}\n")

            # Print the completion message for the batch
            print(f"Completed processing batch from {batch_start} to {batch_end - 1}")

    distance_df = pd.read_csv('data/interim/clustering/mash_distances.csv')
    shutil.rmtree(temp_dir)  # Clean up temporary directory
    print("**MASH clustering completed**")
    return distance_df


if args.distance is None:
    distances = cluster_genomes(genomes, create_sketch=True, batch_size = 10)
    
else:
    distances = pd.read_csv(args.distance)


#TODO: This has not been tested from here. I need to have my full data to be able to test it. It is either
# that or returning to using model_data with is a pain in the ass
# Create a new graph
G = nx.Graph()

# Add edges to the graph based on our distances
for index, row in distances.iterrows():
    distance = 1 - row['Identity'] / 100  # Convert percentage to a float value between 0 and 1
    G.add_edge(row['Genome1'], row['Genome2'], weight=distance)

# Verify that the graph is not empty
if G.number_of_edges() == 0:
    raise ValueError("The graph G has no edges. Check the 'mash_distances.csv' file for correct data.")

# If there are edges, print the number of edges and nodes
print(f"The graph G has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

# Proceed to detect communities only if the graph is not empty
if G.number_of_nodes() > 0:
    communities = community.greedy_modularity_communities(G, weight='weight')
else:
    print("Cannot find communities in an empty graph.")
    

# Creating a NetworkX graph map
if args.network:
    plot_network(df, communities, path_to_save='reports/figures/phage_network.png')



print ("-- Removing isolated clusters for validation...")

#TODO: perhaps remove communities based on distance and not my subjective network inspection
# Initialize the separate variables for storing removed communities
removed_communities = []

# Community numbers to remove (1-indexed as per your enumeration, so subtract 1 for 0-indexed Python lists)
communities_to_remove = [934 - 1, 469 - 1]  # -1 because lists are 0-indexed in Python

# Sort the list in reverse so removing by index doesn't affect the order of unvisited items
for community_index in sorted(communities_to_remove, reverse=True):
    # Remove the community and store it in the 'removed_communities' list if it exists
    try:
        removed_communities.append(communities.pop(community_index)) # pop method removes by index
    except IndexError as e:
        print(f"No community at index: {community_index + 1}")  # Add 1 to match your enumeration

print("\nRemoved communities:")
for i, community_set in enumerate(removed_communities, start=1):
    print(f"Removed Community {i}: {community_set}")



print ("-- Preparing training and test data...")

# Create a mapping dictionary from id code to community index
accession_to_community = {}
for i, community in enumerate(communities):
    for accession in community:
        accession_to_community[accession] = i

# Map the 'id' column to a new 'Community' column
df['Community'] = df['id'].map(accession_to_community)

# Shuffle the community indices to randomly select 80% for training
community_indices = list(range(len(communities)))
random.shuffle(community_indices)
train_community_count = int(0.8 * len(community_indices)) # 80% of data goes to training
train_community_indices = community_indices[:train_community_count]
test_community_indices = community_indices[train_community_count:]

# Now, we can split the dataframe based on these indices
train_df = df[df['Community'].isin(train_community_indices)]
test_df = df[df['Community'].isin(test_community_indices)]

### Checking the split
# Collect community indices in the training and test set
train_communities = set(train_df['Community'].unique())
test_communities = set(test_df['Community'].unique())

# Check for intersection
common_communities = train_communities.intersection(test_communities)

# If the intersection is empty, then the split is correct; otherwise, there's an issue
if common_communities:
    print(f"Error: Communities with indices {common_communities} appear in both training and test sets.")
else:
    print("Success: All communities are exclusively in either the training or the test set.")

# Set 'id' as the index for both train_df and test_df
train_df = train_df.set_index('id')
test_df = test_df.set_index('id')

# Drop the 'Community' column from both dataframes
train_df = train_df.drop(columns='Community')
test_df = test_df.drop(columns='Community')




# Features (independent variables)
#TODO: use real features and not the dummy
# features = ['genome_length_inphared', 'gc_%_inphared',
#         'cds_number_inphared', 'positive_strand_%_inphared',
#         'negative_strand_%_inphared', 'coding_capacity_inphared',
#         'molecule_type_inphared_DNA', 'topology_inphared_linear',
#         'jumbophage_inphared', 'topology_linear_inphared',
#         'topology_circular_inphared', 'molecule_inphared_type_ss-DNA',
#         'molecule_inphared_type_DNA', 'molecule_inphared_type_RNA',
#         'molecule_inphared_type_ss-RNA', 'length', 'gc_perc',
#         'transl_table', 'cds_coding_density', 'CARD_AMR_Genes',
#         'CDS', 'CRISPRs', 'VFDB_Virulence_Factors', 'connector',
#         'head_packaging', 'host_takeover', 'integration and excision', 'lysis',
#         'nucleotide_metabolism', 'other', 'tRNAs', 'tail', 'tmRNAs',
#         'transcription', 'unkown_function', 'frame_negative',
#         'frame_positive'] # removed sequence and complement sequence

features = [
    'genome_length_inphared', 'gc_%_inphared', 'cds_number_inphared', 'positive_strand_%_inphared', 
    'negative_strand_%_inphared', 'coding_capacity_inphared', 'tRNAs', 
    'molecule_inphared_type_DNA', 'molecule_inphared_type_RNA', 
    'molecule_inphared_type_ss-DNA', 'molecule_inphared_type_ss-RNA',
    'jumbophage_inphared', 'topology_circular_inphared',
    'topology_linear_inphared'
]

target = 'staining'

# Extract features and target
X_train = train_df[features]
y_train = train_df[target]
X_test = test_df[features]
y_test = test_df[target]



# Train a model
model = RandomForestClassifier()
model.fit(X_train, y_train)


# save model
pickle.dump(model, open(args.output, "wb"))


# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]  # Probabilities for the positive class

# Predict class labels for F1 score calculation
y_pred = model.predict(X_test)

# Encode labels for mean squared error calculation
label_encoder = LabelEncoder()
y_test_encoded = label_encoder.fit_transform(y_test)
y_pred_encoded = label_encoder.transform(y_pred)


# Training report --------
# Calculate F1 score
positive_label = 'positive'
f1score = f1_score(y_test, y_pred, pos_label=positive_label)

# Calculate log loss using predicted probabilities
logloss = log_loss(y_test, y_pred_proba)

# Calculate mean squared error using encoded labels
mse = mean_squared_error(y_test_encoded, y_pred_encoded)

# Calculate accuracy
accuracy = accuracy_score(y_test_encoded, y_pred_encoded)

# Calculate balanced accuracy
balanced_accuracy = balanced_accuracy_score(y_test_encoded, y_pred_encoded)

# Calculate Matthews Correlation Coefficient (MCC)
mcc, mcc_details = matthews(y_test_encoded, y_pred_encoded)

# Generate the classification report
classification_report_str = classification_report(y_test_encoded, y_pred_encoded)





# Start the training report string
report_str = "\n\n--- Training Report ---\n"

report_str += f"""
Classification Report:
{classification_report_str}

Accuracy: {accuracy}
Balanced Accuracy: {balanced_accuracy}
Log Loss: {logloss}
Mean Squared Error: {mse}
F1 Score: {f1score}
Matthews Correlation Coefficient: {mcc}
Matthews Correlation Coefficient Details: {mcc_details}
"""


## Validation report ----
# Flatten the removed_communities frozensets into a list
removed_accessions = [item for community in removed_communities for item in community]

# Filter the original dataframe to get the validation set based on Accession being in removed_accessions
validation_df = df[df['id'].isin(removed_accessions)].copy().set_index('id')

# Prepare features and true labels from the validation set
X_validation = validation_df[features]
y_validation_true = validation_df[target]

# Encode 'staining' labels to numerical values using the same encoder used during training
y_validation_encoded = label_encoder.transform(y_validation_true)

# Predict the staining for the validation set and get probabilities for both classes
y_validation_pred_proba = model.predict_proba(X_validation)

# Obtain the predicted class labels (as strings, which you already have)
y_validation_pred = model.predict(X_validation)

# Since y_validation_pred contains strings, we don't need to inverse_transform it again
validation_df['Predicted Staining'] = y_validation_pred  # Direct assignment

# Now validation_df will have the predicted staining, add probability columns based on predict_proba output
validation_df['Probability Negative'] = y_validation_pred_proba[:, 0]  # Probability of 'negative'
validation_df['Probability Positive'] = y_validation_pred_proba[:, 1]  # Probability of 'positive'

# Find the entries that the model predicted incorrectly
incorrect_predictions = validation_df[validation_df[target] != validation_df['Predicted Staining']]
# Check if the 'incorrect_predictions' dataframe is empty
if incorrect_predictions.empty:
    incorrect_pred_str = "There were no incorrect predictions: None"
else:
    # If it's not empty, create a string from the DataFrame
    incorrect_pred_str = f"Incorrect Predictions:\n{incorrect_predictions}"

# Evaluate the model on the validation set using the same metrics as before
f1score_validation = f1_score(y_validation_encoded, label_encoder.transform(y_validation_pred), pos_label=label_encoder.transform([positive_label])[0])
logloss_validation = log_loss(y_validation_encoded, y_validation_pred_proba)
mse_validation = mean_squared_error(y_validation_encoded, label_encoder.transform(y_validation_pred))
accuracy_validation = accuracy_score(y_validation_encoded, label_encoder.transform(y_validation_pred))
balanced_accuracy_validation = balanced_accuracy_score(y_validation_encoded, label_encoder.transform(y_validation_pred))
classification_report_validation = classification_report(y_validation_encoded, label_encoder.transform(y_validation_pred))

# Start the validation report string
validation_report_str = "\n\n--- Validation Report ---\n"

# Add validation metrics to the report string
validation_report_str += f"""
Classification Report:
{classification_report_validation}

Accuracy: {accuracy_validation}
Balanced Accuracy: {balanced_accuracy_validation}
Log Loss: {logloss_validation}
Mean Squared Error: {mse_validation}
F1 Score: {f1score_validation}
"""
# Append incorrect predictions string to the validation report
validation_report_str += f"\n{incorrect_pred_str}\n"

# Combine the training and validation reports
full_report_str = report_str + validation_report_str

# Save the report to a text file
with open(args.report, "w") as report_file:
    report_file.write(full_report_str)

print(full_report_str) 



print(f"The model was saved in {args.output}")
print(f"The full report was saved in {args.report}")