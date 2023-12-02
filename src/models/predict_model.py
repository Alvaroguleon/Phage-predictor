import pickle
import argparse
from ..features.build_features import engineer_features
import pandas as pd
import os
import numpy as np
import subprocess


parser = argparse.ArgumentParser(description="Make a gram staining prediction for a set of phages sequences.")
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser.add_argument("-d", "--data", type=str, required=True,
                    help="Path to the data to be predicted. Accepted inputs are FASTA files for Pharokka inputs or a directory path to an already processed Pharokka output.\
                        If the input is a FASTA file, the Pharokka output will be stored locally and predictions will be based on it.")
optional.add_argument("-m","--model_path", dest="model", action='store', default="models/random_forest.pkl", 
                      help="Path to the trained model file (default: models/random_forest.pkl).")
optional.add_argument("-p","--pharokka_directory", dest="pharokka", action='store', 
                          default="/mnt/c/Users/Alvaro/Desktop/pharokka_output/", help="Path where pharokka output will be stored (default: /mnt/c/Users/Alvaro/Desktop/pharokka_output/).")
optional.add_argument("-b","--database_directory", dest="database", action='store', 
                          default="/mnt/c/Users/Alvaro/Desktop/pharokka_database/", help="Path to the Pharokka database (default./mnt/c/Users/Alvaro/Desktop/pharokka_database/).")
optional.add_argument("-o","--output_directory", dest="output", action='store', default='reports/predictions.csv', help="Output path (default: 'reports/predictions.csv').")
args = parser.parse_args()

print("--- Initiating prediction of the input data ---")
fasta_extensions = ["fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"]


def unix_call(command, environment=None):
    if environment:  # If the environment is specified, prepend with 'conda run -n env_name'
        command = f"conda run -n {environment} {command}"
    subprocess.run(command, shell=True, check=True)  

def process_fasta(file_path, output_dir, database_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    print(f"Processing FASTA file: {file_path}")

    # Run Pharokka using the pharokkaENV environment
    unix_call(f'pharokka.py -i {file_path} -o {output_dir} -f -m -t 16 -d {database_dir}', environment='pharokkaENV')

    
# Check if the input data path is a directory or a file
if os.path.isdir(args.data):
    # The input is a directory; process all relevant files within
    print(f"Processing Pharokka output in directory: {args.data}")
    model_data = engineer_features(args.data)

elif os.path.isfile(args.data):
    # The input is a file; check if it's a FASTA file
    file_ext = os.path.splitext(args.data)[1][1:].lower()
    if file_ext in fasta_extensions:
        process_fasta(args.data, args.pharokka, args.database) 

        # The input is a directory; process all relevant files within
        print(f"Processing Pharokka output in directory: {args.pharokka}")
        model_data = engineer_features(args.pharokka)


    else:
        print(f"Unsupported file format: {file_ext}. Supported FASTA formats are: {', '.join(fasta_extensions)}")
        exit(1)

else:
    print(f"The provided data path does not exist: {args.data}")
    exit(1)


features = ['genome_length_inphared', 'gc_%_inphared',
       'cds_number_inphared', 'positive_strand_%_inphared',
       'negative_strand_%_inphared', 'coding_capacity_inphared',
       'tRNAs_inphared', 'molecule_inphared_type_DNA', 'topology_linear', 'jumbophage_inphared',
       'topology_circular', 'molecule_inphared_type_ss-DNA',
       'molecule_inphared_type_RNA', 'molecule_inphared_type_ss-RNA', 'length',
       'gc_perc', 'transl_table', 'cds_coding_density', 'CARD_AMR_Genes',
       'CDS', 'CRISPRs', 'VFDB_Virulence_Factors', 'connector',
       'head_packaging', 'host_takeover', 'integration and excision', 'lysis',
       'nucleotide_metabolism', 'other', 'tRNAs', 'tail', 'tmRNAs',
       'transcription', 'unkown_function', 'frame_positive', 'frame_negative']

model_data = model_data.set_index(model_data.columns[0])


model_data = model_data[features]


# Making the prediction -----

# Load the model using the provided path
rf = pickle.load(open(args.model, "rb"))

# Predictions
new_data_pred = rf.predict(model_data)

# Get the probabilities for the predicted class for each instance
probas = rf.predict_proba(model_data)
predicted_indices = np.argmax(probas, axis=1)  # Get index of max proba for each sample
new_data_pred_proba = [probas[i][predicted_indices[i]] for i in range(len(predicted_indices))]

# Prepare the output DataFrame
output_df = pd.DataFrame({
    'id': model_data.index,
    'prediction': new_data_pred,
    'prediction_probability': new_data_pred_proba
})

# Uncomment the following line if you want to save the output to a CSV
output_df.to_csv(args.output, index=False)

print()
print(f"The output has been saved in {args.output}. The first 5 entries are: ")
print(output_df.head())

# python -m src.models.predict_model -d data/interim/genbank_engineering/50_sequences.gb