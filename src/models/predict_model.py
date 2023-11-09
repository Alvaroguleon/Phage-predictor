import pickle
import argparse
from ..features.build_features import engineer_features
import pandas as pd
import os
import numpy as np

parser = argparse.ArgumentParser(description="Make a gram staining prediction for a given phage.")
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser.add_argument("-d", "--data", type=str, required=True,
                    help="Path to the data to be predicted. Accepted inputs are FASTA files for Pharokka inputs or a directory path to an already processed Pharokka output.")
optional.add_argument("-m","--model_path", dest="model", action='store', default="models/random_forest.pkl", 
                      help="Path to the trained model file (default: models/random_forest.pkl).")

optional.add_argument("-o","--output_directory", dest="output", action='store', default='reports/predictions.csv', help="Output path (default: 'reports/predictions.csv').")
args = parser.parse_args()

print("--- Initiating prediction of the input data ---")
fasta_extensions = ["fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"]


# Function to process a single FASTA file
def process_fasta(file_path):
    # Replace with your code to process the FASTA file
    print(f"Processing FASTA file: {file_path}")
    #TODO: PHAROKKA

# Check if the input data path is a directory or a file
if os.path.isdir(args.data):
    # The input is a directory; process all relevant files within
    print(f"Processing Pharokka output in directory: {args.data}")
    model_data = engineer_features(args.data)

elif os.path.isfile(args.data):
    # The input is a file; check if it's a FASTA file
    file_ext = os.path.splitext(args.data)[1][1:].lower()
    if file_ext in fasta_extensions:
        process_fasta(args.data)  # Call your processing function
    else:
        print(f"Unsupported file format: {file_ext}. Supported FASTA formats are: {', '.join(fasta_extensions)}")
        exit(1)

else:
    print(f"The provided data path does not exist: {args.data}")
    exit(1)


features = ['genome_length_inphared', 'gc_%_inphared',
        'cds_number_inphared', 'positive_strand_%_inphared',
        'negative_strand_%_inphared', 'coding_capacity_inphared',
        'molecule_type_inphared_DNA', 'topology_inphared_linear',
        'jumbophage_inphared', 'topology_linear_inphared',
        'topology_circular_inphared', 'molecule_inphared_type_ss-DNA',
        'molecule_inphared_type_DNA', 'molecule_inphared_type_RNA',
        'molecule_inphared_type_ss-RNA', 'length', 'gc_perc',
        'transl_table', 'cds_coding_density', 'CARD_AMR_Genes',
        'CDS', 'CRISPRs', 'VFDB_Virulence_Factors', 'connector',
        'head_packaging', 'host_takeover', 'integration and excision', 'lysis',
        'nucleotide_metabolism', 'other', 'tRNAs', 'tail', 'tmRNAs',
        'transcription', 'unkown_function', 'frame_negative',
        'frame_positive', 'sequence_inphared', 'reverse_complement_inphared']

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