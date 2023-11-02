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
                    help="Path to the data to be predicted.")
optional.add_argument("-m","--model_path", dest="model", action='store', default="models/random_forest.pkl", 
                      help="Path to the trained model file (default: models/random_forest.pkl).")

optional.add_argument("-o","--output_directory", dest="output", action='store', default='reports/predictions.csv', help="Output path (default: 'reports/predictions.csv').")
args = parser.parse_args()

print("--- Initiating prediction of the input data ---")

# Extract the file extension
file_ext = os.path.splitext(args.data)[1][1:] 

fasta_extensions = ["fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"]
genbank_extensions = ["gb", "gbk"]

# Processing the input data --------

#TODO: pharokka
# Process as a FASTA file
if file_ext in fasta_extensions:
    pass  # Replace PHAROKKA

 # Process as a GenBank file
elif file_ext in genbank_extensions:
    model_data = engineer_features(args.data)

else:
    supported_formats = ', '.join(fasta_extensions + genbank_extensions)
    print(f"Unsupported file format: {file_ext}. Supported formats are: {supported_formats}")
    exit(1)


features = ['genome_length', 'jumbophage', 'gc_%',
       'trna_count', 'cds_number', 'coding_capacity', 'positive_strand_%',
       'negative_strand_%', 'molecule_type_ss-DNA', 'molecule_type_DNA',
       'molecule_type_RNA', 'molecule_type_ss-RNA', 'topology_circular','topology_linear']
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
output_df.to_csv('args.output', index=False)

print()
print(f"The output has been saved in {args.output}. The first 5 entries are: ")
print(output_df.head())

# python -m src.models.predict_model -d data/interim/genbank_engineering/50_sequences.gb