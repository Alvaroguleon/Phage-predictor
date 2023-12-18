import pickle
import argparse
from ..features.build_features import engineer_features, process_pharokka_output
import pandas as pd
import os
import numpy as np
import subprocess
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import matplotlib
matplotlib.use('Agg')

parser = argparse.ArgumentParser(description="Make a gram staining prediction for a set of phages sequences.")
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
parser.add_argument("-d", "--data", type=str, required=True,
                    help="Path to the data to be predicted. Accepted inputs are FASTA files for Pharokka inputs or a directory path to an already processed Pharokka output.\
                        If the input is a FASTA file, the Pharokka output will be stored locally and predictions will be based on it.")
optional.add_argument("-m","--model_path", dest="model", action='store', default="models/random_forest.pkl", 
                      help="Path to the trained model file (default: models/model.pkl).")
optional.add_argument("-p","--pharokka_directory", dest="pharokka", action='store', 
                          default="/mnt/c/Users/Alvaro/Desktop/pharokka_output/", help="Path where pharokka output will be stored (default: /mnt/c/Users/Alvaro/Desktop/pharokka_output/).")
optional.add_argument("-b","--database_directory", dest="database", action='store', 
                          default="/mnt/c/Users/Alvaro/Desktop/pharokka_database/", help="Path to the Pharokka database (default./mnt/c/Users/Alvaro/Desktop/pharokka_database/).")
optional.add_argument("-f","--output_features", dest="output_features", action='store', default='reports/features_10.csv', help="Output path (default: 'reports/features10.csv').")
optional.add_argument("-o","--output_predictions", dest="output_predictions", action='store', default='reports/predictions_user.csv', help="Output path (default: 'reports/predictions_user.csv').")
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

    
    
if os.path.isdir(args.data):
    print(f"Processing Pharokka phages in directory: {args.data}")
    process_pharokka_output(args.data, args.output_features)
    
    # Print completion message
    print(f"Feature engineering completed for Prediction data and data stored to CSV in {args.output_features}.")


elif os.path.isfile(args.data):
    # The input is a file; check if it's a FASTA file
    file_ext = os.path.splitext(args.data)[1][1:].lower()
    if file_ext in fasta_extensions:
        print("The processing of FASTA files with Pharokka is currently unavailable. Please, submit a folder with Pharokka outputs.")
        # run_pharokka(args_data, args_pharokka, args_database) 

        # # The input is a directory; process all relevant files within
        # print(f"Processing Pharokka output in directory: {args_pharokka}")
        # model_data = engineer_features(args_pharokka)
    elif file_ext == 'csv':
            # New code to handle CSV files
            try:
                model_data = pd.read_csv(args.data)
                features = [
                    'genome_length_inphared', 'gc_%_inphared', 'cds_number_inphared', 
                    'positive_strand_%_inphared', 'negative_strand_%_inphared', 
                    'coding_capacity_inphared', 'tRNAs_inphared', 'cds_coding_density',
                    'jumbophage_inphared', 'topology_linear', 'topology_circular', 
                    'transl_table', 'CARD_AMR_Genes', 'CRISPRs', 
                    'VFDB_Virulence_Factors', 'connector', 'head_packaging', 
                    'host_takeover', 'integration and excision', 'lysis', 
                    'nucleotide_metabolism', 'other', 'tail', 'tmRNAs', 
                    'transcription', 'unkown_function', 'frame_positive', 
                    'frame_negative', 'molecule_inphared_type_DNA', 
                    'molecule_inphared_type_ss-DNA', 'molecule_inphared_type_RNA', 
                    'molecule_inphared_type_ss-RNA'
                ]

                if all(column in model_data.columns for column in features):
                    model_data = model_data.set_index(model_data.columns[0])

                    model_data = model_data[features]

                    # Making the prediction -----

                    # Load the model using the provided path
                    model = pickle.load(open(args.model, "rb"))

                    # Load the trained scaler
                    with open("models/scaler.pkl", 'rb') as scaler_file:
                        scaler = pickle.load(scaler_file)


                    model_data_scaled = scaler.transform(model_data)

                    # Predictions
                    new_data_pred = model.predict(model_data_scaled)

                    # Get the probabilities for the predicted class for each instance
                    probas = model.predict_proba(model_data_scaled)
                    predicted_indices = np.argmax(probas, axis=1)  # Get index of max proba for each sample
                    new_data_pred_proba = [probas[i][predicted_indices[i]] for i in range(len(predicted_indices))]

                    # Prepare the output DataFrame
                    output_df = pd.DataFrame({
                        'id': model_data.index,
                        'prediction': new_data_pred,
                        'prediction_probability': new_data_pred_proba
                    })
                    output_df = output_df.sort_values(by='prediction_probability', ascending=False)

                    # Uncomment the following line if you want to save the output to a CSV
                    output_df.to_csv(args.output_predictions, index=False)

                    print()
                    print(f"The output has been saved in {args.output_predictions}. The first 5 entries are: ")
                    print(output_df.head())

                else:
                    print("The CSV file does not contain the expected features.")
                    exit(1)

            except Exception as e:
                print(f"Error reading CSV file: {e}")
                exit(1)

    else:
        print(f"Unsupported file format: {file_ext}. Supported FASTA formats are folders, csv and: {', '.join(fasta_extensions)}")
        exit(1)
else:
    print(f"The provided data path does not exist: {args.data}")
    exit(1)

output_df = pd.read_csv(args.output_predictions)

print("Saving the results as figures")
# KDE Plot
sns.set_style("whitegrid")
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=FutureWarning)

    plt.figure(figsize=(10, 6))
    sns.kdeplot(output_df[output_df['prediction'] == True]['prediction_probability'], fill=False, label='Positive', color="blue", clip=(0, 1))
    sns.kdeplot(output_df[output_df['prediction'] == False]['prediction_probability'], fill=False, label='Negative', color="red", clip=(0, 1))
    sns.kdeplot(output_df['prediction_probability'], fill=False, label='Combined', color='grey', clip=(0, 1))
    plt.xlabel('Prediction Probability')
    plt.ylabel('Density')
    plt.title('Comparative Density of Positive and Negative Predictions')
    plt.legend()

    # Save KDE plot
    plt.tight_layout()
    plt.savefig("reports/" + 'kde_plot.png')
    plt.close()

# Histogram Plot
sns.set_style("whitegrid")
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=FutureWarning)

    plt.figure(figsize=(10, 6))
    sns.histplot(output_df[output_df['prediction'] == True]['prediction_probability'], kde=False, label='Positive', color="blue", bins=20)
    sns.histplot(output_df[output_df['prediction'] == False]['prediction_probability'], kde=False, label='Negative', color="red", bins=20)
    sns.histplot(output_df['prediction_probability'], kde=False, label='Combined', color='grey', bins=20, alpha=0.3)
    plt.xlabel('Prediction Probability')
    plt.ylabel('Frequency')
    plt.title('Distribution of Prediction Probabilities')
    plt.legend()

    # Save histogram plot
    plt.tight_layout()
    plt.savefig("reports/" + 'histogram_plot.png')
    plt.close()

# python -m src.models.predict_model -d data/interim/genbank_engineering/50_sequences.gb