import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import UndefinedSequenceError
import argparse
import os
import subprocess

def process_pharokka_output(args_data, output_csv):
    # Check if output CSV file exists, if not create it
    if not os.path.exists(output_csv):
        with open(output_csv, 'w') as f:
            pass

    # Iterate over each folder in the main folder
    for folder in os.listdir(args_data):
        folder_path = os.path.join(args_data, folder)

        # Check if the path is a directory
        if os.path.isdir(folder_path):
            # Check if "pharokka.gbk" file exists in the folder
            if "pharokka.gbk" in os.listdir(folder_path):
                print(f"Processing {folder_path}...")
                # Pass the full path of "pharokka.gbk" to your function and get the dataframe
                df = engineer_features(folder_path)

                # Append the dataframe to the CSV file
                with open(output_csv, 'a') as f:
                    df.to_csv(f, header=f.tell()==0, index=False)
            else:
                print(f"Skipped {folder_path}: 'pharokka.gbk' file not found or not a directory")

def engineer_features(folder):
    # Read data from the TSV file
    df = pd.read_csv(f'{folder}/pharokka_top_hits_mash_inphared.tsv', sep='\t')

    # Select only the columns you need
    selected_columns = [
        'Accession', 'contig','Genome_Length_(bp)', 'molGC_(%)', 'Molecule',
        'Number_CDS', 'Positive_Strand_(%)', 'Negative_Strand_(%)',
        'Coding_Capacity_(%)', 'tRNAs', 'Host', 
        'Isolation_Host_(beware_inconsistent_and_nonsense_values)'
    ]
    df = df[selected_columns]

    # Rename columns to match your previous DataFrame structure
    df = df.rename(columns={
        'Accession': 'id_inphared',
        'contig': 'contig_inphared',
        'Genome_Length_(bp)': 'genome_length_inphared',
        'molGC_(%)': 'gc_%_inphared',
        'Molecule': 'molecule_inphared_type',
        'Number_CDS': 'cds_number_inphared',
        'Positive_Strand_(%)': 'positive_strand_%_inphared',
        'Negative_Strand_(%)': 'negative_strand_%_inphared',
        'Coding_Capacity_(%)': 'coding_capacity_inphared',
        'tRNAs': 'tRNAs_inphared',
        'Host': 'host_inphared',
        'Isolation_Host_(beware_inconsistent_and_nonsense_values)': 'isolation_host_inphared'
    })

    # Initialize lists for storing data
    ids = []
    topologies = []
    sequences = []

    # Read the GenBank file for topology and sequences
    genbank_file = f'{folder}/pharokka.gbk'
    for record in SeqIO.parse(genbank_file, "genbank"):
        ids.append(record.id)    
        topologies.append(record.annotations.get('topology', 'N/A'))
        sequences.append(str(record.seq))


    # These lines use the "contig" column of inphared tsv when the "Accession" is wrong (does not match the one in Genbank file?)
    genbank_id = ids[0] if ids else None

    # Update Accession value based on the logic described
    for index, row in df.iterrows():
        if row['id_inphared'] != genbank_id and row['contig_inphared'] == genbank_id:
            print(f"Updating Accession for folder {folder}: from {row['id_inphared']} to {genbank_id}")
            df.at[index, 'id_inphared'] = genbank_id



    # Create DataFrames from the topologies and sequences
    topology_df = pd.DataFrame({'id_inphared': ids, 'topology': topologies})
    sequence_df = pd.DataFrame({'id_inphared': ids, 'sequence': sequences})

    # Merge the TSV data with the GenBank topology and sequence data
    df = pd.merge(df, topology_df, on='id_inphared', how='left')
    df = pd.merge(df, sequence_df, on='id_inphared', how='left')
    del(topology_df)
    del(sequence_df)

    # Check for unexpected molecule types
    expected_molecule_types = ['ss-DNA', 'DNA', 'RNA', 'ss-RNA']

    # Check and correct 'cRNA' entries
    cRNA_entries = df[df['molecule_inphared_type'] == 'cRNA']
    if not cRNA_entries.empty:
        for entry_id in cRNA_entries['id_inphared']:
            print(f"Info: Entry with id '{entry_id}' has molecule type 'cRNA'. Changing it to 'RNA'.")
        df.loc[df['molecule_inphared_type'] == 'cRNA', 'molecule_inphared_type'] = 'RNA'

    # Check and correct 'cDNA' entries
    cDNA_entries = df[df['molecule_inphared_type'] == 'cDNA']
    if not cDNA_entries.empty:
        for entry_id in cDNA_entries['id_inphared']:
            print(f"Info: Entry with id '{entry_id}' has molecule type 'cDNA'. Changing it to 'DNA'.")
        df.loc[df['molecule_inphared_type'] == 'cDNA', 'molecule_inphared_type'] = 'DNA'

    unexpected_types = set(df['molecule_inphared_type']) - set(expected_molecule_types)

    if unexpected_types:
        for utype in unexpected_types:
            # Get ids of entries with the unexpected molecule type
            ids_to_exclude = df[df['molecule_inphared_type'] == utype]['id_inphared'].tolist()
            for entry_id in ids_to_exclude:
                print(f"Warning: Entry with id '{entry_id}' has unrecognized molecule type '{utype}'. It will not be considered.")
            df = df[df['molecule_inphared_type'] != utype]
            
    df = pd.get_dummies(df, columns=['molecule_inphared_type'])
    df = pd.get_dummies(df, columns=['topology'])
    expected_columns = ['jumbophage_inphared', 'topology_linear','topology_circular',
    'molecule_inphared_type_ss-DNA', 'molecule_inphared_type_DNA', 'molecule_inphared_type_RNA', 'molecule_inphared_type_ss-RNA']
    
    for col in expected_columns:
        if col not in df.columns:
            df[col] = 0  # Filling with zeros
        df[col] = df[col].astype(bool)  # Convert to boolean

    df['jumbophage_inphared'] = df['genome_length_inphared'].apply(lambda x: x >= 200000)
    df['jumbophage_inphared'] = df['jumbophage_inphared'].astype(int)  # Convert True/False to 1/0
    

    df_length_gc_cds_density = pd.read_csv(f"{folder}/pharokka_length_gc_cds_density.tsv", sep="\t")
    df_length_gc_cds_density['contig'] = df_length_gc_cds_density['contig'].str.slice(0,8)
    df_length_gc_cds_density

    df = pd.merge(df, df_length_gc_cds_density, left_on='id_inphared', right_on='contig', how='outer')
    del(df_length_gc_cds_density)

    df_cds = pd.read_csv(f"{folder}/pharokka_cds_functions.tsv", sep="\t")
    
    # Define a dictionary of replacements
    replacements = {
        "DNA, RNA and nucleotide metabolism":"nucleotide_metabolism",
        "head and packaging": "head_packaging",
        "moron, auxiliary metabolic gene and host takeover": "host_takeover",
        "transcription regulation": "transcription",
        "unknown function": "unkown_function"
    }

    # Replace the values using the dictionary
    df_cds['Description'] = df_cds['Description'].replace(replacements)

    df_cds = df_cds.pivot(index='contig', columns='Description', values='Count').reset_index()
    df_cds['contig'] = df_cds['contig'].str.slice(0,8)

    df = pd.merge(df, df_cds, left_on='id_inphared', right_on='contig', how='outer')
    del(df_cds)


    df_frame = pd.read_csv(f"{folder}/pharokka_cds_final_merged_output.tsv", sep="\t",  low_memory=False)

    frame_counts = df_frame.groupby('contig')['frame'].value_counts().unstack(fill_value=0)

    # Ensure both '+' and '-' columns are present
    frame_counts['+'] = frame_counts.get('+', 0)
    frame_counts['-'] = frame_counts.get('-', 0)

    # Rename columns explicitly
    frame_counts = frame_counts.rename(columns={'+': 'frame_positive', '-': 'frame_negative'})

    # Reset index if needed
    frame_counts = frame_counts.reset_index()
    frame_counts = frame_counts.rename_axis(None, axis=1).reset_index()
    frame_counts['contig'] = frame_counts['contig'].str.slice(0,8)
    df = pd.merge(df, frame_counts, left_on='id_inphared', right_on='contig', how='outer')

    df = df.rename(columns={'id_inphared': 'id'})
    df['dummy_index'] = 0

    # Group by the dummy index and aggregate using 'first'
    df = df.groupby('dummy_index').first()

    # Reset the index
    df = df.reset_index(drop=True)

    columns = ['id','host_inphared','isolation_host_inphared','genome_length_inphared', 'gc_%_inphared', 'cds_number_inphared',
       'positive_strand_%_inphared', 'negative_strand_%_inphared',
       'coding_capacity_inphared', 'tRNAs_inphared', 'host_inphared',
       'isolation_host_inphared', 'molecule_inphared_type_DNA',
       'topology_linear', 'jumbophage_inphared', 'topology_circular',
       'molecule_inphared_type_ss-DNA', 'molecule_inphared_type_RNA',
       'molecule_inphared_type_ss-RNA',  'length', 'gc_perc',
       'transl_table', 'cds_coding_density',  'CARD_AMR_Genes',
       'CDS', 'CRISPRs', 'VFDB_Virulence_Factors', 'connector',
       'head_packaging', 'host_takeover', 'integration and excision', 'lysis',
       'nucleotide_metabolism', 'other', 'tRNAs', 'tail', 'tmRNAs',
        #'contig_x', 'contig', 'contig_y',
       'transcription', 'unkown_function', 'frame_positive',
       'frame_negative', 'sequence']
    
    columns = ['id','host_inphared','isolation_host_inphared','genome_length_inphared', 'gc_%_inphared',
       'cds_number_inphared', 'positive_strand_%_inphared',
       'negative_strand_%_inphared', 'coding_capacity_inphared',
       'tRNAs_inphared', 'cds_coding_density','jumbophage_inphared', 'topology_linear', 
       'topology_circular', 'transl_table',  'CARD_AMR_Genes',
       'CRISPRs', 'VFDB_Virulence_Factors', 'connector',
       'head_packaging', 'host_takeover', 'integration and excision', 'lysis',
       'nucleotide_metabolism', 'other', 'tail', 'tmRNAs',
       'transcription', 'unkown_function', 'frame_positive', 'frame_negative', 
       'molecule_inphared_type_DNA',  'molecule_inphared_type_ss-DNA', 'sequence']

    df = df[columns]
    return df

def staining_feature(staining_df, features_df):
    stain = pd.read_csv(staining_df, index_col = 0)
    stain = stain[['Accession', 'staining']]
    stain = stain.rename(columns={'Accession': 'id'})

    features_df = pd.merge(features_df, stain, on='id', how='left')

    return features_df
    
def main():
    parser = argparse.ArgumentParser(description="Process raw data (filtered data with only valid phage) into features to train a model.")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    parser.add_argument("-d", "--data", type=str, required=True,
                        help="Path to the data. Accepted inputs are FASTA files for Pharokka inputs or a directory path to an already processed Pharokka output.\
                            If the input is a FASTA file, the Pharokka output will be stored locally.")
    optional.add_argument("-s","--staining", dest="staining", action='store', 
                          default='data/interim/gram_staining/staining_assignation.csv', 
                          help="Path to csv file with a customised gram staining class assignation for every entry.")
    optional.add_argument("-o","--output_directory", dest="output", action='store', 
                          default='data/processed/model_data_pharokka.csv', help="Output path (default: data/processed/model_data_pharokka.csv).")

    args = parser.parse_args()

    '''
    Possible inputs:
    - Pharokka output
    - Fasta file on which to build pharokka output

    This output will be the input for the engineer features and staining feature, and its
    output is stored in args.output. This output is used by train_model for training
    '''
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


    model_data = staining_feature(args.staining, model_data)

    model_data.to_csv(args.output, index = False)
    
    print(f"Processed data has been stored in {args.output}")

if __name__ == "__main__":
    main()



