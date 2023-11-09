import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import UndefinedSequenceError
import argparse
import os
import subprocess

def engineer_features(folder):
    # Lists to hold data
    ids = []
    genome_lengths = []
    gc_contents = []
    sequences = []
    reverse_complements = []
    cds_numbers = []
    positive_strands = []
    negative_strands = []
    coding_capacities = []
    molecule_types = []
    topologies = []
    trna_counts = []

    genbank_file = f'{folder}/pharokka.gbk'

    # Read the GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        try:
            # Attempt to access the sequence, which may raise UndefinedSequenceError
            sequence = str(record.seq)
            # print(record.id)
        except UndefinedSequenceError:
            # print(f"Skipping record {record.id} as sequence is undefined.")
            continue  # Skip this record

        # Calculate genome length and GC content
        total_length = len(sequence)
        gc_content = round((sequence.count('C') + sequence.count('G')) / total_length * 100, 3)

        # Initialize counters
        plus = 0
        minus = 0
        coding_count = 0
        trna_count = 0
        seen = set()  # Store seen barcodes

        for feature in record.features:
            start = feature.location.start
            end = feature.location.end
            length = len(FeatureLocation(start, end))
            barcode = f"{start}_{end}_{length}"

            if feature.type != 'source' and barcode not in seen:
                coding_count += length
                seen.add(barcode)

            if feature.type == 'CDS':
                if feature.location.strand == 1:
                    plus += 1
                elif feature.location.strand == -1:
                    minus += 1
            elif feature.type == 'tRNA':
                trna_count += 1

        
        # Calculate total number of CDS
        total_CDS = plus + minus

        # Calculate strand usage as a percentage
        per_plus = round((plus / total_CDS) * 100, 2) if total_CDS != 0 else 0
        per_minus = round((minus / total_CDS) * 100, 2) if total_CDS != 0 else 0

        # Calculate coding capacity as a percentage
        coding_capacity = (coding_count / total_length) * 100

        # Extract molecule_type and topology
        molecule_type = record.annotations.get('molecule_type', 'N/A')
        topology = record.annotations.get('topology', 'N/A')

        # Append data to lists
        ids.append(record.id)
        genome_lengths.append(total_length)
        gc_contents.append(gc_content)
        sequences.append(sequence)
        reverse_complements.append(str(sequence[::-1]))
        cds_numbers.append(total_CDS)
        positive_strands.append(per_plus)
        negative_strands.append(per_minus)
        coding_capacities.append(coding_capacity)
        molecule_types.append(molecule_type)
        topologies.append(topology)
        trna_counts.append(trna_count)
    
    print("Processing the entries...")
    # Convert lists to pandas DataFrame
    df = pd.DataFrame({
        'id_inphared': ids,
        'genome_length_inphared': genome_lengths,
        'gc_%_inphared': gc_contents,
        'sequence_inphared': sequences,
        'reverse_complement_inphared': reverse_complements,
        'cds_number_inphared': cds_numbers,
        'positive_strand_%_inphared': positive_strands,
        'negative_strand_%_inphared': negative_strands,
        'coding_capacity_inphared': coding_capacities,
        'molecule_type_inphared': molecule_types,
        'topology_inphared': topologies
    })


    df['id_inphared'] = df['id_inphared'].str[:-2]

    # Check for unexpected molecule types
    expected_molecule_types = ['ss-DNA', 'DNA', 'RNA', 'ss-RNA']

    # Check and correct 'cRNA' entries
    cRNA_entries = df[df['molecule_type_inphared'] == 'cRNA']
    if not cRNA_entries.empty:
        for entry_id in cRNA_entries['id_inphared']:
            print(f"Info: Entry with id '{entry_id}' has molecule type 'cRNA'. Changing it to 'RNA'.")
        df.loc[df['molecule_type_inphared'] == 'cRNA', 'molecule_type_inphared'] = 'RNA'

    # Check and correct 'cDNA' entries
    cDNA_entries = df[df['molecule_type_inphared'] == 'cDNA']
    if not cDNA_entries.empty:
        for entry_id in cDNA_entries['id_inphared']:
            print(f"Info: Entry with id '{entry_id}' has molecule type 'cDNA'. Changing it to 'DNA'.")
        df.loc[df['molecule_type_inphared'] == 'cDNA', 'molecule_type_inphared'] = 'DNA'

    unexpected_types = set(df['molecule_type_inphared']) - set(expected_molecule_types)

    if unexpected_types:
        for utype in unexpected_types:
            # Get ids of entries with the unexpected molecule type
            ids_to_exclude = df[df['molecule_type_inphared'] == utype]['id'].tolist()
            for entry_id in ids_to_exclude:
                print(f"Warning: Entry with id '{entry_id}' has unrecognized molecule type '{utype}'. It will not be considered.")
            df = df[df['molecule_type'] != utype]
            
    df = pd.get_dummies(df, columns=['molecule_type_inphared'])
    df = pd.get_dummies(df, columns=['topology_inphared'])
    expected_columns = ['jumbophage_inphared', 'topology_linear_inphared','topology_circular_inphared',
    'molecule_inphared_type_ss-DNA', 'molecule_inphared_type_DNA', 'molecule_inphared_type_RNA', 'molecule_inphared_type_ss-RNA']
    
    for col in expected_columns:
        if col not in df.columns:
            df[col] = 0  # Filling with zeros
        df[col] = df[col].astype(bool)  # Convert to boolean

    df['jumbophage_inphared'] = df['genome_length_inphared'].apply(lambda x: x >= 200000)
    df['jumbophage_inphared'] = df['jumbophage_inphared'].astype(int)  # Convert True/False to 1/0
    
    df['id_inphared'] = df['id_inphared'].str.slice(0,8)


    df_length_gc_cds_density = pd.read_csv(f"{folder}pharokka_length_gc_cds_density.tsv", sep="\t")
    df_length_gc_cds_density['contig'] = df_length_gc_cds_density['contig'].str.slice(0,8)
    df_length_gc_cds_density

    df = pd.merge(df, df_length_gc_cds_density, left_on='id_inphared', right_on='contig', how='outer')
    del(df_length_gc_cds_density)

    df_cds = pd.read_csv(f"{folder}pharokka_cds_functions.tsv", sep="\t")
    
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

    # Specify dtype to avoid importing issues 
    dtype_spec = {
        'vfdb_species': 'str',
        'CARD_eVal': 'float',  # Float because it's numeric but can contain NaN
        'CARD_species': 'str',
        'ARO_Accession': 'str',
        'CARD_short_name': 'str',
        'Protein_Accession': 'str',
        'DNA_Accession': 'str',
        'AMR_Gene_Family': 'str',
        'Drug_Class': 'str'
    }

    df_frame = pd.read_csv(f"{folder}pharokka_cds_final_merged_output.tsv", sep="\t", dtype=dtype_spec, low_memory=False)

    frame_counts = df_frame.groupby('contig')['frame'].value_counts().unstack(fill_value=0)
    frame_counts.columns = ['frame_negative', 'frame_positive']
    frame_counts = frame_counts.rename_axis(None, axis=1).reset_index()
    frame_counts['contig'] = frame_counts['contig'].str.slice(0,8)
    df = pd.merge(df, frame_counts, left_on='id_inphared', right_on='contig', how='outer')

    df = df.rename(columns={'id_inphared': 'id'})

    # # List of columns to be ignored in the duplicate check
    # ignore_columns = {
    #     'CARD_AMR_Genes', 'CDS', 'CRISPRs', 'VFDB_Virulence_Factors', 'connector',
    #     'head_packaging', 'host_takeover', 'integration and excision', 'lysis',
    #     'nucleotide_metabolism', 'other', 'tRNAs', 'tail', 'tmRNAs',
    #     'transcription', 'unkown_function'
    # }

    # # Function to check for repeated values excluding 0 and 1 and certain columns
    # def find_repeated_values(row):
    #     # Use defaultdict to associate multiple columns with the same value
    #     repeated_columns = defaultdict(list)
    #     # Use items() instead of iteritems()
    #     for col, value in row.items():
    #         if col in ignore_columns:  # Skip the iteration for specified columns
    #             continue
    #         if (value in [0, 1]) or pd.isnull(value):  # Skip values of 0, 1, or NaN
    #             continue
    #         if row[row == value].count() > 1:  # If value count in row is more than one, it's a duplicate
    #             repeated_columns[value].append(col)
    #     return repeated_columns

    # # Ensure row is a Series before passing to find_repeated_values
    # for index, row in df.iterrows():
    #     repeated_values_dict = find_repeated_values(row)
    #     if repeated_values_dict:  # If there's any repeated value
    #         warning_message = f"Entry {row['id']} shows repeated values in "
    #         # Create a combined message for duplicate columns sharing the same value
    #         message_parts = []
    #         for val, cols in repeated_values_dict.items():
    #             if len(cols) > 1:  # Only add to message if there are indeed duplicates
    #                 combined_cols = ' and '.join(cols)  # Combine column names
    #                 message_parts.append(f"columns {combined_cols} (value {val})")
    #         warning_message += ', '.join(message_parts)
    #         print(warning_message)
    return df

def staining_feature(staining_df, features_df):
    stain = pd.read_csv(staining_df, index_col = 0)
    stain = stain[['Accession', 'staining']]
    stain = stain.rename(columns={'Accession': 'id'})

    features_df = pd.merge(features_df, stain, on='id', how='left')
    features_df = features_df[['id', 'staining', 'genome_length_inphared', 'gc_%_inphared',
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
        'frame_positive', 'sequence_inphared', 'reverse_complement_inphared']]


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
    optional.add_argument("-p","--pharokka_directory", dest="pharokka", action='store', 
                          default="/mnt/c/Users/Alvaro/Desktop/pharokka_output/", help="Path where pharokka output will be stored (default: /mnt/c/Users/Alvaro/Desktop/pharokka_output/).")
    optional.add_argument("-b","--database_directory", dest="database", action='store', 
                          default="/mnt/c/Users/Alvaro/Desktop/pharokka_database/", help="Path to the Pharokka database (default./mnt/c/Users/Alvaro/Desktop/pharokka_database/).")
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



