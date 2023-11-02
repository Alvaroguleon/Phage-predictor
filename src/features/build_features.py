import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import UndefinedSequenceError
import argparse

def engineer_features(genbank_file):
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
        'id': ids,
        'genome_length': genome_lengths,
        'gc_%': gc_contents,
        'sequence': sequences,
        'reverse_complement': reverse_complements,
        'cds_number': cds_numbers,
        'positive_strand_%': positive_strands,
        'negative_strand_%': negative_strands,
        'coding_capacity': coding_capacities,
        'molecule_type': molecule_types,
        'topology': topologies,
        'trna_count': trna_counts
    })


    df['id'] = df['id'].str[:-2]

    # Check for unexpected molecule types
    expected_molecule_types = ['ss-DNA', 'DNA', 'RNA', 'ss-RNA']

    # Check and correct 'cRNA' entries
    cRNA_entries = df[df['molecule_type'] == 'cRNA']
    if not cRNA_entries.empty:
        for entry_id in cRNA_entries['id']:
            print(f"Info: Entry with id '{entry_id}' has molecule type 'cRNA'. Changing it to 'RNA'.")
        df.loc[df['molecule_type'] == 'cRNA', 'molecule_type'] = 'RNA'

    # Check and correct 'cDNA' entries
    cDNA_entries = df[df['molecule_type'] == 'cDNA']
    if not cDNA_entries.empty:
        for entry_id in cDNA_entries['id']:
            print(f"Info: Entry with id '{entry_id}' has molecule type 'cDNA'. Changing it to 'DNA'.")
        df.loc[df['molecule_type'] == 'cDNA', 'molecule_type'] = 'DNA'

    unexpected_types = set(df['molecule_type']) - set(expected_molecule_types)

    if unexpected_types:
        for utype in unexpected_types:
            # Get ids of entries with the unexpected molecule type
            ids_to_exclude = df[df['molecule_type'] == utype]['id'].tolist()
            for entry_id in ids_to_exclude:
                print(f"Warning: Entry with id '{entry_id}' has unrecognized molecule type '{utype}'. It will not be considered.")
            df = df[df['molecule_type'] != utype]
            
    df = pd.get_dummies(df, columns=['molecule_type'])

    expected_columns = ['jumbophage', 'molecule_type_ss-DNA', 'molecule_type_DNA', 'molecule_type_RNA', 'molecule_type_ss-RNA']
    for col in expected_columns:
        if col not in df.columns:
            df[col] = 0  # Filling with zeros
        df[col] = df[col].astype(bool)  # Convert to boolean

    df['jumbophage'] = df['genome_length'].apply(lambda x: x >= 200000)
    df['jumbophage'] = df['jumbophage'].astype(int)  # Convert True/False to 1/0
    df = pd.get_dummies(df, columns=['topology'])
    return df

def staining_feature(staining_df, features_df):
    stain = pd.read_csv(staining_df, index_col = 0)
    stain = stain[['Accession', 'staining']]
    stain = stain.rename(columns={'Accession': 'id'})

    features_df = pd.merge(features_df, stain, on='id', how='left')
    features_df = features_df[['id', 'staining', 'genome_length', 'jumbophage', 'gc_%', 'trna_count','cds_number', 'coding_capacity',
                               'positive_strand_%','negative_strand_%', 'molecule_type_ss-DNA', 'molecule_type_DNA', 'molecule_type_RNA', 'molecule_type_ss-RNA',
                               'topology_circular', 'topology_linear']]

    return features_df

def main():
    parser = argparse.ArgumentParser(description="Load a trained model from a given path.")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    parser.add_argument("-g", "--genbank", type=str, required=True,
                        help="Path to genbank file of the entry.")
    optional.add_argument("-s","--staining", dest="staining", action='store', default='data/interim/gram_staining/staining_assignation.csv', help="Path to csv file with a customised gram staining class assignation for every entry.")
    optional.add_argument("-o","--output_direcotry", dest="output", action='store', default='data/processed/model_data2.csv', help="Output path (default: data/processed/model_data.csv).")

    args = parser.parse_args()

    # The genbank file input should be the output of pharokka 

    print("Reading the genbank file...")
    features_df = engineer_features(args.genbank)

    features_df = staining_feature(args.staining, features_df)

    features_df.to_csv(args.output, index = False)
    print(f"Processed data has been stored in {args.output}")

if __name__ == "__main__":
    main()



