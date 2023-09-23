# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import pandas as pd

@click.command()
@click.argument('input_data_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_data_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')

    print("Processing phage data...")
    df = pd.read_csv(input_data_filepath, sep="\t")
    df.rename(columns={'Isolation Host (beware inconsistent and nonsense values)': 'Isolation Host'}, inplace=True) 

    # Save the phages with no known Host in a separate file
    discarded = df[(df['Host'] == 'Unspecified') | (df['Isolation Host'] == 'Unspecified')]
    discarded = discarded[['Accession','Host', 'Isolation Host']]

    # Filter out phages with unknown hosts from the phage dataframe
    df = df[df['Host'] != "Unspecified"] 
    df = df[df['Isolation Host'] != "Unspecified"]

    # Fixing misspelled genera, taxonomy, etc
    df['Isolation Host'] = df['Isolation Host'].str.split().str[:2].str.join(' ')
    df['Host'] = df['Host'].replace('Enteroccous', 'Enterococcus')
    df.loc[df['Isolation Host'] == 'Salmonella enterica', 'Host'] = 'Salmonella'
    df.loc[df['Isolation Host'] == 'Salmonella typhimurium', 'Host'] = 'Salmonella'
    df.loc[df['Isolation Host'] == 'Escherichia coli', 'Host'] = 'Escherichia'
    df.loc[df['Isolation Host'] == 'E. coli', 'Host'] = 'Escherichia'
    df.loc[df['Isolation Host'] == 'Escherichia coli,', 'Host'] = 'Escherichia'
    df.loc[df['Isolation Host'] == 'Escherichia coli;', 'Host'] = 'Escherichia'
    df.loc[df['Isolation Host'] == 'Shigella flexneri', 'Host'] = 'Shigella'

    print("Processing gram staining data...")
    
    # Load dataset to assign gram staining to each phage based on known host
    gram_information = pd.read_csv("../data/interim/gram_staining_dict.csv")

    # Excluding the elements which are mapped to both stainings - This should be changed to be specific to the data in use
    exclude_genus = ['Clostridium', 'Neobacillus', 'Alteribacter', 'Desulfotomaculum', 'Caloramator', 'Desulforamulus', 'Heyndrickxia', 'Peptoclostridium', 'Thermoanaerobacter', 'Thermoanaerobacterium', 'Aureimonas', 'Actinomadura', 'Alkalibacterium', 'Deinococcus', 'Tepidibacillus', 'Sphingomonas', 'Lysinibacillus', 'Ruminiclostridium', 'Caldicellulosiruptor', 'Pseudomonas', 'Streptococcus', 'Microlunatus', 'Streptomyces', 'Butyricimonas', 'Halalkalibacter', 'Chelativorans', 'Natrinema', 'Ureibacillus', 'Clostridioides', 'Desulfosporosinus', 'Lacibacter', 'Nocardioides', 'Siminovitchia', 'Belliella', 'Tistlia', 'Actinoplanes', 'Paenibacillus', 'Vallitalea', 'Actinotalea', 'Cohnella', 'Rhizobium', 'Anaerotignum', 'Cellulomonas', 'Flavobacterium', 'Bacillus', 'Nesterenkonia']
    exclude_species = ['Ureibacillus massiliensis', 'Tistlia consotensis', 'Clostridioides difficile', 'Vallitalea guaymasensis', 'Belliella pelovolcani', 'Actinotalea ferrariae']

    # Save the discarded dataframe
    excluded_rows = df[df['Host'].isin(exclude_genus) | df['Isolation Host'].isin(exclude_species)]

    # Append the excluded rows to the discarded DataFrame
    discarded = pd.concat([discarded, excluded_rows], ignore_index=True)
    discarded.to_csv(output_filepath + "/discarded_phages.csv")

    df = df[~df['Host'].isin(exclude_genus)]
    df = df[~df['Isolation Host'].isin(exclude_species)]

    # Map staining from Isolation Host and Host
    species_to_stain = dict(zip(all_gram['species'], all_gram['Gram stain']))
    df['staining'] = df['Isolation Host'].map(species_to_stain)

    genus_to_stain = dict(zip(all_gram['Genus'], all_gram['Gram stain']))
    df['staining'] = df['Host'].map(genus_to_stain)

    # Obtaining remaining phages without assignation
    remaining = df[(df["staining"] != "positive") & (df["staining"] != "negative")]

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
