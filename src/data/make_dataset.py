# -*- coding: utf-8 -*-
import argparse
import logging
from logging.handlers import RotatingFileHandler
from datetime import datetime
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import pandas as pd
import os 
import subprocess

# from ..features.build_features import engineer_features

fasta_extensions = ["fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"]


def main(args_data, args_output, args_pharokka_dir, args_pharokka_database):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('-> Start of Make Dataset <-')

    # Take input FASTA and run Pharokka
    if os.path.isfile(args_data):
        # The input is a file; check if it's a FASTA file
        file_ext = os.path.splitext(args_data)[1][1:].lower()
        if file_ext in fasta_extensions:
            logger.error("The processing of FASTA files with Pharokka is currently unavailable. Please, submit a folder with Pharokka outputs.")
            # run_pharokka(args_data, args_pharokka_dir, args_pharokka_database) 

            # # The input is a directory; process all relevant files within
            # logger.info(f"Processing Pharokka output in directory: {args_pharokka_dir}")
            # model_data = engineer_features(args_pharokka_dir)

        else:
            logger.error(f"Unsupported file format: {file_ext}. Supported FASTA formats are: {', '.join(fasta_extensions)}")
            exit(1)

    # Take Pharokka output folder as input
    if os.path.isdir(args_data):
        # The input is a directory; process all relevant files within
        logger.info(f"Processing Pharokka output in directory: {args_data}")

        # Define the command to run the external script
        command = ["python","-m","src.features.build_features", "-d", args_data]

        # Run the script
        subprocess.run(command, check=True)
        logger.info(f"Pharokka output in directory {args_data} has been processed.")


    else:
        logger.error(f"The provided data path does not exist: {args_data}")
        exit(1)


if __name__ == '__main__':
    # Set up basic configuration for logging
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # Get current date
    current_date = datetime.now().strftime('%Y-%m-%d')

    # Create a logger
    logger = logging.getLogger(__name__)

    # Create a file handler that logs even debug messages with the current date in filename
    # Directory where you want to save the log files
    log_directory = '/mnt/c/Users/Alvaro/Desktop/projects/phage/reports/logs/'

    # Create the directory if it doesn't exist
    Path(log_directory).mkdir(parents=True, exist_ok=True)

    # Filename for the log file
    log_filename = f'{log_directory}log_{current_date}.log'
    file_handler = RotatingFileHandler(log_filename, maxBytes=10000, backupCount=1)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(log_fmt))

    # Add the file handler to the logger
    logger.addHandler(file_handler)

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    # load_dotenv(find_dotenv())

    # Create the parser
    parser = argparse.ArgumentParser(description="Data processing script.")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional.add_argument("-d", "--data", dest="data", action='store', default="/mnt/c/Users/Alvaro/Desktop/projects/phage/data/interim/pharokka/10_folders/", 
                    help="Path to the raw data. Accepted input is a FASTA file with phage sequences (and temporarily, also a Pharokak output folder).")

    optional.add_argument("-p","--pharokka_directory", dest="pharokka_dir", action='store', 
                            default="/mnt/c/Users/Alvaro/Desktop/projects/phage/data/interim/pharokka/pharokka_full_output/", help="Path where pharokka output will be stored (default: /mnt/c/Users/Alvaro/Desktop/projects/phage/data/interim/pharokka/pharokka_full_output/).")
    optional.add_argument("-b","--database_directory", dest="pharokka_database", action='store', 
                            default="/mnt/c/Users/Alvaro/Desktop/pharokka_database/", help="Path to the Pharokka database (default./mnt/c/Users/Alvaro/Desktop/pharokka_database/).")
    optional.add_argument("-o", "--output", dest="output", action='store', default="/mnt/c/Users/Alvaro/Desktop/projects/phage/data/processed/test_engineer.csv", 
                    help="Path to the raw data. Accepted input is a FASTA file with phage sequences.")    
    args = parser.parse_args()

    # Call main function with parsed arguments
    main(args.data, args.output, args.pharokka_dir, args.pharokka_database)
