import pandas as pd

#########################################################################
###---------------------------1.DATA CLEANING-------------------------###
#########################################################################

def data_cleaning():

    # The dataset is stored in a file named 'E-MTAB-8626_tpms.tsv'
    file_path = 'E-MTAB-8626_tpms.tsv'

    # Read the dataset into a pandas DataFrame, considering the irregularities in the number of columns
    df = pd.read_csv(file_path, sep='\t', skiprows=5)

    # Display the original dataset
    print("Original Dataset:")
    print(df)

    # Data Cleaning: Remove rows with missing values
    df_cleaned = df.dropna()

    # Display the cleaned dataset
    print("\nCleaned Dataset:")
    print(df_cleaned)

    # Save the cleaned dataset to a new file
    cleaned_file_path = 'cleaned_data.csv'
    df_cleaned.to_csv(cleaned_file_path, index=False)
    print(f"\nCleaned data saved to {cleaned_file_path}")