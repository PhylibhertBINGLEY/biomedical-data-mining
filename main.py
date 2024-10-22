import pandas as pd
import numpy as np
from custom_packages import preprocessing, clustering_kmeans, motif_discovery


#########################################################################
###------------------------------MAIN---------------------------------###
#########################################################################

def main():

    ###------------------------Data cleaning------------------------##

    preprocessing.data_cleaning()
    print('data cleaning done')
    print('----------------------------------')


    ###-------------------------Clustering--------------------------##
    #print('clustering already performed with k=4')

    # Load cleaned data from the CSV file
    file_path = 'cleaned_data.csv'
    df_cleaned = pd.read_csv(file_path)

    # Extract numeric data for clustering (assuming expression data starts from column 2)
    data = df_cleaned.iloc[:, 2:].values

    # Number of clusters (we can adjust the value)
    k = 4

    ## Run K-means clustering
    labels, centroids = clustering_kmeans.k_means(data, k)

    ## Display the cluster labels and centroids
    print("Cluster Labels:")
    print(labels)
    print("\nCentroids:")
    print(centroids)
    

    # Count the number of members in each cluster
    unique_labels, counts = np.unique(labels, return_counts=True)
    cluster_sizes = dict(zip(unique_labels, counts))
    # Print the number of members in each cluster
    print("\nNumber of Members in Each Cluster:")
    for cluster, size in cluster_sizes.items():
        print(f"Cluster {cluster}: {size} members")
        

    ## Add the cluster assignment column to the cleaned DataFrame
    df_cleaned['Cluster_Label'] = labels
    # Save the DataFrame with the cluster assignment to a new file
    df_cleaned.to_csv('cleaned_data_with_clusters.csv', index=False, header=False)


    print('Clustering done')
    print('----------------------------------')

    ###-------------------------Genes extraction--------------------------##

    print('Genes extraction done: check the file gene_names_output.txt')

    ## choose the cluster (choosen size after performing some tests on the value of k: 102)
    # and then extract the gene names from that cluster
    # Load data with clustering labels from the CSV file
    file_path_bis = 'cleaned_data_with_clusters.csv'
    df_cleaned_with_labels = pd.read_csv(file_path_bis, header=None)
    # Specify the output file path
    output_file_path = 'gene_names_output.txt'
    # Call the function
    clustering_kmeans.cluster_gene_extraction(df_cleaned_with_labels, labels, output_file_path)

    print('Genes extraction done')
    print('----------------------------------')

    ###-------------------------Motif discovery--------------------------##

    # Replace 'your_sequence_files' with the actual paths to your 102 sequence files
    sequence_paths = ['./genes_sequences/S288C_YAL003W_EFB1_flanking.fsa', './genes_sequences/S288C_YAL038W_CDC19_flanking.fsa', './genes_sequences/S288C_YBL027W_RPL19B_flanking.fsa','./genes_sequences/S288C_YBR011C_IPP1_flanking.fsa','./genes_sequences/S288C_YBR084C-A_RPL19A_flanking.fsa','./genes_sequences/S288C_YBR189W_RPS9B_flanking.fsa','./genes_sequences/S288C_YBR191W_RPL21A_flanking.fsa','./genes_sequences/S288C_YCR012W_PGK1_flanking.fsa','./genes_sequences/S288C_YCR024C-A_PMP1_flanking.fsa','./genes_sequences/S288C_YCR024C-B_YCR024C-B_flanking.fsa','./genes_sequences/S288C_YCR031C_RPS14A_flanking.fsa','./genes_sequences/S288C_YDL061C_RPS29B_flanking.fsa','./genes_sequences/S288C_YDL075W_RPL31A_flanking.fsa','./genes_sequences/S288C_YDL081C_RPP1A_flanking.fsa','./genes_sequences/S288C_YDL083C_RPS16B_flanking.fsa','./genes_sequences/S288C_YDL130W_RPP1B_flanking.fsa','./genes_sequences/S288C_YDR025W_RPS11A_flanking.fsa','./genes_sequences/S288C_YDR050C_TPI1_flanking.fsa','./genes_sequences/S288C_YDR064W_RPS13_flanking.fsa','./genes_sequences/S288C_YDR077W_SED1_flanking.fsa','./genes_sequences/S288C_YDR133C_YDR133C_flanking.fsa','./genes_sequences/S288C_YDR155C_CPR1_flanking.fsa','./genes_sequences/S288C_YDR382W_RPP2B_flanking.fsa','./genes_sequences/S288C_YDR447C_RPS17B_flanking.fsa','./genes_sequences/S288C_YDR450W_RPS18A_flanking.fsa','./genes_sequences/S288C_YDR500C_RPL37B_flanking.fsa','./genes_sequences/S288C_YDR502C_SAM2_flanking.fsa','./genes_sequences/S288C_YER091C_MET6_flanking.fsa','./genes_sequences/S288C_YER117W_RPL23B_flanking.fsa','./genes_sequences/S288C_YER131W_RPS26B_flanking.fsa','./genes_sequences/S288C_YFL039C_ACT1_flanking.fsa','./genes_sequences/S288C_YFR032C-A_RPL29_flanking.fsa','./genes_sequences/S288C_YGL008C_PMA1_flanking.fsa','./genes_sequences/S288C_YGL030W_RPL30_flanking.fsa','./genes_sequences/S288C_YGL031C_RPL24A_flanking.fsa','./genes_sequences/S288C_YGL123W_RPS2_flanking.fsa','./genes_sequences/S288C_YGL189C_RPS26A_flanking.fsa','./genes_sequences/S288C_YGR027C_RPS25A_flanking.fsa','./genes_sequences/S288C_YGR034W_RPL26B_flanking.fsa','./genes_sequences/S288C_YGR060W_ERG25_flanking.fsa','./genes_sequences/S288C_YGR118W_RPS23A_flanking.fsa','./genes_sequences/S288C_YGR148C_RPL24B_flanking.fsa','./genes_sequences/S288C_YGR254W_ENO1_flanking.fsa','./genes_sequences/S288C_YHR010W_RPL27A_flanking.fsa','./genes_sequences/S288C_YHR021C_RPS27B_flanking.fsa','./genes_sequences/S288C_YHR141C_RPL42B_flanking.fsa','./genes_sequences/S288C_YIL052C_RPL34B_flanking.fsa','./genes_sequences/S288C_YIL053W_GPP1_flanking.fsa','./genes_sequences/S288C_YIL148W_RPL40A_flanking.fsa','./genes_sequences/S288C_YJL136C_RPS21B_flanking.fsa','./genes_sequences/S288C_YJL159W_HSP150_flanking.fsa','./genes_sequences/S288C_YJL190C_RPS22A_flanking.fsa','./genes_sequences/S288C_YJR104C_SOD1_flanking.fsa','./genes_sequences/S288C_YKL056C_TMA19_flanking.fsa','./genes_sequences/S288C_YKL152C_GPM1_flanking.fsa','./genes_sequences/S288C_YKL180W_RPL17A_flanking.fsa','./genes_sequences/S288C_YKR094C_RPL40B_flanking.fsa','./genes_sequences/S288C_YLR029C_RPL15A_flanking.fsa','./genes_sequences/S288C_YLR048W_RPS0B_flanking.fsa','./genes_sequences/S288C_YLR058C_SHM2_flanking.fsa','./genes_sequences/S288C_YLR075W_RPL10_flanking.fsa','./genes_sequences/S288C_YLR150W_STM1_flanking.fsa','./genes_sequences/S288C_YLR185W_RPL37A_flanking.fsa','./genes_sequences/S288C_YLR249W_YEF3_flanking.fsa','./genes_sequences/S288C_YLR264C-A_YLR264C-A_flanking.fsa','./genes_sequences/S288C_YLR303W_MET17_flanking.fsa','./genes_sequences/S288C_YLR325C_RPL38_flanking.fsa','./genes_sequences/S288C_YLR355C_ILV5_flanking.fsa','./genes_sequences/S288C_YLR388W_RPS29A_flanking.fsa','./genes_sequences/S288C_YLR441C_RPS1A_flanking.fsa','./genes_sequences/S288C_YML024W_RPS17A_flanking.fsa','./genes_sequences/S288C_YML026C_RPS18B_flanking.fsa','./genes_sequences/S288C_YML063W_RPS1B_flanking.fsa','./genes_sequences/S288C_YMR116C_ASC1_flanking.fsa','./genes_sequences/S288C_YMR122W-A_NCW1_flanking.fsa','./genes_sequences/S288C_YMR142C_RPL13B_flanking.fsa','./genes_sequences/S288C_YMR143W_RPS16A_flanking.fsa','./genes_sequences/S288C_YMR230W_RPS10B_flanking.fsa','./genes_sequences/S288C_YMR242C_RPL20A_flanking.fsa','./genes_sequences/S288C_YMR251W-A_HOR7_flanking.fsa','./genes_sequences/S288C_YNL069C_RPL16B_flanking.fsa','./genes_sequences/S288C_YNL145W_MFA2_flanking.fsa','./genes_sequences/S288C_YNL178W_RPS3_flanking.fsa','./genes_sequences/S288C_YNL302C_RPS19B_flanking.fsa','./genes_sequences/S288C_YOL039W_RPP2A_flanking.fsa','./genes_sequences/S288C_YOL040C_RPS15_flanking.fsa','./genes_sequences/S288C_YOL058W_ARG1_flanking.fsa','./genes_sequences/S288C_YOL109W_ZEO1_flanking.fsa','./genes_sequences/S288C_YOL120C_RPL18A_flanking.fsa','./genes_sequences/S288C_YOL127W_RPL25_flanking.fsa','./genes_sequences/S288C_YOR063W_RPL3_flanking.fsa','./genes_sequences/S288C_YOR096W_RPS7A_flanking.fsa','./genes_sequences/S288C_YOR167C_RPS28A_flanking.fsa','./genes_sequences/S288C_YOR182C_RPS30B_flanking.fsa','./genes_sequences/S288C_YOR293W_RPS10A_flanking.fsa','./genes_sequences/S288C_YOR312C_RPL20B_flanking.fsa','./genes_sequences/S288C_YOR369C_RPS12_flanking.fsa','./genes_sequences/S288C_YOR375C_GDH1_flanking.fsa','./genes_sequences/S288C_YPL079W_RPL21B_flanking.fsa','./genes_sequences/S288C_YPL131W_RPL5_flanking.fsa','./genes_sequences/S288C_YPR036W-A_SPO24_flanking.fsa','./genes_sequences/S288C_YPR132W_RPS23B_flanking.fsa']

    # Read gene sequences from files
    sequences = [motif_discovery.read_sequence_from_file(path) for path in sequence_paths]

    # Set motif length and iterations
    motif_length = 7
    iterations = 2000 # multiplied by 10

    # Create background probabilities using all sequences
    background_probabilities = motif_discovery.create_profile_matrix(sequences)

    # Now we can use the 'sequences' list, 'motif_length', 'iterations', and 'background_probabilities' in the Gibbs sampler function
    best_motifs = motif_discovery.gibbs_sampler(sequences, motif_length, iterations)

    # Display the best motifs found
    for i, motif in enumerate(best_motifs):
        print(f"Motif {i + 1}: {motif}")

    # Replace 'T' with 'U' in the motifs
    best_motifs_rna = [motif.replace('T', 'U') for motif in best_motifs]

    # Calculate and display the final score using the new scoring method
    final_score = motif_discovery.score_motifs(best_motifs, background_probabilities)
    print(f"Final Score: {final_score}")


    # Write the best motifs to a text file
    output_file_path = 'best_motifs.txt'
    with open(output_file_path, 'w') as output_file:
        for motif in best_motifs:
            output_file.write(f"{motif}\n\n")

    print(f"Best motifs have been saved to {output_file_path}")


    # Write the best motifs to a text file with 'T' replaced by 'U'
    output_file_path = 'best_motifs_rna.txt'
    with open(output_file_path, 'w') as output_file:
        for motif in best_motifs_rna:
            output_file.write(f"{motif}\n\n")

    print(f"Best motifs (RNA) have been saved to {output_file_path}")
    print('\n')
    print('END OF THE PROGRAM')




if __name__ == '__main__':
    main()







