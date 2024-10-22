import random
import math
from collections import Counter

#########################################################################
###---------------3.MOTIF DISCOVERY - GIBBS SAMPLING------------------###
#########################################################################


def read_sequence_from_file(file_path):
    with open(file_path, 'r') as file:
        # Skip the first line (metadata)
        file.readline()
        # Remove newline characters and read the remaining content
        sequence = file.read().replace('\n', '').strip()
    return sequence



def initialize_motifs(sequences, motif_length):
    """
    Initializes motifs randomly from input sequences.

    Parameters:
    - sequences: List of DNA sequences.
    - motif_length: Length of the motifs to be initialized.

    Returns:
    List of randomly initialized motifs.
    """
    motifs = []

    # Iterate over each sequence to randomly select a motif
    for sequence in sequences:
        # Generate a random start position within the sequence
        start = random.randint(0, len(sequence) - motif_length)

        # Extract the motif based on the random start position and motif length
        motif = sequence[start:start + motif_length]

        # Append the motif to the list
        motifs.append(motif)

    return motifs



def create_profile_matrix(motifs):
    """
    Creates a profile matrix from a list of motifs.

    Parameters:
    - motifs: List of motifs (DNA sequences).

    Returns:
    Profile matrix represented as a list of dictionaries.
    """
    motif_length = len(motifs[0]) if motifs else 0
    profile_matrix = []

    # Iterate over each position in the motifs
    for i in range(motif_length):
        # Extract the column (nucleotide at the current position) from all motifs
        column = [motif[i] for motif in motifs if len(motif) > i]

        # Count occurrences of each nucleotide in the column
        counts = Counter(column)

        # Calculate probabilities for each nucleotide
        probabilities = {base: count / len(column) for base, count in counts.items()}

        # Append the probabilities to the profile matrix
        profile_matrix.append(probabilities)

    return profile_matrix



def calculate_probabilities(sequence, profile_matrix, motif_length, pseudocount=1):
    """
    Calculates probabilities for each position in a sequence based on a profile matrix.

    Parameters:
    - sequence: Input DNA sequence.
    - profile_matrix: Profile matrix generated from motifs.
    - motif_length: Length of the motifs.
    - pseudocount: Pseudocount value to handle zero probabilities.

    Returns:
    List of probabilities for each position in the sequence.
    """
    probabilities = []

    # Iterate over possible starting positions in the sequence
    for i in range(len(sequence) - motif_length + 1):
        motif = sequence[i:i + motif_length]
        probability = 1.0

        # Iterate over each position in the motif
        for j in range(motif_length):
            nucleotide = motif[j]

            # Check if the nucleotide is in the profile matrix
            if nucleotide in profile_matrix[j]:
                probability *= profile_matrix[j][nucleotide]
            else:
                # Pseudocount to handle zero probability
                probability *= pseudocount / (len(profile_matrix[j]) + pseudocount * 4)

        probabilities.append(probability)

    total_probability = sum(probabilities)

    # Check if total_probability is zero before division
    return [prob / total_probability if total_probability != 0 else 0.0 for prob in probabilities]



def score_motifs(motifs, background_probabilities):
    """
    Calculates the score of a set of motifs using a background probability matrix.

    Parameters:
    - motifs: List of motifs to be scored.
    - background_probabilities: Background probability matrix.

    Returns:
    The calculated score for the motifs.
    """
    score = 0
    motif_length = len(motifs[0])
    num_motifs = len(motifs)

    # Iterate over each position in the motifs
    for i in range(motif_length):
        column = [motif[i] for motif in motifs]
        counts = Counter(column)

        # Iterate over unique bases in the column
        for base, count in counts.items():
            probability = count / num_motifs
            background_probability = background_probabilities[i][base]  # Use integer index

            # Calculate score using information content formula
            if probability > 0:
                score += count * math.log2(probability / background_probability)

    return score



def gibbs_sampler(sequences, motif_length, iterations):
    """
    Performs Gibbs sampling for motif discovery.

    Parameters:
    - sequences: List of DNA sequences.
    - motif_length: Length of the motifs to be discovered.
    - iterations: Number of iterations for the Gibbs sampling.

    Returns:
    The best motifs discovered by the Gibbs sampler.
    """
    motifs = initialize_motifs(sequences, motif_length)
    best_motifs = motifs.copy()

    # Create background probabilities using all sequences
    background_probabilities = create_profile_matrix(sequences)

    best_score = score_motifs(best_motifs, background_probabilities)

    # Perform Gibbs sampling iterations
    for _ in range(iterations * 10):
        i = random.randint(0, len(sequences) - 1)
        motif_i = motifs[i]
        other_motifs = motifs[:i] + motifs[i + 1:]

        # Update the profile matrix using other motifs
        profile_matrix = create_profile_matrix(other_motifs)
        probabilities = calculate_probabilities(sequences[i], profile_matrix, motif_length)

        # Check if all probabilities are zero
        if all(prob == 0 for prob in probabilities):
            continue  # Skip this iteration if all probabilities are zero

        # Choose a new motif start position based on probabilities
        new_motif_start = random.choices(range(len(sequences[i]) - motif_length + 1), weights=probabilities)[0]
        new_motif = sequences[i][new_motif_start:new_motif_start + motif_length]

        motifs[i] = new_motif

        current_score = score_motifs(motifs, background_probabilities)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs.copy()

    return best_motifs
