from sklearn.cluster import *
from sklearn.metrics import *
import pandas as pd
import numpy as np

# /**
# * Calculates and prints evaluation metrics for a given clustering model.
# *
# * @param model_name The name of the clustering model to use ('KMeans', 'SpectralClustering', 'BisectingKMeans').
# * @param features The feature matrix (distances between molecules) to apply clustering on.
# *
# * Models:
# * - KMeans: A clustering algorithm that tries to minimize the distance of points from the centroid of each cluster.
# * - SpectralClustering: A method based on spectral decomposition that identifies non-linear clusters.
# * - BisectingKMeans: A variant of KMeans that applies iterative splitting of clusters.
# *
# * Metrics:
# * - Silhouette Score: Measures how well each point lies within its cluster compared to other clusters.
# * - Davies-Bouldin Index: Measures the compactness and separation of clusters. Lower values indicate better-defined clusters.
# */
def metrics_calculator(model_name, features):
    model = None
    match model_name:
        case 'KMeans':
            # Apply the KMeans algorithm with 20 clusters
            model = KMeans(n_clusters=20)
        case 'SpectralClustering':
            # Apply Spectral Clustering with 20 clusters and a precomputed affinity matrix
            model = SpectralClustering(n_clusters=20, affinity='precomputed')
        case 'BisectingKMeans':
            # Apply the Bisecting KMeans algorithm with 20 clusters
            model = BisectingKMeans(n_clusters=20)

    # Perform clustering and get the cluster labels for each point
    labels = model.fit_predict(features)

    # Calculate Silhouette Score
    score = silhouette_score(features, labels)

    # Calculate Davies-Bouldin Index
    db_index = davies_bouldin_score(features, labels)

    # Print the results
    print(f'Model: {model_name}')
    print(f'Silhouette Score: {score}')
    print(f'Davies-Bouldin Index: {db_index}\n\n\n')


# /**
# * Builds a matrix of feature vectors based on a CSV file containing distances between pairs of molecules.
# *
# * @param file_name The name of the CSV file containing molecule distances.
# * @return A list of vectors where each vector represents the distances of a molecule to other molecules.
# *
# * The CSV is composed of three columns: mol1 (first molecule), mol2 (second molecule), and distance (the distance between the two molecules).
# * A symmetric matrix is built where rows and columns represent molecules, and the values represent the distances between them.
# */
def vectors_builder(file_name):
    # Read the CSV file and create a DataFrame with three columns: mol1, mol2, and distance
    molecules = pd.read_csv('Pdb_ARTS_alignments.csv', header=None)
    molecules.columns = ["mol1", "mol2", "distance"]

    # Get a unique list of all molecules
    unique_molecules = pd.unique(molecules[['mol1', 'mol2']].values.ravel())

    # Create a distance matrix (initialized with NaN)
    distance_matrix = pd.DataFrame(index=unique_molecules, columns=unique_molecules, data=np.nan)

    # Populate the matrix with distance values
    for k in range(len(molecules)):
        mol1 = molecules.loc[k, 'mol1']
        mol2 = molecules.loc[k, 'mol2']
        dist = molecules.loc[k, 'distance']

        # Fill the symmetric matrix
        distance_matrix.loc[mol1, mol2] = dist
        distance_matrix.loc[mol2, mol1] = dist

    # Replace NaN values with 0 (assume missing distances are zero)
    distance_matrix.fillna(0, inplace=True)

    # print(distance_matrix)

    # Convert each row into a vector, where each row represents a molecule
    features_list = distance_matrix.values  # Each row represents a vector for one molecule
    # print(vectors)
    return features_list


# List of CSV files containing molecular alignment data
csv_list = ['Pdb_ARTS_alignments.csv', 'Pdb_CLICK_alignments.csv', 'Pdb_RMalign_alignments.csv',
            'Pdb_USalign_alignments.csv', 'Pdb_STalign_alignments.csv', 'Pdb_SARA_alignments.csv']

# List of clustering models to be tested
model_names = ['KMeans', 'SpectralClustering', 'BisectingKMeans']

# List that will contain feature vectors for each dataset
vectors_list = []

# For each CSV file, build the feature vectors and add them to the list
for csv_name in csv_list:
    vectors_list.append(vectors_builder(csv_name))

# Perform clustering for each CSV file and for each model
index = 0
for vectors in vectors_list:
    print(csv_list[index] + ':\n')
    index += 1
    for name in model_names:
        # Calculate and print the metrics for the current model on the dataset's feature vectors
        metrics_calculator(name, vectors)
