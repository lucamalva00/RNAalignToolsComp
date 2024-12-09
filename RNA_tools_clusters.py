from sklearn.cluster import KMeans, BisectingKMeans, SpectralClustering, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def elbow_plot_generator(data, max_k, dataset_name):
    means = []
    inertias = []

    for k in range(1, max_k):
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(data)
        means.append(k)
        inertias.append(kmeans.inertia_)

    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(dataset_name, fontsize=16)
    plt.plot(means, inertias)
    plt.xlabel('Number of clusters')
    plt.ylabel('Inertia')
    plt.grid(True)
    plt.show()

# Calculate and print evaluation metrics for a given clustering model.

# @param model_name: The name of the clustering model to use ('KMeans', 'GaussianMixture', 'BisectingKMeans').
# @param features: The feature matrix (distances between molecules) on which to apply clustering.

# Models:
# - KMeans: A clustering algorithm that seeks to minimize the distance of points from the centroid of each cluster.
# - AgglomerativeClustering: a hierarchical, bottom-up, clustering method.
# - BisectingKMeans: A variant of KMeans that applies an iterative split of the clusters.

# Metrics:
# - Silhouette Score: Measures how well each point is located within its own cluster compared to other clusters.
# - Davies-Bouldin Index: Measures the compactness and separation of clusters. Lower values indicate better defined clusters.
def metrics_calculator(model_name, features, n_clusters):
    clustering_model = None
    match model_name:
        case 'KMeans':
            clustering_model = KMeans(n_clusters=n_clusters)
        case 'AgglomerativeClustering':
            clustering_model = AgglomerativeClustering(n_clusters=n_clusters)
        case 'BisectingKMeans':
            clustering_model = BisectingKMeans(n_clusters=n_clusters)

    # Obtain the cluster labels
    labels = clustering_model.fit_predict(features)

    score = silhouette_score(features, labels)
    db_index = davies_bouldin_score(features, labels)


    print(f'Model: {model_name}')
    print(f'Silhouette Score: {score}')
    print(f'Davies-Bouldin Index: {db_index}\n\n\n')

    return clustering_model

# Constructs a features vector matrix based on a CSV file containing distances between pairs of molecules.

# @param file_name: The name of the CSV file containing the molecular distances.
# @return: A list of vectors, where each vector represents the distances of a molecule to other molecules.

# The CSV is composed of three columns: mol1 (first molecule), mol2 (second molecule) and distance (distance between the two molecules).
# A symmetric matrix is constructed where rows and columns represent molecules, and values represent the distances between them.

def vectors_builder(file_name):
    scaler = StandardScaler()

    # Read the CSV file and create a DataFrame
    molecules_df = pd.read_csv(file_name, header=None)

    # Print the DataFrame's columns to debug
    print(f"Columns in {file_name}: {molecules_df.columns.tolist()}")

    # If there are more than 3 columns, drop the extra columns
    if molecules_df.shape[1] > 3:
        print(f"Warning: Found {molecules_df.shape[1]} columns, dropping the extra columns.")
        molecules_df = molecules_df.iloc[:, :3]  # Keep only the first three columns

    # Ensure the correct number of columns before assigning names
    molecules_df.columns = ["mol1", "mol2", "distance"]

    # Get a unique list of all molecules
    unique_molecules = pd.unique(molecules_df[['mol1', 'mol2']].values.ravel())

    # Create a distance matrix (initialized with NaN)
    distance_matrix = pd.DataFrame(index=unique_molecules, columns=unique_molecules, data=np.nan)

    # Populate the matrix with distance values
    for k in range(len(molecules_df)):
        mol1 = molecules_df.loc[k, 'mol1']
        mol2 = molecules_df.loc[k, 'mol2']
        dist = molecules_df.loc[k, 'distance']

        # Fill the symmetric matrix
        distance_matrix.loc[mol1, mol2] = dist
        distance_matrix.loc[mol2, mol1] = dist

    # Replace NaN values with 0 (missing distances are the distances from one molecule from itself)
    distance_matrix.fillna(0, inplace=True)

    # Convert each row into a vector, where each row represents a molecule and use scaler
    # to standardize the data
    features_list = scaler.fit_transform(distance_matrix.values)
    return features_list


csv_list = ['Pdb_ARTS_alignments.csv', 'Pdb_CLICK_alignments.csv', 'Pdb_RMalign_alignments.csv',
            'Pdb_USalign_alignments.csv', 'Pdb_STalign_alignments.csv', 'Pdb_SARA_alignments.csv']

model_names = ['KMeans', 'AgglomerativeClustering', 'BisectingKMeans']

# Contains features vectors for all datasets
vectors_list = []


for csv_name in csv_list:
    vectors_list.append(vectors_builder(csv_name))

index = 0
for vectors in vectors_list:
    elbow_plot_generator(vectors, 10, csv_list[index])
    index += 1


num_of_clusters = int(input('Select number of clusters:\n'))

index = 0
for vectors in vectors_list:
    print(csv_list[index] + ':\n')
    index += 1

    for model_name in model_names:
        model = metrics_calculator(model_name, vectors, num_of_clusters)

        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(vectors)
        principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

        labels = model.labels_

        # Create a separate plot for each tool and method
        plt.figure(figsize=(8, 8))
        for label in set(labels):
            indices = labels == label
            plt.scatter(principalDf.loc[indices, 'PC1'], principalDf.loc[indices, 'PC2'], label=f'Cluster {label}')

        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title(f'Tool: {csv_list[index - 1]}, Model: {model_name}')
        plt.legend()
        plt.grid(True)
        plt.show()
