from sklearn.cluster import KMeans, BisectingKMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture

#
# Calcola e stampa le metriche di valutazione per un dato modello di clustering.
#
# @param model_name: Il nome del modello di clustering da utilizzare ('KMeans', 'GaussianMixture', 'BisectingKMeans').
# @param features: La matrice delle feature (distanze tra le molecole) su cui applicare il clustering.
#
# Modelli:
# - KMeans: Algoritmo di clustering che cerca di minimizzare la distanza dei punti dal centroide di ciascun cluster.
# - GaussianMixture: Un modello statistico che assume che i dati siano generati da una combinazione di distribuzioni gaussiane.
# - BisectingKMeans: Una variante di KMeans che applica uno split iterativo dei cluster.
#
# Metriche:
# - Silhouette Score: Misura quanto bene ogni punto si trova all'interno del proprio cluster rispetto agli altri cluster.
# - Davies-Bouldin Index: Misura la compattezza e la separazione dei cluster. Valori più bassi indicano cluster meglio definiti.
#
def metrics_calculator(model_name, features):
    model = None
    match model_name:
        case 'KMeans':
            # Applica l'algoritmo KMeans con 20 cluster
            model = KMeans(n_clusters=20)
        case 'GaussianMixture':
            # Applica il modello Gaussian Mixture con 20 componenti
            model = GaussianMixture(n_components=20)
        case 'BisectingKMeans':
            # Applica l'algoritmo Bisecting KMeans con 20 cluster
            model = BisectingKMeans(n_clusters=20)

    # Esegui il clustering e ottieni le etichette dei cluster per ciascun punto
    labels = model.fit_predict(features)

    # Calcola il Silhouette Score
    score = silhouette_score(features, labels)

    # Calcola l'Indice di Davies-Bouldin
    db_index = davies_bouldin_score(features, labels)

    # Stampa i risultati
    print(f'Model: {model_name}')
    print(f'Silhouette Score: {score}')
    print(f'Davies-Bouldin Index: {db_index}\n\n\n')


#
# Costruisce una matrice di vettori di feature basata su un file CSV contenente le distanze tra coppie di molecole.
#
# @param file_name: Il nome del file CSV contenente le distanze delle molecole.
# @return: Una lista di vettori, dove ciascun vettore rappresenta le distanze di una molecola rispetto ad altre molecole.
#
# Il CSV è composto da tre colonne: mol1 (prima molecola), mol2 (seconda molecola) e distance (distanza tra le due molecole).
# Viene costruita una matrice simmetrica in cui righe e colonne rappresentano le molecole, e i valori rappresentano le distanze tra di esse.
#
def vectors_builder(file_name):
    # Read the CSV file and create a DataFrame
    molecules = pd.read_csv(file_name, header=None)

    # Print the DataFrame's columns to debug
    print(f"Columns in {file_name}: {molecules.columns.tolist()}")

    # If there are more than 3 columns, drop the extra columns
    if molecules.shape[1] > 3:
        print(f"Warning: Found {molecules.shape[1]} columns, dropping the extra columns.")
        molecules = molecules.iloc[:, :3]  # Keep only the first three columns

    # Ensure the correct number of columns before assigning names
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

    # Replace NaN values with 0 (assuming missing distances are zero)
    distance_matrix.fillna(0, inplace=True)

    # Convert each row into a vector, where each row represents a molecule
    features_list = distance_matrix.values  # Each row represents a vector for one molecule
    return features_list


# Lista di file CSV contenenti dati di allineamento molecolare
csv_list = ['Pdb_ARTS_alignments.csv', 'Pdb_CLICK_alignments.csv', 'Pdb_RMalign_alignments.csv',
            'Pdb_USalign_alignments.csv', 'Pdb_STalign_alignments.csv', 'Pdb_SARA_alignments.csv']

# Lista dei modelli di clustering da testare
model_names = ['KMeans', 'GaussianMixture', 'BisectingKMeans']

# Lista che conterrà i vettori di feature per ciascun dataset
vectors_list = []

# Per ciascun file CSV, costruisci i vettori di feature e aggiungili alla lista
for csv_name in csv_list:
    vectors_list.append(vectors_builder(csv_name))

# Esegui il clustering per ciascun file CSV e per ciascun modello
index = 0
for vectors in vectors_list:
    print(csv_list[index] + ':\n')
    index += 1
    for name in model_names:
        # Calcola e stampa le metriche per il modello corrente sui vettori di feature del dataset
        metrics_calculator(name, vectors)
