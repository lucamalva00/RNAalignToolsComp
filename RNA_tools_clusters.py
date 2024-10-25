from sklearn.cluster import *
from sklearn.metrics import *
import pandas as pd
import numpy as np

#Nel modello k-means, i cluster vengono calcolati attraverso
# un processo iterativo che cerca di minimizzare la distanza tra i punti
# dati e il centroide di ciascun cluster. Il procedimento generale è composto da
# pochi passi principali

def metrics_calculator(model_name, features):
    model = None
    match model_name:
        case 'KMeans':
            model = KMeans(n_clusters=20)
        case 'SpectralClustering':
            model = SpectralClustering(n_clusters=20, affinity='precomputed')
        case 'BisectingKMeans':
            model = BisectingKMeans(n_clusters=20)
    labels = model.fit_predict(features)
    # SIlhouette Score
    score = silhouette_score(features, labels)
    # Davies - Bouldin score
    db_index = davies_bouldin_score(features, labels)
    print(f'Model: {model_name}\n')
    print(f'Silhouette Score: {score}')
    print(f'Davies-Bouldin Index: {db_index}\n')

def vectors_builder(file_name):
    molecules = pd.read_csv('Pdb_ARTS_alignments.csv', header=None)
    molecules.columns = ["mol1", "mol2", "distance"]

    # Ottieni una lista unica di tutte le molecole
    unique_molecules = pd.unique(molecules[['mol1', 'mol2']].values.ravel())

    # Crea una matrice di distanze (inizializzata a NaN)
    distance_matrix = pd.DataFrame(index=unique_molecules, columns=unique_molecules, data=np.nan)

    # Popola la matrice con i valori di distanza
    for k in range(len(molecules)):
        mol1 = molecules.loc[k, 'mol1']
        mol2 = molecules.loc[k, 'mol2']
        dist = molecules.loc[k, 'distance']

        # Riempie la matrice simmetrica
        distance_matrix.loc[mol1, mol2] = dist
        distance_matrix.loc[mol2, mol1] = dist

    # Sostituisci NaN con 0
    distance_matrix.fillna(0, inplace=True)

    #print(distance_matrix)

    # Converti ogni riga in un vettore
    features_list = distance_matrix.values  # Ogni riga rappresenta un vettore per una molecola
    #print(vectors)
    return features_list


csv_list = ['Pdb_ARTS_alignments.csv', 'Pdb_CLICK_alignments.csv', 'Pdb_RMalign_alignments.csv', 'Pdb_USalign_alignments.csv', 'Pdb_STalign_alignments.csv', 'Pdb_SARA_alignments.csv']
model_names = ['KMeans', 'SpectralClustering', 'BisectingKMeans']
vectors_list = []

for csv_name in  csv_list:
    vectors_list.append(vectors_builder(csv_name))


index = 0
for vectors in vectors_list:
    print(csv_list[index] + ':\n')
    index+=1
    for name in model_names:
        metrics_calculator(name, vectors)



#for cluster_num in range(model.n_clusters):
#    print(f"Cluster {cluster_num}:")
#    print(vectors[model.labels_==cluster_num])




