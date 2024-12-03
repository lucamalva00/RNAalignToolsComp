# RNAalignToolsComp - A Comparison of RNA 3D structural alignment tools

This project is developed for the Group Project exam at UNICAM (Università degli Studi di Camerino). The purpose of this project is comparing six RNA structural alignment tools (RMalign, USalign, STAlign, SARA, CLICK and ARTS) and evaluate the differences on alignment results using unsupervised ML methods.

## Project Overview

We designed and conducted an experiment to classify various tertiary structure alignment tools using unsupervised Machine Learning algorithms. To achieve this, we selected several tools and created a dataset to test them. We then used the results to perform our classification. This type of comparison is particularly relevant because, despite the wide availability of structural alignment software, different tools use different algorithms and can produce varying results depending on the molecules and algorithms involved. The main goal of the project is to analyze the behavior of STAlign (developed by the University of Camerino) in relation to the other tools selected for the experiment. We aim to assess how different alignment approaches, based on various algorithms, influence clustering results.

## Configuration

Since all the alignment tools support Linux Operating System, a version of Linux OS needs to be installed to run this tool.
Before using RNAalignToolsComp, the following needs to be installed to be able to run all the alignment tools:

- Python 3 ([Python 3.12.4](https://www.python.org/downloads/release/python-3124/))
- Python 2.7 ([Python 2.7.0](https://www.python.org/download/releases/2.7/))
- Java SE Runtime Environment 8 or later ([Java 21](https://www.oracle.com/it/java/technologies/downloads/#java21))
- Biopython 1.76 and 1.84 ([Biopython download](https://biopython.org/wiki/Download))
- Numpy 2.0.1 and 1.16 ([Numpy releases](https://numpy.org/doc/stable/release.html))
- scikit-learn 1.5 ([scikit-learn installation](https://scikit-learn.org/stable/install.html#install-official-release))

## Inputs and Outputs

All the tools operate their alignments on a given dataset consisting of any number of molecules: the molecule file type supported by this tools is the **PDB format (.pdb)**.
The all-to-all alignment, performed on the given dataset by all the structural alignment tools, produces a csv file representing the alignments results: each row follows this format: 
```
molecule a, molecule b, distance
```
After performing the all-to-all alignments, the csv files are used to produce clusters for each tool: the result plots obtained from the Principal Component Analysis (PCA) for each alignment tool are displayed.

## Usage
```
> python3 build_dataset.py pdb_IDs_file.txt
```

This command is used to produce a dataset from a text file containg molecules IDs (a comma separated list). The `Dataset` folder containing all the PDB files obtained from the IDs on the given file is created.

> [!WARNING]
> The `Dataset` folder created must be manually moved under the path `./Datasets/`, otherwise other scripts won't work. If the `Datasets` folders contains other Datasets, could be necessary to rename the newly created dataset folder.

```
> python3 secondary_structures_from_dataset.py dataset_name
```

This command is used to derive 2D structures of PDB files in a given dataset: 2D structures are needed for performing alignments using the ARTS tool.

```
> python3 basepairs_from_secondary_structures.py dataset_name
```

Formats the 2D structures of all molecules in the given dataset. Generates `.pairs` (2D structure file format required by ARTS tool) files for all PDBs in the dataset.

> [!WARNING]
> This command must be issued only after the generation of 2D structures.

```
> python3 run_tool_alignments.py tool_name dataset_name
```

This command is used to perform all-to-all alignment on the given dataset and using the given tool. Generates a `.csv` file that contains the results of the alignment.

```
> python3 RNA_tools_clusters.py
```

Performs clustering on alignments results produced by the alignments tools: uses three clustering methods (KMeans, GaussianMixture and BisectingKMeans) and calculates, for each method, two metrics: the Silhouette Score and the Davies-Bouldin Index. It also performs PCA (Principal Component Analysis) for each combination of tool and clustering method, producing 2D plots representing visualy the cluster obtained.

# Credits

Authors:
- Antonio Addesa
- Luca Malvatani

Supervisors:
- Luca Tesei
- Michela Quadrini
