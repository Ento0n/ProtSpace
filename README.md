# RostSpace

RostSpace visualizes embeddings in 3D and 2D space and colors them by groups, provided throug a CSV file. Dimensionality reduction can be performed with UMAP, PCA or t-SNE. Additional features for biologists are implemented, to display molecule structures by providing PDB files.

## Installation
RostSpace can be installed with `pip`:
```shell
# install it with pip
pip install rostspace
```

## Getting started
After installation, RostSpace can be run using the `rostspace` command.

The required arguments are:

    -o      Path to the output directory where generated files are stored
    --hdf   Path to HDF5-file containing the per protein embeddings as key-value pair
    --csv   Path to the .csv-file containing the metadata

Optional arguments are:

    --pdb   Path to the directory that holds the .pdb files
    --fasta Path to the .fasta-file
    -v      Verbose -- prints internal process to the terminal

### Download example data
Example data can be downloaded from [here](https://nextcloud.in.tum.de/index.php/s/BPWWA9tiXTawjjW).

### Run RostSpace
RostSpace can be run by providing arguments on the command line, by giving it a YAML file, or a combination of both.

1. Giving arguments on the command line:
    ```shell
    rostspace -o data/KLK --hdf data/KLK/KLK_esm2.h5 --csv data/KLK/KLK.csv
    ```

2. Creating a YAML file:

    Pla2g2.yaml
    ```yaml
    o: data/Pla2g2
    hdf: data/Pla2g2/emb_esm2.h5
    csv: data/Pla2g2/Pla2g2.csv
    ```

    ```bash
    rostspace -conf conf/Pla2g2.yaml
    ```

3. Using a YAML file and giving extra arguments:
    ```bash
    rostspace -conf conf/Pla2g2.yaml --pdb data/Pla2g2/colabfold/pdb
    ```


For more information to the arguments run
```shell
rostspace --help
```
