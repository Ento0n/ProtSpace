# RostSpace

RostSpace is an embedding visualisation tool that allows to display any embeddings in a 3D or 2D graph. 
For the dimensionality reduction UMAP, PCA and t-SNE are selectable.
Additionally, features for biologists are implemented. It allows to display molecule structures or 
sequences by providing a corresponding .pdb- or .fasta-file.

## Running the script

After installing the tool with pypi, it can be started with

```shell
rostspace
```

The required arguments are:

    -o      Path to the output directory where all generated files are stored
    --hdf   Path to HDF5-file containing the per protein embeddings as a key-value pair
    --csv   Path to the .csv-file containing the metadata

Optional arguments are:

    --fasta Path to the .fasta-file
    --pdb   Path to the directory that holds the .pdb-files
    -v      Takes no value, if set the tool prints its internal process in the terminal

Example:

```shell
rostspace -o data/Pla2g2 --hdf data/Pla2g2/Pla2g2_prott5.h5 --csv data/Pla2g2/Pla2g2.csv
```

Alternatively the arguments can be provided with a configuration file in PyYaml format:

Pla2g2.yaml
```shell
o: data/Pla2g2
hdf: data/Pla2g2/emb_esm2.h5
csv: data/Pla2g2/Pla2g2.csv
```

Example:

```shell
rostspace -conf conf/Pla2g2.yaml
```
