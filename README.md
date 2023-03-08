# RostSpace

RostSpace is a embedding visualisation tool specialized for biologists, but not limited by this. It
allows to display embeddings in a 3D or 2D graph.

## Running the script

After installing the tool with pypi, it can be started with

```shell
rostspace
```

The required arguments are:

    -o      Path to the output directory where all generated files are stored
    --hdf   Path to HDF5-file containing the per protein embeddings as a key-value pair
    --csv   Path to the csv-file containing the metadata

Optional arguments are:

    --pdb   Path to the directory that holds the .pdb-files
    -v      Takes no value, if set the tool prints its internal process in the terminal

Example:

```shell
rostspace -o data/Pla2g2 --hdf data/Pla2g2/Pla2g2_prott5.h5 --csv data/Pla2g2/Pla2g2.csv
```

Alternatively the arguments can be provided with a configuration file in yaml format.

Example:

```shell
rostspace -conf conf/Pla2g2.yaml
```
