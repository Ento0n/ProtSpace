# ProtSpace3D

This is a bachelor thesis project, 
Development of protein-embedding visualization tool.

## Installing dependencies

Python-poetry(https://python-poetry.org/) is used for installing the dependencies. Follow the instruction 
on the website to install poetry.
After that run

```shell
poetry install
```
to install the dependencies for this project.

## Running the script

The script to be executed is processing.py with the arguments:

    ->  -d          Name of the folder which holds the required data, .h5 .csv & .fasta (String)
    ->  -b          Name of the files which are in the data folder, requires equal names (String)
    ->  --sep       The character which seperates the columns in the .csv file (Character)
    ->  --uid_col   The column number which holds the unique ID, starting from 0 (Integer)
    ->  --html_col  If set, html file of the selected column is saved in data directory (Integer)

Example:

```shell
python app.py -d data -b VA --sep , --uid_col 0
```
