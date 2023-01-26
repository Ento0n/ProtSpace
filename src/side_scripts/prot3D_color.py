"""
Sample code for residue coloring in a dash-bio 3d molecule view

Requirements:
- dash
- dash-bio
- numpy
- biopython

pip install dash dash-bio numpy biopython
"""

import dash
import dash_bio
import numpy as np

from dash_bio.utils import PdbParser, create_mol3d_style
from Bio.PDB.Polypeptide import protein_letters_3to1


def rgb2hex(r, g, b):
    return '#{:02X}{:02X}{:02X}'.format(r, g, b)


# Load and parse PDB file from url
parser = PdbParser('https://alphafold.ebi.ac.uk/files/AF-L8XZM1-F1-model_v4.pdb')
data = parser.mol3d_data()

# Create amino acid sequence as one-letter code string
prot_seq = ''.join(protein_letters_3to1[residue.name] for residue in parser.structure.residues)

# Create 3D mol view for our data
styles = create_mol3d_style(
    data['atoms'], visualization_type='sphere', color_element='atom'
)

# Coloring per residue as uint8 RGB. Transparency does not seem to be supported yet.
# We use random colors here for demonstration purposes.
color_per_residue = np.random.randint(0, 256, (len(prot_seq), 3))

for i, (style, atom) in enumerate(zip(styles, data['atoms'])):
    # Note: we iterate over all atoms but have one color per residue.
    # We therefore get the residue index for every atom.
    residue_idx = atom['residue_index']

    style['color'] = rgb2hex(*color_per_residue[residue_idx])


app = dash.Dash(__name__)
app.layout = dash_bio.Molecule3dViewer(
    id='mol-viewer',
    modelData=data,
    styles=styles
)


if __name__ == '__main__':
    app.run_server()