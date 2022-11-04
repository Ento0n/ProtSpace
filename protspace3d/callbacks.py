#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from visualization import Visualizator

from dash import Input, Output
from dash.exceptions import PreventUpdate
import dash
import dash_bio.utils.ngl_parser as ngl_parser
import pandas as pd

from pandas import DataFrame


def to_mapped_id(original_seq_ids: list, original_id_col: list, df: DataFrame):
    seq_ids = list()

    for id in original_seq_ids:
        index_num = original_id_col.index(id)
        mapped_id = df.index[index_num]
        seq_ids.append(mapped_id)

    return seq_ids


def to_original_id(mapped_seq_ids: list, original_id_col: list, df: DataFrame):
    seq_ids = list()

    for id in mapped_seq_ids:
        index_num = df.index.get_indexer_for([id])[0]
        original_id = original_id_col[index_num]
        seq_ids.append(original_id)

    return seq_ids


def get_callbacks_pdb(app, df, struct_container, original_id_col):
    @app.callback(
        Output("graph", "figure"),
        Input("dd_menu", "value"),
    )
    def update_graph(selected_value: str):
        """
        Renders new graph for selected drop down menu value
        :param selected_value: selected column of dropdown menu
        :return: graph to be displayed
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        fig = Visualizator.render(
            df, selected_column=selected_value, original_id_col=original_id_col
        )
        return fig

    @app.callback(
        Output("ngl_molecule_viewer", "data"),
        Output("molecules_dropdown", "value"),
        Output("range_start", "disabled"),
        Output("range_end", "disabled"),
        Output("selected_atoms", "disabled"),
        Output("range_start", "options"),
        Output("range_end", "options"),
        Output("selected_atoms", "options"),
        Output("selected_atoms", "value"),
        Input("graph", "clickData"),
        Input("molecules_dropdown", "value"),
        Input("range_start", "value"),
        Input("range_end", "value"),
        Input("selected_atoms", "value"),
    )
    def display_molecule(
        clickdata,
        dd_molecules: list,
        range_start: int,
        range_end: int,
        selected_atoms: list,
    ):
        """
        callback function to handle the displaying of the molecule
        :param clickdata: given data by clicking on a datapoint in the 3D plot
        :param dd_molecules: selected molecules in the dropdown menu
        :param range_start: in the dropdown menu selected start
        :param range_end: in the dropdown menu selected end
        :param selected_atoms: selected values of the dropdown menu for highlighted atoms
        :return:
        """

        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # convert original to mapped IDs
        seq_ids = list()
        saved_seq_ids = list()
        if dd_molecules is not None:
            seq_ids = to_mapped_id(dd_molecules, original_id_col, df)

            # save unedited seq ids to replace edited seq ids later
            # (editing for range selection and highlighting)
            saved_seq_ids = list(seq_ids)

        # triggered by click on graph
        if ctx.triggered_id == "graph":
            # dict with data of clickdata
            points = clickdata["points"][0]
            # class_index value and it's index number
            index_num = int(points["pointNumber"])
            class_index = points["curveNumber"]

            # extract df_row of selected protein
            class_df = df[df["class_index"] == class_index]
            df_row = class_df.iloc[index_num]
            # add missing name to series
            name = pd.Series(class_df.index[index_num])
            name.index = ["Name"]
            df_row = pd.concat([name, df_row])

            # extract sequence ID
            seq_id = df_row["Name"]

            # add selected sequence ID to already selected IDs
            seq_ids.append(seq_id)

        # path to .pdb file
        struct_path = str(struct_container.get_structure_dir()) + "/"

        # disable range and highlighting selection in case more than 1 atom is selected
        start_disabled = False
        end_disabled = False
        atoms_disabled = False

        range_for_start = list()
        range_for_end = list()
        selectable_atoms = list()

        strand = None

        if len(seq_ids) > 1 or len(seq_ids) == 0:
            # disable the dropdown menus for range selection and highlighting if more than 1 molecule is selected
            start_disabled = True
            end_disabled = True
            atoms_disabled = True

        # only one molecule selected
        else:
            molecule_range, strand = struct_container.get_range(seq_ids[0])
            range_for_start = molecule_range
            range_for_end = molecule_range
            selectable_atoms = molecule_range

            # start of range selected
            if range_start is not None:
                # reset values
                range_for_end = []

                # remove numbers below for selection of end of range
                for num in molecule_range:
                    if num > range_start:
                        range_for_end.append(num)

                # set selectable atoms accordingly
                selectable_atoms = range_for_end

            # end of range selected
            if range_end is not None:
                # reset values
                range_for_start = []

                # remove numbers above for selection of start of range
                for num in molecule_range:
                    if num < range_end:
                        range_for_start.append(num)

                # set selectable atoms accordingly
                selectable_atoms = range_for_start

            # both, start and end selected
            if range_start is not None and range_end is not None:
                # reset values
                selectable_atoms = []

                for num in molecule_range:
                    if range_start <= num <= range_end:
                        selectable_atoms.append(num)

                # remove selected atom if not in range of selectable atoms
                if selected_atoms is not None:
                    for atom in selected_atoms:
                        if atom not in selectable_atoms:
                            selected_atoms.remove(atom)

            # seq id has to be edited if range start and end are selected
            seq_id = seq_ids[0]
            strand_set = False

            # append the range selection to the seq id accordingly
            if range_start is not None and range_end is not None:
                # append strand to seq id string
                seq_id = seq_id + f".{strand}"
                strand_set = True

                seq_id = seq_id + f":{range_start}-{range_end}"

            # append selected atoms to the string accordingly
            if selected_atoms is not None:
                if len(selected_atoms) > 0:
                    if not strand_set:
                        # append strand to seq id string
                        seq_id = seq_id + f".{strand}"

                    # bring selected atoms into string format comma separated
                    atoms = ""
                    for atom in selected_atoms:
                        atoms = atoms + f"{atom},"

                    # remove last comma
                    atoms = atoms[:-1]

                    seq_id = seq_id + f"@{atoms}"

            # replace seq id with new edited seq id
            seq_ids[0] = seq_id

        # data format for molecule viewer
        data_list = [
            ngl_parser.get_data(
                data_path=struct_path,
                pdb_id=seq_id,
                color="blue",
                reset_view=True,
                local=True,
            )
            for seq_id in seq_ids
        ]

        if ctx.triggered_id != "graph":
            # replace edited seq ids with saved unedited seq ids
            seq_ids = saved_seq_ids

        # back to original IDs
        seq_ids = to_original_id(seq_ids, original_id_col, df)

        # prevent updating if data list is empty since molecule viewer gives error otherwise
        if not data_list:
            raise PreventUpdate

        return (
            data_list,
            seq_ids,
            start_disabled,
            end_disabled,
            atoms_disabled,
            range_for_start,
            range_for_end,
            selectable_atoms,
            selected_atoms,
        )

    @app.callback(
        Output("ngl_molecule_viewer", "molStyles"),
        Input("representation_dropdown", "value"),
    )
    def set_mol_style(selected_representation):
        molstyles_dict = {
            "representations": selected_representation,
            "chosenAtomsColor": "white",
            "chosenAtomsRadius": 0.5,
            "molSpacingXaxis": 30,
            "sideByside": True,
        }

        return molstyles_dict


def get_callbacks(app, df):
    @app.callback(
        Output("graph", "figure"),
        Input("dd_menu", "value"),
    )
    def update_graph(selected_value: str, original_id_col: list):
        """
        Renders new graph for selected drop down menu value
        :param selected_value: selected column of dropdown menu
        :return: graph to be displayed
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        fig = Visualizator.render(
            df, selected_column=selected_value, original_id_col=original_id_col
        )
        return fig
