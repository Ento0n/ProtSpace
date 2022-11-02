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
        Input("graph", "clickData"),
        Input("molecules_dropdown", "value"),
        Input("range_start", "value"),
    )
    def display_molecule(clickdata, dd_molecules: list, range_start):
        """
        callback function to handle the displaying of the molecule
        :param clickdata: given data by clicking on a datapoint in the 3D plot
        :param dd_molecules: selected molecules in the dropdown menu
        :return:
        """

        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # convert original to mapped IDs
        seq_ids = list()
        if dd_molecules is not None:
            seq_ids = to_mapped_id(dd_molecules, original_id_col, df)

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

        # disable range and highlighting selection in case more than 1 atom is selected
        start_disabled = False
        end_disabled = False
        atoms_disabled = False

        range = list()
        if len(seq_ids) > 1 or len(seq_ids) == 0:
            start_disabled = True
            end_disabled = True
            atoms_disabled = True

        # only one molecule selected
        else:
            range = struct_container.get_range(seq_ids[0])

            # start of range selected
            if range_start is not None:
                # remove numbers below for selection of end of range
                range_for_end = list

                for num in range:
                    if num > range_start:
                        range_for_end.append(num)

        # back to original IDs
        seq_ids = to_original_id(seq_ids, original_id_col, df)

        print(f"data_list: {data_list}")

        # prevent updating if data list is empty since molecule viewer gives error otherwise
        if not data_list:
            raise PreventUpdate

        return (
            data_list,
            seq_ids,
            start_disabled,
            end_disabled,
            atoms_disabled,
            range,
            range,
            range,
        )

    @app.callback(
        Output("ngl_molecule_viewer", "molStyles"),
        Input("representation_dropdown", "value"),
    )
    def set_mol_style(selected_representation):
        molstyles_dict = {
            "representations": selected_representation,
            "chosenAtomsColor": "white",
            "chosenAtomsRadius": 1,
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
