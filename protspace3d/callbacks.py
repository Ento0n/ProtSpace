#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from visualization import Visualizator

from dash import Input, Output
from dash.exceptions import PreventUpdate
import dash
import dash_bio.utils.ngl_parser as ngl_parser
import pandas as pd


def get_callbacks(app, df, struct_container):
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

        fig = Visualizator.render(df, selected_column=selected_value)
        return fig

    @app.callback(
        Output("ngl_molecule_viewer", "data"),
        Output("ngl_molecule_viewer", "molStyles"),
        Input("graph", "clickData"),
    )
    def display_molecule(clickdata):
        """
        callback function to handle the displaying of the molecule
        :param clickdata: given data by clicking on a datapoint in the 3D plot
        :return:
        """

        if clickdata is None:
            raise PreventUpdate

        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

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

        # set structure container ID accordingly
        struct_container.set_structure_ids(seq_id)

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
        ]

        print(data_list)

        molstyles_dict = {
            "representations": ["cartoon", "axes+box"],
            "chosenAtomsColor": "white",
            "chosenAtomsRadius": 1,
            "molSpacingXaxis": 100,
        }

        return data_list, molstyles_dict
