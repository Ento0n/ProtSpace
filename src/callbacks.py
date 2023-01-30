#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from visualization.visualizator import Visualizator
from structurecontainer import StructureContainer
from preprocessing import DataPreprocessor

from pathlib import Path
from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate
import dash
import dash_bio.utils.ngl_parser as ngl_parser
import plotly.graph_objects as go

import dash_bootstrap_components as dbc

from pandas import DataFrame
import numpy as np
import json
from statistics import mean


def to_mapped_id(sel_original_seq_ids: list, original_id_col: list, df: DataFrame):
    """
    Converts IDs from original to mapped
    :param sel_original_seq_ids: selected original sequence IDs
    :param original_id_col: list of original IDs
    :param df: Dataframe with all data
    :return: sequence IDs in mapped
    """
    seq_ids = list()

    for seq_id in sel_original_seq_ids:
        index_num = original_id_col.index(seq_id)
        mapped_id = df.index[index_num]
        seq_ids.append(mapped_id)

    return seq_ids


def to_original_id(sel_mapped_seq_ids: list, original_id_col: list, df: DataFrame):
    """
    Converts IDs from mapped to original
    :param sel_mapped_seq_ids: selected mapped sequence IDs
    :param original_id_col: list of original IDs
    :param df: Dataframe with all data
    :return: sequence IDs in mapped
    """
    seq_ids = list()

    for seq_id in sel_mapped_seq_ids:
        index_num = df.index.get_indexer_for([seq_id])[0]
        original_id = original_id_col[index_num]
        seq_ids.append(original_id)

    return seq_ids


def handle_highlighting(
    seq_ids: list,
    struct_container: StructureContainer,
    range_start: int,
    range_end: int,
    selected_atoms: list,
):
    """
    Takes selected sequence IDs and if one is selected and range, highlighting is set up, convert
    the ID into needed string representation.
    :param seq_ids: selected sequence IDs
    :param struct_container: structure container handling files
    :param range_start: selected range start
    :param range_end: selected range end
    :param selected_atoms: selected atoms to be highlighted
    :return: highlighting parameters
    """
    # disable range and highlighting selection in case more than 1 atom is selected
    start_disabled = False
    end_disabled = False
    atoms_disabled = False

    range_for_start = list()
    range_for_end = list()
    selectable_atoms = list()

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

    return (
        start_disabled,
        end_disabled,
        atoms_disabled,
        range_for_start,
        range_for_end,
        selectable_atoms,
    )


def clickdata_to_seqid(click_data: dict):
    """
    takes clickdata recieved from graph and converts it into the sequence ID
    :param click_data: graph click data
    :return: retrieved sequence ID
    """
    # dict with data of clickdata
    values = click_data["points"][0]
    seq_id = values["text"]

    return seq_id


def get_callbacks_pdb(
    app: dash.Dash,
    df: DataFrame,
    struct_container: StructureContainer,
    original_id_col: list,
):
    """
    Holds callbacks needed for pdb layout
    :param app: dash application
    :param df: dataframe with all data
    :param struct_container: structure container handling files
    :param original_id_col: list with original IDs
    :return: None
    """

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
        Output("download_molecule_button", "disabled"),
        Output("mol_name_storage", "data"),
        Output("clicked_mol_storage", "data"),
        Output("no_pdb_toast", "is_open"),
        Output("molecules_dropdown_save", "data"),
        Input("graph", "clickData"),
        Input("molecules_dropdown", "value"),
        Input("range_start", "value"),
        Input("range_end", "value"),
        Input("selected_atoms", "value"),
        Input("reset_view_button", "n_clicks"),
        Input("clicked_mol_storage", "data"),
        Input("molecules_dropdown_save", "data"),
    )
    def display_molecule(
        click_data: dict,
        dd_molecules: list,
        range_start: int,
        range_end: int,
        selected_atoms: list,
        reset_view_clicks: int,
        last_clicked_mol: str,
        dd_molecules_save: list,
    ):
        """
        callback function to handle the displaying of the molecule
        :param click_data: given data by clicking on a datapoint in the 3D plot
        :param dd_molecules: selected molecules in the dropdown menu
        :param range_start: in the dropdown menu selected start
        :param range_end: in the dropdown menu selected end
        :param selected_atoms: selected values of the dropdown menu for highlighted atoms
        :param reset_view_clicks: button to reset the view of the molecule viewer.
        :param last_clicked_mol: sequence ID of the last clicked molecule
        :param dd_molecules_save: molecule values in dd menu before new interaction
        :return: molecule viewer variables etc.
        """

        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Make redundant variable used
        if reset_view_clicks:
            pass

        # convert original to mapped IDs
        seq_ids = list()
        saved_seq_ids = list()
        if dd_molecules is not None:
            if original_id_col is not None:
                seq_ids = to_mapped_id(dd_molecules, original_id_col, df)
            else:
                seq_ids = dd_molecules

            # save unedited seq ids to replace edited seq ids later
            # (editing for range selection and highlighting)
            saved_seq_ids = list(seq_ids)

        # triggered by click on graph
        clicked_seq_id = last_clicked_mol
        if ctx.triggered_id == "graph":
            clicked_seq_id = clickdata_to_seqid(click_data)

            if original_id_col is not None:
                clicked_seq_id = to_mapped_id([clicked_seq_id], original_id_col, df)[0]

            # Add to seq ids or replace last clicked molecule
            if last_clicked_mol is None:
                seq_ids.append(clicked_seq_id)
            else:
                for idx, value in enumerate(seq_ids):
                    if value == last_clicked_mol:
                        seq_ids[idx] = clicked_seq_id

        # triggered by the molecule dropdown menu
        if ctx.triggered_id == "molecules_dropdown":
            # set clicked sequence id to None if the selection was deleted in the dropdown menu
            if clicked_seq_id not in seq_ids:
                clicked_seq_id = None

        # path to .pdb file
        struct_path = str(struct_container.get_structure_dir()) + "/"

        # needed for displaying molecules in dd menu
        orig_seq_ids = seq_ids
        if original_id_col is not None:
            orig_seq_ids = to_original_id(seq_ids, original_id_col, df)

        # check whether .pdb file is present
        for seq in seq_ids:
            if seq not in dd_molecules_save:
                file_path = Path(struct_path + seq + ".pdb")
                if not file_path.is_file():
                    return (
                        dash.no_update,
                        orig_seq_ids,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        dash.no_update,
                        clicked_seq_id,
                        True,
                        seq_ids,
                    )

        # handle the range and atom selection
        (
            start_disabled,
            end_disabled,
            atoms_disabled,
            range_for_start,
            range_for_end,
            selectable_atoms,
        ) = handle_highlighting(
            seq_ids, struct_container, range_start, range_end, selected_atoms
        )

        # data format for molecule viewer
        data_list = [
            ngl_parser.get_data(
                data_path=struct_path,
                pdb_id=seq_id,
                color="black",
                reset_view=True,
                local=True,
            )
            for seq_id in seq_ids
        ]

        if ctx.triggered_id != "graph":
            # replace edited seq ids with saved unedited seq ids
            seq_ids = saved_seq_ids

        mapped_seq_ids = seq_ids
        if original_id_col is not None:
            # back to original IDs
            seq_ids = to_original_id(seq_ids, original_id_col, df)

        # enable download button at first protein selection
        download_disabled = False

        # convert sequence ids list into string format
        mol_names = "_".join(seq_ids)

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
            download_disabled,
            mol_names,
            clicked_seq_id,
            dash.no_update,
            mapped_seq_ids,
        )

    @app.callback(
        Output("ngl_molecule_viewer", "molStyles"),
        Input("representation_dropdown", "value"),
        Input("spacing_slider", "value"),
    )
    def set_mol_style(selected_representation: list, spacing_slider_value: int):
        """
        Updates the representation of the molecule viewer
        :param selected_representation: in the dropdown menu selected representation(s)
        :param spacing_slider_value: value of the slider, space between molecules
        :return: filled dictionary
        """
        molstyles_dict = {
            "representations": selected_representation,
            "chosenAtomsColor": "white",
            "chosenAtomsRadius": 0.5,
            "molSpacingXaxis": spacing_slider_value,
            "sideByside": True,
        }

        return molstyles_dict

    @app.callback(
        Output("molecules_offcanvas", "is_open"),
        Input("molecules_settings_button", "n_clicks"),
    )
    def handle_molecules_canvas(button_click: int):
        """
        Opens the molecule viewer settings offcanvas
        :param button_click: settings button clicked
        :return: True for open offcanvas
        """
        if button_click:
            return True

    # molecule viewer sizing
    app.clientside_callback(
        """
        function(href) {
            var w = window.innerWidth;
            var h = window.innerHeight;
            return [h, w];
        }
        """,
        Output("molviewer_sizing_div", "children"),
        Input("url", "href"),
        Input("recal_size_button", "n_clicks"),
    )

    @app.callback(
        Output("ngl_molecule_viewer", "height"),
        Output("ngl_molecule_viewer", "width"),
        Output("height_slider", "value"),
        Output("width_slider", "value"),
        Output("moleculeviewer_div", "style"),
        Input("molviewer_sizing_div", "children"),
        Input("height_slider", "value"),
        Input("width_slider", "value"),
        Input("moleculeviewer_div", "style"),
        Input("distribution_slider", "value"),
    )
    def set_molviewer_size(
        sizing: list,
        slider_height: int,
        slider_width: int,
        div_style_dict: dict,
        distribution: int,
    ):
        """
        Calculates and or simply changes the molecule viewer height and width, also
        depending on space distribution.
        :param sizing: calculated height and width
        :param slider_height: height value of the height slider
        :param slider_width: width value of the width slider
        :param div_style_dict: Border of the molecule viewer
        :param distribution: space distribution slider value
        :return: height and width of the components
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # sliders are used
        if ctx.triggered_id == "height_slider" or ctx.triggered_id == "width_slider":
            # set style of div accordingly
            div_style_dict["height"] = str(slider_height + 1) + "px"
            div_style_dict["width"] = str(slider_width + 1) + "px"

            return (
                slider_height,
                slider_width,
                slider_height,
                slider_width,
                div_style_dict,
            )

        # automatic sizing
        height = (sizing[0] - 150) / 1.05
        width = sizing[1] / 2.1

        # fit sizing to distribution
        if distribution == 3:
            height = (sizing[0] - 150) / 1.05
            width = sizing[1] / 1.4
        elif distribution == 4:
            height = (sizing[0] - 150) / 1.05
            width = sizing[1] / 1.6
        elif distribution == 5:
            height = (sizing[0] - 150) / 1.05
            width = sizing[1] / 1.8
        elif distribution == 6:
            height = (sizing[0] - 150) / 1.05
            width = sizing[1] / 2.1
        elif distribution == 7:
            height = (sizing[0] - 150) / 1.05
            width = sizing[1] / 2.55
        elif distribution == 8:
            height = (sizing[0] - 180) / 1.05
            width = sizing[1] / 3.3
        elif distribution == 9:
            height = (sizing[0] - 180) / 1.05
            width = sizing[1] / 4.5

        # set style of div accordingly
        div_style_dict["height"] = str(height + 1) + "px"
        div_style_dict["width"] = str(width + 1) + "px"

        return height, width, height, width, div_style_dict

    @app.callback(
        Output("ngl_molecule_viewer", "downloadImage"),
        Output("ngl_molecule_viewer", "imageParameters"),
        Input("download_molecule_button", "n_clicks"),
        Input("filename_input", "value"),
        Input("mol_name_storage", "data"),
    )
    def download_molecule(button_clicks: int, filename_input: str, mol_names: str):
        """
        download option for molecule viewer, takes in selected sequence IDs in string format,
        a filename input and names the file accordingly
        :param button_clicks: button click of download button
        :param filename_input: text of the filename input
        :param mol_names: selected sequence IDs in string format
        :return: needed download parameters
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Make redundant variable used
        if button_clicks:
            pass

        if filename_input is None or filename_input == "":
            filename_input = mol_names

        image_parameters = {
            # Antialiasing, makes image smoother
            "antialias": True,
            # change the background from black to transparent
            "transparent": True,
            # trim background to edges of molecule, only works with transparent
            "trim": True,
            "defaultFilename": filename_input,
        }

        download_image = False

        if ctx.triggered_id == "download_molecule_button":
            download_image = True

        return download_image, image_parameters

    @app.callback(
        Output("left_col", "width"),
        Output("right_col", "width"),
        Input("distribution_slider", "value"),
    )
    def set_space_distribution(left_width: int):
        right_width = 12 - left_width

        return left_width, right_width


def get_callbacks(
    app,
    df: DataFrame,
    original_id_col: list,
    umap_paras: dict,
    output_d: Path,
    csv_header: list[str],
    embeddings: np.stack,
    embedding_uids: list,
    distance_dic: dict,
    umap_paras_dict: dict,
    fasta_dict: dict,
    struct_container: StructureContainer,
):
    """
    General callbacks needed for application
    :param app: application
    :param df: dataframe with all data
    :param original_id_col: list of original IDs
    :param umap_paras: UMAP parameters in dictionary
    :param output_d: output directory
    :param csv_header: the csv headers
    :param embeddings: the embeddings in a numpy stack
    :param embedding_uids: the unique IDs of the embeddings
    :param umap_paras_dict: already calculated UMAP parameters and their coordinates
    :param fasta_dict: fasta file in dictionary format
    :param struct_container: the structure container handling files
    :return:
    """

    @app.callback(
        Output("graph", "figure"),
        Output("n_neighbours_input", "disabled"),
        Output("min_dist_input", "disabled"),
        Output("metric_input", "disabled"),
        Output("last_umap_paras_dd", "options"),
        Output("last_umap_paras_dd", "disabled"),
        Output("umap_recalculation_button", "disabled"),
        Output("n_neighbours_input", "value"),
        Output("min_dist_input", "value"),
        Output("metric_input", "value"),
        Output("last_umap_paras_dd", "value"),
        Output("highlighting_bool", "data"),
        Output("relayoutData_save", "data"),
        Output("load_umap_spinner", "children"),
        Input("dd_menu", "value"),
        Input("dim_red_radio", "value"),
        Input("n_neighbours_input", "value"),
        Input("min_dist_input", "value"),
        Input("metric_input", "value"),
        Input("umap_recalculation_button", "n_clicks"),
        Input("last_umap_paras_dd", "value"),
        Input("graph", "clickData"),
        Input("highlighting_bool", "data"),
        Input("relayoutData_save", "data"),
        Input("dim_radio", "value"),
        State("graph", "figure"),
        State("graph", "relayoutData"),
    )
    def update_graph(
        selected_value: str,
        dim_red: str,
        n_neighbours: int,
        min_dist: float,
        metric: str,
        recal_button_clicks: int,
        umap_paras_dd_value: str,
        click_data: dict,
        highlighting_bool: bool,
        relayout_data_save: dict,
        dim: str,
        fig: go.Figure,
        relayout_data: dict,
    ):
        """
        Handles updating the graph when another group is selected, dimensionality reduction has changed,
        new parameters for UMAP should be calculated and applied, or a trace is removed/added for highlighting
        :param selected_value: selected group in dropdown menu
        :param dim_red: selected dimensionality reduction
        :param n_neighbours: UMAP value n_neighbours
        :param min_dist: UMAP value min_dist
        :param metric: UMAP value metric
        :param recal_button_clicks: Button for recalculating UMAP with new values and applying these.
        :param umap_paras_dd_value: selected UMAP parameters of already calculated ones
        :param click_data: data received from clicking the graph
        :param highlighting_bool: boolean indicating whether highlighting circle is already displayed or not
        :param relayout_data_save: relayout dict of the last click, needed for buggy plotly
        :param fig: graph Figure
        :param relayout_data: scene data of the graph
        :param dim: chosen dimension, 2D or 3D
        :return: Output variables
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Prevent constant resetting of the graph
        if ctx.triggered_id in ["n_neighbours_input", "min_dist_input", "metric_input"]:
            raise PreventUpdate

        # Make redundant variable used
        if recal_button_clicks:
            pass

        # In case dropdown menu is being cleared
        if selected_value is None:
            raise PreventUpdate

        # Prevent constant resetting of the graph
        if ctx.triggered_id in ["n_neighbours_input", "min_dist_input", "metric_input"]:
            raise PreventUpdate

        # load df into inner scope so that it can be modified
        nonlocal df

        # If plotly relayoutData is bugged and dict is empty for a reason, use last saved relayoutData
        if not relayout_data or "scene.camera" not in relayout_data:
            relayout_data = relayout_data_save

        # convert dictionary state of graph figure into go object
        fig = go.Figure(fig)

        umap_axis_names = ["x_umap", "y_umap", "z_umap"]

        # If umap parameters are selected in the dropdown menu
        if ctx.triggered_id == "last_umap_paras_dd":

            splits = umap_paras_dd_value.split(" ; ")
            umap_paras["n_neighbours"] = splits[0]
            umap_paras["min_dist"] = splits[1]
            umap_paras["metric"] = splits[2]

            coords = umap_paras_dict[umap_paras_dd_value]
            df.drop(labels=umap_axis_names, axis="columns", inplace=True)
            df["x_umap"] = coords["x_umap"]
            df["y_umap"] = coords["y_umap"]
            df["z_umap"] = coords["z_umap"]

        # If UMAP parameters are changed and accepted
        if ctx.triggered_id == "umap_recalculation_button":
            umap_paras["n_neighbours"] = n_neighbours
            umap_paras["min_dist"] = min_dist
            umap_paras["metric"] = metric

            # String representation of the current UMAP parameters
            umap_paras_string = (
                str(n_neighbours) + " ; " + str(min_dist) + " ; " + metric
            )

            df.drop(labels=umap_axis_names, axis="columns", inplace=True)

            df_umap = DataPreprocessor.generate_umap(embeddings, umap_paras)
            df_umap.index = embedding_uids

            df = df.join(df_umap, how="outer")

            coords_dict = dict(
                x_umap=df_umap["x_umap"].to_list(),
                y_umap=df_umap["y_umap"].to_list(),
                z_umap=df_umap["z_umap"].to_list(),
            )

            umap_paras_dict[umap_paras_string] = coords_dict
        # String representation of UMAP parameters still to be created if not button used
        else:
            umap_paras_string = (
                str(umap_paras["n_neighbours"])
                + " ; "
                + str(umap_paras["min_dist"])
                + " ; "
                + umap_paras["metric"]
            )

        if dim == "2D":
            two_d = True
        else:
            two_d = False

        if (
            ctx.triggered_id == "dd_menu"
            or ctx.triggered_id == "umap_recalculation_button"
            or ctx.triggered_id == "dim_red_radio"
            or ctx.triggered_id == "last_umap_paras_dd"
            or ctx.triggered_id == "dim_radio"
        ):
            fig = Visualizator.render(
                df,
                selected_column=selected_value,
                original_id_col=original_id_col,
                dim_red=dim_red,
                umap_paras=umap_paras,
                two_d=two_d,
            )
            # Set highlighting_bool to False since new graph is displayed and highlighting circle is removed
            highlighting_bool = False

        # Add trace that highlights the selected molecule with a circle
        # get seq id from click data
        if ctx.triggered_id == "graph":
            seq_id = clickdata_to_seqid(click_data)

            if original_id_col is not None:
                seq_id = to_mapped_id([seq_id], original_id_col, df)[0]

            if dim_red == "UMAP":
                x = df.at[seq_id, "x_umap"]
                y = df.at[seq_id, "y_umap"]
                z = df.at[seq_id, "z_umap"]
            elif dim_red == "PCA":
                x = df.at[seq_id, "x_pca"]
                y = df.at[seq_id, "y_pca"]
                z = df.at[seq_id, "z_pca"]
            else:
                x = df.at[seq_id, "x_tsne"]
                y = df.at[seq_id, "y_tsne"]
                z = df.at[seq_id, "z_tsne"]

            # Remove last trace (highlighting circle) if present
            if highlighting_bool:
                data_list = list(fig.data)
                data_list.pop(-1)
                fig.data = data_list

            if dim == "3D":
                fig.add_trace(
                    go.Scatter3d(
                        x=[x],
                        y=[y],
                        z=[z],
                        mode="markers",
                        marker=dict(
                            size=15,
                            color="black",
                            symbol="circle-open",
                        ),
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )
            else:
                fig.add_trace(
                    go.Scatter(
                        x=[x],
                        y=[y],
                        mode="markers",
                        marker=dict(
                            size=15,
                            color="black",
                            symbol="circle-open",
                        ),
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )

            # set highlighting_bool to true since highlighting circle is now displayed
            highlighting_bool = True

            # set camera to old settings so that camera stays in its position and doesn't reset
            if relayout_data:
                fig.update_layout(scene_camera=relayout_data["scene.camera"])

            # save relayoutData for buggy plotly
            relayout_data_save = relayout_data

        # Disable UMAP parameter input or not?
        disabled = False
        if not dim_red == "UMAP":
            disabled = True

        return (
            fig,
            disabled,
            disabled,
            disabled,
            list(umap_paras_dict.keys()),
            disabled,
            disabled,
            umap_paras["n_neighbours"],
            umap_paras["min_dist"],
            umap_paras["metric"],
            umap_paras_string,
            highlighting_bool,
            relayout_data_save,
            "Recalculate UMAP",
        )

    @app.callback(
        Output("disclaimer_modal", "is_open"),
        Input("disclaimer_modal_button", "n_clicks"),
    )
    def close_disclaimer_modal(button: int):
        """
        Handles the closing of the disclaimer modal on button click
        :param button: "Agree" button on disclaimer modal
        :return: open state of disclaimer modal
        """
        if button:
            return False

        return True

    @app.callback(
        Output("help_modal", "is_open"),
        Input("help_button", "n_clicks"),
    )
    def open_help_modal(button: int):
        """
        handles the closing of the help modal on button click
        :param button: "Close" button on help modal
        :return: open state of help modal
        """
        if button:
            return True

        return False

    @app.callback(
        Output("graph_offcanvas", "is_open"), Input("graph_settings_button", "n_clicks")
    )
    def handle_graph_canvas(button_click: int):
        """
        Opens settings offcanvas for graph on button click
        :param button_click: settings button click in graph container
        :return: open settings
        """
        if button_click:
            return True

    @app.callback(
        Output("html_download_toast", "is_open"),
        Input("dd_menu", "value"),
        Input("html_download_button", "n_clicks"),
        Input("button_html_all", "n_clicks"),
        Input("dim_red_radio", "value"),
    )
    def create_html(dd_value: str, button: int, all_button: int, dim_red: str):
        """
        Creates html of selected group on button click and indicates this with an download toast
        :param dd_value: selected group in the dropdown menu
        :param button: html download button
        :param all_button: button indicating all groups should be downloaded
        :param dim_red: selected dimensionality reduction
        :return: open download toast
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Make redundant variable used
        if button or all_button:
            pass

        if ctx.triggered_id == "html_download_button":
            fig = Visualizator.render(
                df, dd_value, original_id_col, umap_paras, dim_red
            )
            fig.write_html(output_d / f"3Dspace_{dd_value}_{dim_red}.html")

            return True

        if ctx.triggered_id == "button_html_all":
            for header in csv_header:
                fig = Visualizator.render(
                    df, header, original_id_col, umap_paras, dim_red
                )
                fig.write_html(output_d / f"3Dspace_{header}_{dim_red}.html")

            return True

    def get_expand_sequence_tooltip(button_id: str):
        """
        Returns the tooltip for the expand sequence button
        :param button_id: target button id the tooltip is to be displayed
        :return: tooltip
        """
        tooltip = dbc.Tooltip(
            "Expand sequence",
            target=button_id,
            placement="right",
        )
        return tooltip

    def get_collapse_sequence_tooltip(button_id: str):
        """
        Returns the tooltip for the collapse sequence button
        :param button_id: target button id the tooltip is to be displayed
        :return: tooltip
        """
        tooltip = dbc.Tooltip(
            "Collapse sequence",
            target=button_id,
            placement="top",
        )
        return tooltip

    def get_json_info(seq_id: str, info_text: list):
        """
        Creates json info in needed format for the info toast
        :param seq_id: selected sequence ID
        :param info_text: info text list to be filled
        :return: info_text list filled or not
        """
        if struct_container.json_flag:
            json_file = struct_container.get_json_dir() / f"{seq_id}.json"
            if json_file.is_file():
                with open(json_file) as f:
                    json_dict = json.load(f)

                ptm = json_dict["ptm"]
                plddt = mean(json_dict["plddt"])

                info_text.append(
                    dbc.ListGroupItem(
                        [
                            html.B("plDDT mean: "),
                            html.P(f"{round(plddt, 2)}"),
                            html.B("pTM: "),
                            html.P(f"{ptm}"),
                        ]
                    )
                )

    def get_fasta_info(seq_id: str, info_text: list):
        """
        Creates fasta info in needed format for the info toast
        :param seq_id: selected sequence ID
        :param info_text: info text list to be filled
        :return: info_text list filled or not
        """
        if fasta_dict is not None:
            if seq_id in fasta_dict.keys():
                sequence = str(fasta_dict[seq_id].seq)
                info_text.append(
                    dbc.ListGroupItem(
                        [
                            get_expand_sequence_tooltip(button_id="expand_seq_button"),
                            get_collapse_sequence_tooltip(
                                button_id="collapse_seq_button"
                            ),
                            html.B("Sequence:"),
                            html.Div(
                                id="collapsed_seq_div",
                                hidden=False,
                                children=[
                                    dbc.Row(
                                        [
                                            dbc.Col(
                                                [
                                                    html.P(sequence[:5]),
                                                ],
                                                width=6,
                                            ),
                                            dbc.Col(
                                                [
                                                    html.Button(
                                                        "...",
                                                        id="expand_seq_button",
                                                        style={
                                                            "padding": 0,
                                                            "border": "none",
                                                            "background": "none",
                                                        },
                                                    )
                                                ],
                                                width=6,
                                            ),
                                        ]
                                    ),
                                ],
                            ),
                            html.Div(
                                id="expanded_seq_div",
                                hidden=True,
                                children=[
                                    html.Button(
                                        "...",
                                        id="collapse_seq_button",
                                        style={
                                            "padding": 0,
                                            "border": "none",
                                            "background": "none",
                                        },
                                    ),
                                    html.P(f"{sequence}"),
                                ],
                            ),
                            html.B("Seq. length:"),
                            html.P(f"{len(sequence)}"),
                        ]
                    )
                )

    def get_group_info(seq_id: str, info_text: list):
        """
        Creates group info in needed format for the info toast
        :param seq_id: selected sequence ID
        :param info_text: info text list to be filled
        :return: info_text list filled or not
        """
        # group info in children format
        group_info_children = [
            dbc.Button(
                children=[
                    "Group info",
                    html.I(className="bi bi-arrow-up"),
                ],
                id="group_info_collapse_button",
                color="dark",
                outline=True,
                style={"border": "none"},
            ),
        ]
        for header in csv_header:
            group_info_children.append(html.B(f"{header}:"))
            group_info_children.append(html.P(f"{df.at[seq_id, header]}"))

        info_text.append(
            dbc.ListGroupItem(
                [
                    html.Div(
                        id="group_info_collapsed_div",
                        children=[
                            dbc.Button(
                                children=[
                                    "Group info",
                                    html.I(className="bi bi-arrow-down"),
                                ],
                                id="group_info_expand_button",
                                color="dark",
                                outline=True,
                                style={
                                    "border": "none",
                                },
                            ),
                        ],
                        hidden=False,
                    ),
                    html.Div(
                        id="group_info_expanded_div",
                        children=group_info_children,
                        hidden=True,
                    ),
                ]
            )
        )

        # Don't show toast if no information is present
        if len(info_text) == 0:
            raise PreventUpdate

    @app.callback(
        Output("info_toast", "header"),
        Output("info_toast", "children"),
        Output("info_toast", "is_open"),
        Input("graph", "clickData"),
    )
    def show_info_toast(click_data: dict):
        """
        Opens info toast and fills it with information for selected molecule in graph
        :param click_data: data received from clicking on graph
        :return: open info toast and give text to display
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # get seq id from click data
        actual_seq_id = clickdata_to_seqid(click_data)

        seq_id = actual_seq_id
        if original_id_col is not None:
            seq_id = to_mapped_id([seq_id], original_id_col, df)[0]

        info_header = actual_seq_id

        info_text = []

        get_json_info(seq_id, info_text)

        get_fasta_info(seq_id, info_text)

        get_group_info(seq_id, info_text)

        info_text = dbc.ListGroup(info_text, flush=True)

        open_now = True

        return info_header, info_text, open_now

    @app.callback(
        Output("collapsed_seq_div", "hidden"),
        Output("expanded_seq_div", "hidden"),
        Input("expand_seq_button", "n_clicks"),
        Input("collapse_seq_button", "n_clicks"),
    )
    def expand_sequence(expand_button: int, collapse_button: int):
        """
        Expands or collapses the sequence in the info toast
        :param expand_button: button indicating sequence should expand
        :param collapse_button: button indicating sequence should be collapsed
        :return: which sequence to display
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Make redundant variable used
        if expand_button or collapse_button:
            pass

        collapse_hidden = False
        expand_hidden = True
        if ctx.triggered_id == "expand_seq_button":
            collapse_hidden = True
            expand_hidden = False

        return collapse_hidden, expand_hidden

    @app.callback(
        Output("group_info_expanded_div", "hidden"),
        Output("group_info_collapsed_div", "hidden"),
        Input("group_info_expand_button", "n_clicks"),
        Input("group_info_collapse_button", "n_clicks"),
    )
    def handle_group_info(expand_button: int, collapse_button: int):
        """
        Expand or collapse group info on button click
        :param expand_button: button indicating group info should be shown
        :param collapse_button: button indicating group info should be hidden.
        :return: Hide or show group info
        """
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        # Make redundant variable used
        if expand_button or collapse_button:
            pass

        collapse_hidden = False
        expand_hidden = True

        if ctx.triggered_id == "group_info_expand_button":
            collapse_hidden = True
            expand_hidden = False

        return expand_hidden, collapse_hidden

    @app.callback(
        Output("neighbour_toast", "header"),
        Output("neighbour_toast", "children"),
        Output("neighbour_toast", "is_open"),
        Input("graph", "clickData"),
    )
    def show_neighbour_toast(clickData):
        # Check whether an input is triggered
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        return "Test", "Test", True
