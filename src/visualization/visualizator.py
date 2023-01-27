#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from colorsys import hls_to_rgb

import numpy as np
import plotly.graph_objects as go
from matplotlib import cm
from pandas import DataFrame

from .base import init_app
from .pdb import init_app_pdb


class Visualizator:
    SYMBOLS = [
        "circle",
        "square",
        "diamond",
        "cross",
        "circle-open",
        "diamond-open",
        "square-open",
    ]

    def __init__(self, fig: go.Figure, csv_header: list[str], dim_red: str):
        self.fig = fig
        self.csv_header = csv_header
        self.dim_red = dim_red

    @staticmethod
    def n_symbols_equation(n: int):
        """
        Output is the number of symbols used based on the number of groups
        :param n: number of groups
        :return: number of symbols
        """
        if n <= 8:
            return 1
        elif n <= 11:
            return 2
        elif n <= 15:
            return 3
        elif n <= 20:
            return 4
        elif n <= 26:
            return 5
        elif n <= 33:
            return 6
        else:
            return 7

    @staticmethod
    def n_samples_equation(n: int, val_range: int) -> float:
        """
        Output is the range in which saturation and luminositiy of the coloring are picked,
        dependent on the number of groups. The underlying equation is (x-10)/(x+10).
        :param n: number of groups
        :param val_range: maximal range
        :return: calculated range
        """
        numerator = n - 10
        denominator = n + 10
        res = numerator / denominator

        if res < 0:
            res = 0

        res_range = res * val_range

        return res_range

    @staticmethod
    def gen_distinct_colors(n: int, sort: bool = True):
        """
        Creates are list in rgb format that holds the most distinct colors for the number of groups.
        :param n: number of groups
        :param sort: Whether colors should be sorted or not.
        :return: list with colors
        """
        color_list = list()
        np.random.seed(42)
        hues = np.arange(0, 360, 360 / n)
        hues = hues[np.random.permutation(hues.size)]
        for hue in hues:
            # default saturation is 100 and range from 100-50
            standard_sat = 100
            sat_range = Visualizator.n_samples_equation(n, val_range=50)
            saturation = standard_sat - np.random.ranf() * sat_range
            # default luminosity is 50 and range from 35-65
            standard_lum = 50
            lum_range = Visualizator.n_samples_equation(n, val_range=30)
            luminosity = (standard_lum - lum_range / 2) + np.random.ranf() * lum_range
            # color list in hls style
            color_list.append(tuple([hue / 360, luminosity / 100, saturation / 100]))
        if sort:
            color_list.sort()

        # round small values, otherwise dash has difficulties to display
        rgb_list = []
        for h, l, s in color_list:
            # Also convert to rgb values
            rgb = hls_to_rgb(round(h, 2), round(l, 2), round(s, 2))

            # change value range from 0-1 to 0-255, otherwise saturation=100 is black
            rgb = list(rgb)
            for idx, value in enumerate(rgb):
                rgb[idx] = value * 255
            rgb = tuple(rgb)

            rgb_list.append(rgb)
        return rgb_list

    @staticmethod
    def handle_colorbar(
        col_groups: list, fig: go.Figure, n_symbols: int, color_list: list, two_d: bool
    ):
        """
        Creates a colorbar trace for the plot if only numeric values are in the group.
        :param col_groups: the column groups
        :param fig: the graph figure
        :param n_symbols: number of different symbols to be displayed.
        :param color_list: list with the different colors to be displayed.
        :param two_d: if True graph should be displayed in 2D
        :return: flag whether only numeric values are in the column and number of symbols to be used
        """
        numeric_flag = False
        if all(
            [
                isinstance(item, int) or isinstance(item, float) or item == "NA"
                for item in col_groups
            ]
        ):
            # Only numeric values and "NA" in the group
            numeric_flag = True

            # Only one symbol to be used
            n_symbols = 1

            # Find min and max value of the group, excluding "NA"
            min_val = sys.float_info.max
            max_val = sys.float_info.min
            for item in col_groups:
                if isinstance(item, int) or isinstance(item, float):
                    if min_val > item:
                        min_val = item
                    if max_val < item:
                        max_val = item

            if "NA" in col_groups:
                no_na_len = len(col_groups) - 1
            else:
                no_na_len = len(col_groups)

            viridis = cm.get_cmap("viridis", no_na_len)
            color_list = list()
            for i in range(len(col_groups)):
                # change rgba to rgb and range of values from 0-1 to 0-255
                rgba_tuple = viridis(i)
                rgba_list = list(rgba_tuple)
                red = rgba_list[0] * 255
                green = rgba_list[1] * 255
                blue = rgba_list[2] * 255
                rgb_list = [red, green, blue]
                rgb_tuple = tuple(rgb_list)

                color_list.append(rgb_tuple)

            # Use plotly Viridis as colorscale
            colorscale = list()
            for col in color_list:
                colorscale.append(f"rgb{col}")

            # create colorbar
            colorbar = dict(
                title="Colorbar",
                lenmode="fraction",
                len=0.5,
                yanchor="bottom",
                ypad=50,
            )

            # create and add a dummy trace that holds the colorbar
            if not two_d:
                color_trace = go.Scatter3d(
                    x=[None],
                    y=[None],
                    z=[None],
                    mode="markers",
                    marker=dict(
                        colorscale=colorscale,
                        showscale=True,
                        colorbar=colorbar,
                        cmin=min_val,
                        cmax=max_val,
                    ),
                    showlegend=False,
                )
            else:
                color_trace = go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        colorscale=colorscale,
                        showscale=True,
                        colorbar=colorbar,
                        cmin=min_val,
                        cmax=max_val,
                    ),
                    showlegend=False,
                )

            fig.add_trace(color_trace)

        return numeric_flag, n_symbols, color_list

    @staticmethod
    def customize_axis_titles(dim_red: str, fig: go.Figure, df: DataFrame):
        """
        Axis titles are edited dependent on dimensionality reduction (UMAP or PCA)
        :param dim_red: to be displayed dimensionality reduction
        :param fig: graph figure
        :param df: dataframe with all the data
        :return: None
        """
        if dim_red == "UMAP" or dim_red == "TSNE":
            fig.update_layout(
                # Remove axes ticks and labels as they are usually not informative
                scene=dict(
                    xaxis=dict(showticklabels=False, showspikes=False, title=""),
                    yaxis=dict(showticklabels=False, showspikes=False, title=""),
                    zaxis=dict(showticklabels=False, showspikes=False, title=""),
                ),
            )
        else:
            # extract variance column
            unique_variance_column = df["variance"].unique().tolist()

            pca_variance = list()
            for value in unique_variance_column:
                if value != "NA":
                    pca_variance.append(value)

            # Sort descending since the first component of PCA has more variance than the second and so on
            pca_variance.sort(reverse=True)

            fig.update_layout(
                # Remove axes ticks and labels as they are usually not informative
                scene=dict(
                    xaxis=dict(
                        showticklabels=False,
                        showspikes=False,
                        title=f"PC1 ({float(pca_variance[0]):.1f}%)",
                    ),
                    yaxis=dict(
                        showticklabels=False,
                        showspikes=False,
                        title=f"PC2 ({float(pca_variance[1]):.1f}%)",
                    ),
                    zaxis=dict(
                        showticklabels=False,
                        showspikes=False,
                        title=f"PC3 ({float(pca_variance[2]):.1f}%)",
                    ),
                ),
            )

    @staticmethod
    def handle_title(dim_red: str, umap_paras: dict, fig: go.Figure):
        """
        Sets up the title of the graph depending on the dimensionality reduction (UMAP or PCA)
        :param dim_red: to be displayed dimensionality reduction
        :param umap_paras: parameters of the UMAP calculation
        :param fig: graph figure
        :return: None
        """
        # Update title
        if dim_red == "UMAP":
            title = (
                "UMAP" + f"<br>n_neighbours: {umap_paras['n_neighbours']},"
                f" min_dist: {umap_paras['min_dist']}, metric: {umap_paras['metric']}"
            )
        elif dim_red == "PCA":
            title = "PCA"
        else:
            title = "TSNE"

        fig.update_layout(
            title={
                "text": title,
                "y": 0.98,
                "x": 0.4,
            }
        )

    @staticmethod
    def update_layout(fig: go.Figure):
        """
        Updates the layout of the graph figure. Margin, hoverinfo and legend positioning
        :param fig: graph Figure
        :return: None
        """
        # Set hoverinfo
        fig.update_traces(
            hoverinfo=["name", "text"],
            hoverlabel=dict(namelength=-1),
            hovertemplate="%{text}",
        )

        # Set legend in right upper corner of the plot
        fig.update_layout(legend=dict(yanchor="top", y=0.97, xanchor="right", x=0.99))

        # change margins of the graph
        fig.update_layout(margin=dict(l=1, r=1, t=1, b=1))

    @staticmethod
    # https://github.com/sacdallago/bio_embeddings/blob/develop/bio_embeddings/visualize/plotly_plots.py
    def render(
        df: DataFrame,
        selected_column: str,
        original_id_col: object,
        umap_paras: dict,
        dim_red: str = "UMAP",
        two_d: bool = False,
    ):
        """
        Renders the plotly graph with the selected column in the dataframe df
        :param df: dataframe
        :param selected_column: column of the dataframe
        :param original_id_col: the colum "original id" of the mapped csv file
        :param umap_paras: dictionary holding the parameters of UMAP
        :param dim_red: to be displayed dimensionality reduction
        :param two_d: if True plot should be in 2D
        :return: plotly graphical object
        """

        # custom separator to sort str, int and float (str case-insensitive)
        # order: 1. int and float 2. str 3. rest 4. NA
        def my_comparator(val):
            if isinstance(val, float) or isinstance(val, int):
                return 0, val
            elif val == "NA":
                return 3, val
            elif isinstance(val, str):
                val = val.lower()
                return 1, val
            else:
                return 2, val

        mapped_index = None
        if original_id_col is not None:
            # swap index
            mapped_index = df.index
            df.index = original_id_col

        col_groups = df[selected_column].unique().tolist()
        col_groups.sort(key=my_comparator)

        color_list = Visualizator.gen_distinct_colors(n=len(col_groups))

        # Figure out how many symbols to use depending on number of column groups
        n_symbols = Visualizator.n_symbols_equation(n=len(col_groups))

        fig = go.Figure()

        numeric_flag, n_symbols, color_list = Visualizator.handle_colorbar(
            col_groups, fig, n_symbols, color_list, two_d
        )

        df["class_index"] = np.ones(len(df)) * -100

        if dim_red == "UMAP":
            x = "x_umap"
            y = "y_umap"
            z = "z_umap"
        elif dim_red == "PCA":
            x = "x_pca"
            y = "y_pca"
            z = "z_pca"
        else:
            x = "x_tsne"
            y = "y_tsne"
            z = "z_tsne"

        # iterate over different values of the selected column
        for group_idx, group_value in enumerate(col_groups):
            # Show only NA in legend if colorbar is shown
            show_legend = True
            if numeric_flag:
                if group_value != "NA":
                    show_legend = False

            # set up color depending on numeric processing
            if not numeric_flag:
                color = f"rgb{color_list[group_idx]}"
            else:
                if group_value == "NA":
                    color = f"rgb(255,0,0)"
                else:
                    color = f"rgb{color_list[group_idx]}"

            # set up opacity dependent on NA or not
            opacity = 1.0
            if group_value == "NA":
                opacity = 0.4

            # extract df with only group value
            df_group = df[df[selected_column] == group_value]

            if not two_d:
                trace = go.Scatter3d(
                    x=df_group[x],
                    y=df_group[y],
                    z=df_group[z],
                    mode="markers",
                    name=group_value,
                    opacity=opacity,
                    marker=dict(
                        size=10,
                        color=color,
                        symbol=Visualizator.SYMBOLS[group_idx % n_symbols],
                        line=dict(color="black", width=1),
                    ),
                    text=df_group.index.to_list(),
                    showlegend=show_legend,
                )
            else:
                trace = go.Scatter(
                    x=df_group[x],
                    y=df_group[y],
                    mode="markers",
                    name=group_value,
                    opacity=opacity,
                    marker=dict(
                        size=10,
                        color=color,
                        symbol=Visualizator.SYMBOLS[group_idx % n_symbols],
                        line=dict(color="black", width=1),
                    ),
                    text=df_group.index.to_list(),
                    showlegend=show_legend,
                )
            fig.add_trace(trace)
            # Give the different group values a number
            df.loc[df[selected_column] == group_value, "class_index"] = group_idx

        Visualizator.update_layout(fig)

        Visualizator.handle_title(dim_red, umap_paras, fig)

        Visualizator.customize_axis_titles(dim_red, fig, df)

        # swap index again
        if original_id_col is not None:
            df.index = mapped_index

        return fig

    def get_base_app(self, umap_paras: dict):
        """
        Initializes the dash app in base.py
        :param umap_paras: Parameters of the UMAP calculation
        :return: the application layout
        """
        return init_app(umap_paras, self.csv_header, self.fig, self.dim_red)

    def get_pdb_app(self, orig_id_col: list[str], umap_paras: dict):
        """
        Initializes the dash app in pdb.py
        :param orig_id_col: List of the original IDs
        :param umap_paras: Parameters of the UMAP calculation
        :return: the application layout
        """
        return init_app_pdb(
            orig_id_col, umap_paras, self.csv_header, self.fig, self.dim_red
        )
