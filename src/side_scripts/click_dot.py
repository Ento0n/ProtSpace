import json
from textwrap import dedent as d

import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import numpy as np

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

styles = {"pre": {"border": "thin lightgrey solid", "overflowX": "scroll"}}


fig = go.Figure(
    go.Scatter3d(
        x=np.arange(10),
        y=np.arange(10),
        z=np.arange(10),
        marker_size=8 * np.ones(10),
        mode="markers",
    )
)

app.layout = html.Div(
    [
        dcc.Graph(id="basic-interactions", figure=fig),
        html.Div(
            className="row",
            children=[
                html.Div(
                    [
                        dcc.Markdown(
                            d(
                                """
                **Hover Data**

                Mouse over values in the graph.
            """
                            )
                        ),
                        html.Pre(id="hover-data", style=styles["pre"]),
                    ],
                    className="three columns",
                ),
                html.Div(
                    [
                        dcc.Markdown(
                            d(
                                """
                **Click Data**

                Click on points in the graph.
            """
                            )
                        ),
                        html.Pre(id="click-data", style=styles["pre"]),
                    ],
                    className="three columns",
                ),
                html.Div(
                    [
                        dcc.Markdown(
                            d(
                                """
                **Zoom and Relayout Data**

                Click and drag on the graph to zoom or click on the zoom
                buttons in the graph's menu bar.
                Clicking on legend items will also fire
                this event.
            """
                            )
                        ),
                        html.Pre(id="relayout-data", style=styles["pre"]),
                    ],
                    className="three columns",
                ),
            ],
        ),
    ]
)


@app.callback(
    Output("hover-data", "children"), [Input("basic-interactions", "hoverData")]
)
def display_hover_data(hoverData):
    return json.dumps(hoverData, indent=2)


@app.callback(
    [Output("click-data", "children"), Output("basic-interactions", "figure")],
    [Input("basic-interactions", "clickData")],
    [State("basic-interactions", "relayoutData")],
)
def display_click_data(clickData, relayoutData):
    if not clickData:
        return dash.no_update, dash.no_update
    point = clickData["points"][0]
    # Do something only for a specific trace
    if point["curveNumber"] > 0:
        return dash.no_update, dash.no_update
    else:
        fig = go.Figure(
            go.Scatter3d(
                x=np.arange(10), y=np.arange(10), z=np.arange(10), mode="markers"
            )
        )
        sizes = 8 * np.ones(10)
        sizes[point["pointNumber"]] = 15
        colors = [
            "blue",
        ] * 10
        colors[point["pointNumber"]] = "red"
        fig.update_traces(marker_size=sizes, marker_color=colors)
    # Make sure the view/angle stays the same when updating the figure
    if relayoutData and "scene.camera" in relayoutData:
        fig.update_layout(scene_camera=relayoutData["scene.camera"])
    return json.dumps(clickData, indent=2), fig


@app.callback(
    Output("relayout-data", "children"), [Input("basic-interactions", "relayoutData")]
)
def display_relayout_data(relayoutData):
    return json.dumps(relayoutData, indent=2)


if __name__ == "__main__":
    app.run_server(debug=True, port=8089)
