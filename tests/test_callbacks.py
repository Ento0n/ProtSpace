from pathlib import Path

from src.callbacks import get_callbacks, get_callbacks_pdb
from src.preprocessing import  DataPreprocessor
from src.structurecontainer import StructureContainer
from src.visualization.visualizator import Visualizator

from contextvars import copy_context
from dash._callback_context import context_value
from dash._utils import AttributeDict


def setup():
    """
    Handles the process of the application
    :return: app & html_flag
    """

    output_d = Path("data/VA")
    hdf_path = Path("data/VA/VA.h5")
    csv_path = Path("data/VA/VA.csv")
    fasta_path = Path("data/VA/VA.fasta")
    csv_sep = ","
    uid_col = 0
    html_cols = None
    pdb_d = None
    json_d = None
    reset = False
    conf = None
    pca_flag = False
    tsne_flag = False
    iterations = 1000
    perplexity = 30
    learning_rate = 10
    tsne_metric = "euclidean"
    n_neighbours = 25
    min_dist = 0.5
    metric = "euclidean"
    port = 8050
    verbose = True

    dim_red = "UMAP"
    if pca_flag:
        dim_red = "PCA"
    elif tsne_flag:
        dim_red = "TSNE"

    # put UMAP parameters in dictionary
    umap_paras = dict()
    umap_paras["n_neighbours"] = n_neighbours
    umap_paras["min_dist"] = min_dist
    umap_paras["metric"] = metric

    # Put TSNE parameters in dictionary
    tsne_paras = dict()
    tsne_paras["iterations"] = iterations
    tsne_paras["perplexity"] = perplexity
    tsne_paras["learning_rate"] = learning_rate
    tsne_paras["tsne_metric"] = tsne_metric

    # Create data preprocessor object
    data_preprocessor = DataPreprocessor(
        output_d,
        hdf_path,
        csv_path,
        fasta_path,
        csv_sep,
        uid_col,
        html_cols,
        reset,
        dim_red,
        umap_paras,
        tsne_paras,
        verbose,
    )

    # Preprocessing
    (
        df,
        fig,
        csv_header,
        original_id_col,
        embeddings,
        embedding_uids,
        distance_dic,
        fasta_dict,
    ) = data_preprocessor.data_preprocessing()

    # initialize structure container if flag set
    structure_container = StructureContainer(pdb_d, json_d)

    # Create visualization object
    visualizator = Visualizator(fig, csv_header, dim_red)

    # get ids of the proteins
    if original_id_col is not None:
        ids = original_id_col
    else:
        ids = df.index.to_list()

    umap_paras_dict = data_preprocessor.get_umap_paras_dict(df)
    tsne_paras_dict = data_preprocessor.get_tsne_paras_dict(df)

    # --- APP creation ---
    if structure_container.pdb_flag:
        application = visualizator.get_pdb_app(ids, umap_paras, tsne_paras)
    else:
        application = visualizator.get_base_app(umap_paras, tsne_paras, ids)

    download_graph, expand_sequence, handle_graph_canvas = get_callbacks(application, df, original_id_col, umap_paras, tsne_paras, output_d, csv_header, embeddings, embedding_uids, distance_dic, umap_paras_dict, tsne_paras_dict, fasta_dict, structure_container)

    return download_graph, expand_sequence, handle_graph_canvas


def test_download_graph():
    def run_callback():
        context_value.set(AttributeDict(**{"triggered_inputs": [{"prop_id": "graph_download_button.n_clicks"}, {"prop_id": "button_graph_all.n_clicks"}]}))
        return download_graph("Assigned group", 1, 1, "UMAP", "2D")

    download_graph, two, three = setup()
    ctx = copy_context()
    output = ctx.run(run_callback)
    assert output is True


def test_expand_sequence():
    def run_callback_1():
        context_value.set(AttributeDict(**{"triggered_inputs": [{"prop_id": "expand_seq_button.n_clicks"}]}))
        return expand_sequence(1, 1)

    one, expand_sequence, three = setup()
    ctx = copy_context()
    output = ctx.run(run_callback_1)
    assert output == (True, False)


def test_collapse_sequence():
    def run_callback():
        context_value.set(AttributeDict(**{"triggered_inputs": [{"prop_id": "dummy_prop_id.n_clicks"}]}))
        return expand_sequence(1, 1)

    one, expand_sequence, three = setup()
    ctx = copy_context()
    output = ctx.run(run_callback)
    assert output == (False, True)


def test_open_graph_settings():
    one, two, handle_graph_canvas = setup()
    output = handle_graph_canvas(1)
    assert output is True

