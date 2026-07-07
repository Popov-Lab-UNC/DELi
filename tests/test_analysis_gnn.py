"""Tests for GNN analysis helpers and classifier orchestration."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest
import torch

torch = pytest.importorskip("torch")
pytest.importorskip("torch_geometric")

from deli.analysis.analyzers.gnn_analyzer import GNNAnalyzer
from deli.analysis.cube_class import DELi_Cube
from deli.analysis.gnn import Final_Network, MyConv, MyDataset, smi_to_pyg

pytestmark = pytest.mark.ml

ANALYSIS_INDEXES = {"sel": ["corrected_idx1", "corrected_idx2", "corrected_idx3"]}
ANALYSIS_RAW_INDEXES = {"sel": ["raw_idx1", "raw_idx2", "raw_idx3"]}
ANALYSIS_CONTROL_COLS = {"sel": ["control_sel"], "ntc": ["control_sel"]}


@pytest.fixture
def gnn_cube(ml_cube_df: pd.DataFrame) -> DELi_Cube:
    return DELi_Cube(
        ml_cube_df,
        "DEL_ID",
        ANALYSIS_INDEXES,
        control_cols=ANALYSIS_CONTROL_COLS,
        lib_size=1000,
        raw_indexes=ANALYSIS_RAW_INDEXES,
    )


def test_myconv_supports_gin_and_gat():
    gin = MyConv(node_dim=8, edge_dim=4, arch="GIN")
    gat = MyConv(node_dim=8, edge_dim=4, arch="GAT")
    assert gin is not None
    assert gat is not None


def test_myconv_rejects_unknown_architecture():
    with pytest.raises(ValueError, match="Unknown architecture"):
        MyConv(node_dim=8, edge_dim=4, arch="GCN")


def test_smi_to_pyg_builds_graph_for_valid_smiles():
    graph = smi_to_pyg("CCO", y=1, encoding="embedding")
    assert graph is not None
    assert graph.x.ndim == 2
    assert graph.edge_index.ndim == 2
    assert graph.y.item() == 1


def test_smi_to_pyg_returns_none_for_invalid_smiles():
    assert smi_to_pyg("not_a_smiles", y=0) is None


def test_my_dataset_filters_invalid_graphs():
    dataset = MyDataset(["CCO", "not_a_smiles", "CCC"], [1, 0, 1], encoding="embedding")
    assert len(dataset) == 2


def test_final_network_forward_pass():
    graph = smi_to_pyg("CCO", y=1, encoding="embedding")
    model = Final_Network(
        node_dim=graph.x.shape[1],
        edge_dim=graph.edge_attr.shape[1],
        arch="GAT",
        num_layers=2,
        encoding="embedding",
    )
    batch = torch.zeros(graph.num_nodes, dtype=torch.long)
    model.eval()
    out = model(graph.x, graph.edge_index, graph.edge_attr, batch)
    assert out.shape == (1, 2)


def test_gnn_analyzer_prepare_experiment_data(gnn_cube: DELi_Cube):
    analyzer = GNNAnalyzer(gnn_cube.data, ANALYSIS_INDEXES)
    subset = analyzer._prepare_experiment_data("sel", ANALYSIS_INDEXES["sel"], threshold=2)
    assert "graphs" in subset.columns
    assert "target" in subset.columns
    assert subset["graphs"].notna().all()
    assert len(subset) == 300


def test_gnn_analyzer_create_pytorch_datasets(gnn_cube: DELi_Cube):
    analyzer = GNNAnalyzer(gnn_cube.data, ANALYSIS_INDEXES)
    subset = analyzer._prepare_experiment_data("sel", ANALYSIS_INDEXES["sel"], threshold=2)
    loaders = analyzer._create_pytorch_datasets(subset)
    train_loader, val_loader, test_loader, node_dim, edge_dim = loaders
    assert len(train_loader.dataset) > 0
    assert len(test_loader.dataset) > 0
    assert node_dim > 0
    assert edge_dim > 0


def test_gnn_analyzer_plot_confusion_matrix(tmp_path: Path):
    GNNAnalyzer._plot_confusion_matrix(
        y_true=[0, 1, 1, 0],
        y_pred=[0, 1, 0, 0],
        acc=0.75,
        exp_name="sel",
        arch="GAT",
        threshold=2,
        output_dir=str(tmp_path),
    )
    assert list(tmp_path.glob("sel_GAT_classifier.svg"))


def test_gnn_analyzer_classify_with_mocked_training(gnn_cube: DELi_Cube, tmp_path: Path):
    analyzer = GNNAnalyzer(gnn_cube.data, ANALYSIS_INDEXES)
    subset = analyzer._prepare_experiment_data("sel", ANALYSIS_INDEXES["sel"], threshold=2)
    _, _, test_loader, node_dim, edge_dim = analyzer._create_pytorch_datasets(subset)
    device = torch.device("cpu")
    model = Final_Network(node_dim, edge_dim, arch="GAT", num_layers=1, encoding="embedding").to(
        device
    )

    with (
        patch.object(GNNAnalyzer, "_train_model", return_value=model),
        patch.object(GNNAnalyzer, "_evaluate_model", return_value=([0, 1], [0, 1], 1.0)),
        patch.object(GNNAnalyzer, "_plot_confusion_matrix") as plot_mock,
    ):
        analyzer.classify(threshold=2, output_dir=str(tmp_path), arch="GAT", num_layers=1)

    plot_mock.assert_called_once()


def test_deli_cube_gnn_classifier_delegates_to_analyzer(gnn_cube: DELi_Cube, tmp_path: Path):
    with patch("deli.analysis.analyzers.gnn_analyzer.GNNAnalyzer.classify") as classify_mock:
        gnn_cube.gnn_classifier(threshold=2, output_dir=str(tmp_path), arch="GAT", num_layers=1)
    classify_mock.assert_called_once()
