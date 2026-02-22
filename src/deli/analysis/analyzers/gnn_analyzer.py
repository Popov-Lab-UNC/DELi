import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import train_test_split
from torch_geometric.loader import DataLoader
from tqdm import tqdm

from deli.analysis.gnn import Final_Network, smi_to_pyg


class GNNAnalyzer:
    """Graph Neural Network analysis for DEL data"""

    def __init__(self, data, indexes):
        """
        Initialize the GNN analyzer.

        Parameters
        ----------
        data : pd.DataFrame
            DataFrame with DEL data including SMILES column
        indexes : Dict[str, list]
            Dictionary mapping experiment names to column lists
        """
        self.data = data
        self.indexes = indexes

    def classify(
        self,
        threshold: int = 2,
        output_dir: str = ".",
        arch: str = "GAT",
        num_layers: int = 3,
        encoding: str = "embedding",
    ) -> None:
        """
        Train a Graph Neural Network classifier to predict enrichment status.

        Parameters
        ----------
        threshold : int, optional
            Threshold value for enrichment status, by default 2
        output_dir : str, optional
            Directory to save the classifier image, by default "."
        arch : str, optional
            GNN architecture (GAT or GCN), by default "GAT"
        num_layers : int, optional
            Number of GNN layers, by default 3
        encoding : str, optional
            Node encoding method (embedding or onehot), by default "embedding"
        """
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        for exp_name, indices in self.indexes.items():
            subset_data = self._prepare_experiment_data(exp_name, indices, threshold)

            train_loader, val_loader, test_loader, node_dim, edge_dim = (
                self._create_pytorch_datasets(subset_data)
            )

            model = self._train_model(
                train_loader=train_loader,
                val_loader=val_loader,
                node_dim=node_dim,
                edge_dim=edge_dim,
                arch=arch,
                num_layers=num_layers,
                encoding=encoding,
                output_dir=output_dir,
                exp_name=exp_name,
                device=device,
            )

            y_true, y_pred, acc = self._evaluate_model(
                model=model, test_loader=test_loader, device=device
            )

            self._plot_confusion_matrix(
                y_true=y_true,
                y_pred=y_pred,
                acc=acc,
                exp_name=exp_name,
                arch=arch,
                threshold=threshold,
                output_dir=output_dir,
            )

    def _prepare_experiment_data(
        self, exp_name: str, indices: list, threshold: int
    ) -> pd.DataFrame:
        """
        Prepare data for a single experiment.

        Parameters
        ----------
        exp_name : str
            Name of the experiment
        indices : list
            List of column names for this experiment
        threshold : int
            Threshold value for enrichment status

        Returns
        -------
        pd.DataFrame
            Processed subset data with target labels
        """
        self.data[f"{exp_name}_average_enrichment"] = self.data[indices].mean(axis=1)

        top_100 = self.data.nlargest(100, f"{exp_name}_average_enrichment")
        remaining_data = self.data.drop(top_100.index)
        random_200 = remaining_data.sample(n=200, random_state=42)

        subset_data = pd.concat([top_100, random_200]).copy()
        subset_data["graphs"] = [
            smi_to_pyg(smiles, 0, encoding="embedding") for smiles in tqdm(subset_data["SMILES"])
        ]
        subset_data.dropna(subset=["graphs"], inplace=True)
        subset_data["target"] = (subset_data[f"{exp_name}_average_enrichment"] > threshold).astype(
            int
        )
        return subset_data

    @staticmethod
    def _create_pytorch_datasets(subset_data: pd.DataFrame) -> tuple:
        """
        Convert processed data into PyTorch datasets and loaders.

        Parameters
        ----------
        subset_data : pd.DataFrame
            Processed subset data with graphs and targets

        Returns
        -------
        tuple
            (train_loader, val_loader, test_loader, node_dim, edge_dim)
        """
        graphs = list(subset_data["graphs"])
        labels = torch.tensor(subset_data["target"].values, dtype=torch.long)

        train_graphs, test_graphs, train_labels, test_labels = train_test_split(
            graphs, labels, test_size=0.2, random_state=42
        )
        train_graphs, val_graphs, train_labels, val_labels = train_test_split(
            train_graphs, train_labels, test_size=0.1, random_state=42
        )

        for i, g in enumerate(train_graphs):
            g.y = train_labels[i]
        for i, g in enumerate(test_graphs):
            g.y = test_labels[i]

        train_loader = DataLoader(train_graphs, batch_size=10, shuffle=True)
        val_loader = DataLoader(val_graphs, batch_size=10, shuffle=False)
        test_loader = DataLoader(test_graphs, batch_size=10, shuffle=False)

        node_dim, edge_dim = train_graphs[0].x.shape[1], train_graphs[0].edge_attr.shape[1]

        return train_loader, val_loader, test_loader, node_dim, edge_dim

    @staticmethod
    def _train_model(
        train_loader,
        val_loader,
        node_dim: int,
        edge_dim: int,
        arch: str,
        num_layers: int,
        encoding: str,
        output_dir: str,
        exp_name: str,
        device,
    ) -> torch.nn.Module:
        """
        Train the GNN model.

        Parameters
        ----------
        train_loader : DataLoader
            Training data loader
        val_loader : DataLoader
            Validation data loader
        node_dim : int
            Node feature dimension
        edge_dim : int
            Edge feature dimension
        arch : str
            GNN architecture (GAT or GCN)
        num_layers : int
            Number of GNN layers
        encoding : str
            Node encoding method
        output_dir : str
            Directory to save model
        exp_name : str
            Experiment name for saving
        device : torch.device
            Device to train on

        Returns
        -------
        torch.nn.Module
            Trained model
        """
        model = Final_Network(
            node_dim, edge_dim, arch=arch, num_layers=num_layers, encoding=encoding
        ).to(device)
        optimizer = optim.AdamW(model.parameters(), lr=1e-5)
        scheduler = optim.lr_scheduler.OneCycleLR(
            optimizer, max_lr=1e-3, steps_per_epoch=len(train_loader), epochs=200
        )
        criterion = nn.CrossEntropyLoss()

        train_losses = []
        val_losses = []

        best_val_loss = float("inf")

        for epoch in range(200):
            model.train()
            epoch_loss = 0
            for batch in train_loader:
                batch = batch.to(device)
                optimizer.zero_grad()
                out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                loss = criterion(out, batch.y)
                loss.backward()
                optimizer.step()
                scheduler.step()
                epoch_loss += loss.item()
            train_losses.append(epoch_loss / len(train_loader))

            model.eval()
            val_loss = 0
            with torch.no_grad():
                for batch in val_loader:
                    batch = batch.to(device)
                    out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                    loss = criterion(out, batch.y)
                    val_loss += loss.item()
            val_loss /= len(val_loader)
            val_losses.append(val_loss)

            # Early stopping (optional)
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save(model.state_dict(), f"{output_dir}/{exp_name}_{arch}_best_model.pth")

        # # Plot Training Loss Curve
        # plt.figure(figsize=(6, 5))
        # plt.plot(range(1, 201), train_losses, marker='o', linestyle='-')
        # plt.xlabel("Epoch")
        # plt.ylabel("Training Loss")
        # plt.title(f"{exp_name} {arch} Training Loss (Fold {fold+1})")
        # plt.grid()
        # plt.savefig(f'{output_dir}/{exp_name}_{arch}_loss_fold{fold+1}.png')
        # plt.close()

        model.load_state_dict(torch.load(f"{output_dir}/{exp_name}_{arch}_best_model.pth"))
        return model

    @staticmethod
    def _evaluate_model(model, test_loader, device) -> tuple:
        """
        Evaluate the trained model on test data.

        Parameters
        ----------
        model : torch.nn.Module
            Trained model
        test_loader : DataLoader
            Test data loader
        device : torch.device
            Device for evaluation

        Returns
        -------
        tuple
            (y_true, y_pred, accuracy)
        """
        model.eval()
        y_true, y_pred = [], []

        with torch.no_grad():
            for batch in test_loader:
                batch = batch.to(device)
                out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                preds = torch.argmax(out, dim=1)
                y_true.extend(batch.y.cpu().numpy())
                y_pred.extend(preds.cpu().numpy())
        acc = accuracy_score(y_true, y_pred)
        return y_true, y_pred, acc

    @staticmethod
    def _plot_confusion_matrix(
        y_true, y_pred, acc, exp_name: str, arch: str, threshold: int, output_dir: str
    ) -> None:
        """
        Create and save confusion matrix plot.

        Parameters
        ----------
        y_true : list
            True labels
        y_pred : list
            Predicted labels
        acc : float
            Model accuracy score
        exp_name : str
            Experiment name
        arch : str
            GNN architecture
        threshold : int
            Classification threshold
        output_dir : str
            Directory to save plot
        """
        cm = confusion_matrix(y_true, y_pred)

        plt.figure(figsize=(8, 6))
        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap="Blues",
            xticklabels=["Not Enriched", "Enriched"],
            yticklabels=["Not Enriched", "Enriched"],
        )
        plt.title(f"{exp_name} {arch} Classifier (Acc = {acc:.2f}), Threshold = {threshold}")
        plt.xlabel("Predicted")
        plt.ylabel("True")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{exp_name}_{arch}_classifier.svg", dpi=300)
        plt.close()
