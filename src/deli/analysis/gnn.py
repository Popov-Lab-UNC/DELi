import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GINEConv, GATv2Conv, global_add_pool
from torch_geometric.data import Data, Dataset
from rdkit import Chem
from tqdm import tqdm

# --- GNN code implemented from itakigawa's resources on GitHub ---

# Define the convolution layer based on the architecture (GIN or GAT)
def MyConv(node_dim, edge_dim, arch='GIN'):
    """Creates a graph convolutional layer based on the chosen architecture."""
    if arch == 'GIN':
        h = nn.Sequential(nn.Linear(node_dim, node_dim, bias=True))
        return GINEConv(h, edge_dim=edge_dim)
    elif arch == 'GAT':
        return GATv2Conv(node_dim, node_dim, edge_dim=edge_dim)
    else:
        raise ValueError(f"Unknown architecture: {arch}")

class MyGNN(nn.Module):
    """Stacked GNN layers with activation."""
    def __init__(self, node_dim, edge_dim, arch='GIN', num_layers=3):
        super().__init__()
        self.convs = nn.ModuleList([MyConv(node_dim, edge_dim, arch) for _ in range(num_layers)])

    def forward(self, x, edge_index, edge_attr):
        for conv in self.convs:
            x = conv(x, edge_index, edge_attr)
            x = F.leaky_relu(x)
        return x


class Final_Network(nn.Module):
    """Full GNN-based classifier with feature encoding and prediction head."""
    def __init__(self, node_dim, edge_dim, arch='GIN', num_layers=3, encoding='onehot'):
        super().__init__()

        self.encoding = encoding

        if encoding != 'onehot':
            self.atom_encoder = nn.Embedding(num_embeddings=118+1, embedding_dim=64)
            self.bond_encoder = nn.Embedding(num_embeddings=21+1, embedding_dim=8)
            node_dim = (node_dim - 1) + 64
            edge_dim = (edge_dim - 1) + 8
        else:
            node_dim = (node_dim - 1) + 118 + 1
            edge_dim = (edge_dim - 1) + 21 + 1

        self.gnn = MyGNN(node_dim, edge_dim, arch, num_layers=num_layers)

        embed_dim = node_dim // 2
        self.head = nn.Sequential(
            nn.BatchNorm1d(node_dim),
            nn.Dropout(p=0.5),
            nn.Linear(node_dim, embed_dim, bias=True),
            nn.ReLU(),
            nn.BatchNorm1d(embed_dim),
            nn.Dropout(p=0.5),
            nn.Linear(embed_dim, 2)
        )

    def forward(self, x, edge_index, edge_attr, batch):
        """Forward pass through encoding, GNN, and prediction head."""
        if self.encoding == 'onehot':
            x0 = F.one_hot(x[:, 0].to(torch.int64), num_classes=118 + 1)
            edge_attr0 = F.one_hot(edge_attr[:, 0].to(torch.int64), num_classes=21 + 1)
        else:
            x0 = self.atom_encoder(x[:, 0].int())
            edge_attr0 = self.bond_encoder(edge_attr[:, 0].int())

        # Combine embeddings with remaining node/edge features
        x = torch.cat([x0, x[:, 1:]], dim=1)
        edge_attr = torch.cat([edge_attr0, edge_attr[:, 1:]], dim=1)

        # Pass through GNN
        node_out = self.gnn(x, edge_index, edge_attr)
        graph_out = global_add_pool(node_out, batch)  # Graph-level pooling

        # Classification head
        return self.head(graph_out)


# Feature extraction functions
def atom_feature(atom):
    return [atom.GetAtomicNum(), atom.GetDegree(), atom.GetNumImplicitHs(), atom.GetIsAromatic()]

def bond_feature(bond):
    return [bond.GetBondTypeAsDouble(), bond.GetStereo()]


# Function to convert SMILES to PyG data format
def smi_to_pyg(smi, y, encoding='embedding'):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    id_pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    atom_pairs = [z for (i, j) in id_pairs for z in [(i, j), (j, i)]]
    
    bonds = [mol.GetBondBetweenAtoms(i, j) for (i, j) in atom_pairs]

    if encoding == 'onehot':
        atom_features = [[a.GetAtomicNum()] for a in mol.GetAtoms()]
        bond_features = [[b.GetBondTypeAsDouble()] for b in bonds]
    else:
        atom_features = [a.GetAtomicNum() for a in mol.GetAtoms()]
        bond_features = [b.GetBondTypeAsDouble() for b in bonds]

    return Data(
        edge_index=torch.LongTensor(list(zip(*atom_pairs))),
        x=torch.LongTensor(atom_features).unsqueeze(1),  # Ensures shape compatibility
        edge_attr=torch.LongTensor(bond_features).unsqueeze(1),
        y=torch.LongTensor([y]),
        mol=mol,
        smiles=smi
    )

class MyDataset(Dataset):
    def __init__(self, smiles, response, encoding='embedding'):
        mols = [smi_to_pyg(smi, y, encoding) for smi, y in tqdm(zip(smiles, response), total=len(smiles))]
        self.X = [m for m in mols if m]

    def __getitem__(self, idx):
        return self.X[idx]

    def __len__(self):
        return len(self.X)
