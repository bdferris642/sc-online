from __future__ import annotations
from typing import Tuple, Dict, Any, List
import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import LabelEncoder
from .metrics import classification_metrics


class DropoutMLP(nn.Module):
    def __init__(self, in_dim: int, n_classes: int, hidden: List[int], dropout: float):
        super().__init__()
        layers: List[nn.Module] = []
        prev = in_dim
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(p=dropout)]
            prev = h
        layers += [nn.Linear(prev, n_classes)]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def _to_tensor(X):
    if sp.issparse(X):
        X = X.toarray()
    return torch.from_numpy(np.asarray(X, dtype=np.float32))


def set_seed(seed: int):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def train_eval_single_split(
    X_cxg,
    y_labels: np.ndarray,
    groups: np.ndarray,
    train_subjects: List[str],
    hidden: List[int],
    dropout: float,
    epochs: int,
    batch_size: int,
    lr: float,
    seed: int = 42,
    device: str | None = None,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Any, Any]:
    """Train on `train_subjects` and evaluate on held-out subjects.

    Returns (cv_summary, cv_per_fold, test_metrics, label_encoder, best_model_state_dict)
    """
    set_seed(seed)
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(y_labels)
    class_names = le.classes_.tolist()

    in_train = np.isin(groups, train_subjects)
    X_train = X_cxg[in_train]
    y_train = y[in_train]
    g_train = groups[in_train]

    in_test = ~in_train
    X_test = X_cxg[in_test]
    y_test = y[in_test]

    X_all_tensor = _to_tensor(X_train)
    y_all_tensor = torch.from_numpy(y_train.astype(np.int64))

    n_classes = len(class_names)
    in_dim = X_all_tensor.shape[1]

    # Cross-validation (GroupKFold on subjects within train)
    gkf = GroupKFold(n_splits=min(5, len(np.unique(g_train))))
    fold_acc, fold_f1, fold_auc = [], [], []

    for tr_idx, va_idx in gkf.split(X_train, y_train, groups=g_train):
        Xtr, Xva = X_all_tensor[tr_idx], X_all_tensor[va_idx]
        ytr, yva = y_all_tensor[tr_idx], y_all_tensor[va_idx]

        train_ds = TensorDataset(Xtr, ytr)
        val_ds = TensorDataset(Xva, yva)
        train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)

        model = DropoutMLP(in_dim=in_dim, n_classes=n_classes, hidden=hidden, dropout=dropout).to(device)

        # class weights (handle imbalance)
        binc = np.bincount(ytr.numpy(), minlength=n_classes)
        binc[binc == 0] = 1
        weights = torch.tensor((binc.sum() / binc), dtype=torch.float32, device=device)
        criterion = nn.CrossEntropyLoss(weight=weights)
        optimizer = torch.optim.AdamW(model.parameters(), lr=lr)

        model.train()
        for _ in range(epochs):
            for xb, yb in train_loader:
                xb, yb = xb.to(device), yb.to(device)
                optimizer.zero_grad()
                logits = model(xb)
                loss = criterion(logits, yb)
                loss.backward()
                optimizer.step()

        # Validate
        model.eval()
        with torch.no_grad():
            logits = []
            ys = []
            for xb, yb in val_loader:
                xb = xb.to(device)
                lg = model(xb).cpu().numpy()
                logits.append(lg)
                ys.append(yb.numpy())
        yva_true = np.concatenate(ys)
        yva_log = np.vstack(logits)
        yva_prob = softmax(yva_log)
        yva_pred = yva_prob.argmax(axis=1)

        m = classification_metrics(yva_true, yva_pred, yva_prob, class_names)
        fold_acc.append(m["accuracy"]) ; fold_f1.append(m["f1_macro"]) ; fold_auc.append(m["auc_ovr_macro"])

    cv_summary = {
        "accuracy_mean": float(np.mean(fold_acc)) if len(fold_acc) else float("nan"),
        "accuracy_sd": float(np.std(fold_acc, ddof=1)) if len(fold_acc) > 1 else 0.0,
        "f1_macro_mean": float(np.mean(fold_f1)) if len(fold_f1) else float("nan"),
        "f1_macro_sd": float(np.std(fold_f1, ddof=1)) if len(fold_f1) > 1 else 0.0,
        "auc_ovr_macro_mean": float(np.mean(fold_auc)) if len(fold_auc) else float("nan"),
        "auc_ovr_macro_sd": float(np.std(fold_auc, ddof=1)) if len(fold_auc) > 1 else 0.0,
    }
    cv_per_fold = {
        "accuracy": fold_acc,
        "f1_macro": fold_f1,
        "auc_ovr_macro": fold_auc,
    }

    # Final fit on all training subjects
    train_loader = DataLoader(TensorDataset(X_all_tensor, y_all_tensor), batch_size=batch_size, shuffle=True)
    model = DropoutMLP(in_dim=in_dim, n_classes=n_classes, hidden=hidden, dropout=dropout).to(device)
    binc = np.bincount(y_all_tensor.numpy(), minlength=n_classes)
    binc[binc == 0] = 1
    weights = torch.tensor((binc.sum() / binc), dtype=torch.float32, device=device)
    criterion = nn.CrossEntropyLoss(weight=weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)
    model.train()
    for _ in range(epochs):
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            logits = model(xb)
            loss = criterion(logits, yb)
            loss.backward()
            optimizer.step()

    # Test on held-out subjects
    X_test_tensor = _to_tensor(X_test)
    test_loader = DataLoader(TensorDataset(X_test_tensor, torch.from_numpy(y_test.astype(np.int64))), batch_size=batch_size, shuffle=False)
    model.eval()
    with torch.no_grad():
        logits = []
        ys = []
        for xb, yb in test_loader:
            xb = xb.to(device)
            lg = model(xb).cpu().numpy()
            logits.append(lg)
            ys.append(yb.numpy())
    y_true = np.concatenate(ys) if len(ys) else np.array([], dtype=int)
    y_log = np.vstack(logits) if len(logits) else np.zeros((0, n_classes), dtype=np.float32)
    y_prob = softmax(y_log) if y_log.size else None
    y_pred = y_prob.argmax(axis=1) if y_prob is not None else np.array([], dtype=int)

    test_metrics = classification_metrics(y_true, y_pred, y_prob, class_names)

    return cv_summary, cv_per_fold, test_metrics, le, model.state_dict()


def softmax(logits: np.ndarray) -> np.ndarray:
    # numerically stable softmax
    z = logits - logits.max(axis=1, keepdims=True)
    exp = np.exp(z)
    return exp / exp.sum(axis=1, keepdims=True)
