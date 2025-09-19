from typing import Dict, Any
from sklearn.metrics import (
    accuracy_score,
    precision_recall_fscore_support,
    roc_auc_score,
)
import numpy as np


def classification_metrics(y_true, y_pred, y_prob, class_names) -> Dict[str, Any]:
    """Compute a bundle of metrics (macro + weighted) and multiclass AUC (OVR macro).
    y_prob: (n_samples, n_classes) or None -> AUC only if probabilities provided.
    """
    out: Dict[str, Any] = {}
    out["accuracy"] = float(accuracy_score(y_true, y_pred))

    # macro / weighted precision, recall, f1
    for avg in ("macro", "weighted"):
        p, r, f1, _ = precision_recall_fscore_support(y_true, y_pred, average=avg, zero_division=0)
        out[f"precision_{avg}"] = float(p)
        out[f"recall_{avg}"] = float(r)
        out[f"f1_{avg}"] = float(f1)

    # AUC (macro OVR) if probabilities available and >1 class
    if y_prob is not None and len(class_names) > 1:
        try:
            out["auc_ovr_macro"] = float(roc_auc_score(y_true, y_prob, multi_class="ovr", average="macro"))
        except Exception:
            out["auc_ovr_macro"] = float("nan")
    else:
        out["auc_ovr_macro"] = float("nan")

    return out