import json
from pathlib import Path
from typing import Dict, Any, List
import torch


def save_artifacts(
    out_dir: str,
    model_state: Dict[str, Any],
    class_names: List[str],
    genes_used: List[str],
    train_subjects: List[str],
    test_subjects: List[str],
    args_dict: Dict[str, Any],
    cv_summary: Dict[str, Any],
    cv_per_fold: Dict[str, Any],
    test_metrics: Dict[str, Any],
):
    p = Path(out_dir)
    p.mkdir(parents=True, exist_ok=True)

    torch.save({
        "state_dict": model_state,
        "class_names": class_names,
        "genes_used": genes_used,
        "train_subjects": train_subjects,
        "test_subjects": test_subjects,
        "args": args_dict,
    }, p / "classifier.pt")

    with open(p / "metadata.json", "w") as f:
        json.dump({
            "class_names": class_names,
            "genes_used_count": len(genes_used),
            "train_subjects": train_subjects,
            "test_subjects": test_subjects,
            "args": args_dict,
        }, f, indent=2)

    with open(p / "metrics_cv.json", "w") as f:
        json.dump({"summary": cv_summary, "per_fold": cv_per_fold}, f, indent=2)

    with open(p / "metrics_test.json", "w") as f:
        json.dump(test_metrics, f, indent=2)


def load_classifier(model_path: str):
    pkg = torch.load(model_path, map_location="cpu")
    return pkg