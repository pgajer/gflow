#!/usr/bin/env python3
"""Generate deterministic PHATE Phase 6 baseline fixtures.

This script is a development-only fixture generator. It uses the original
Python PHATE package to freeze behavior for small examples. The R package tests
consume the resulting CSV/JSON fixtures and do not run Python.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, Dict, Iterable, Tuple

import numpy as np
import phate
import phate.mds
import phate.tree
import phate.vne
import scipy
import scipy.spatial.distance


ArrayPair = Tuple[np.ndarray, Dict[str, object]]


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def _standardize(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1
    return (X - mu) / sd


def _sample_ball(rng: np.random.Generator, n: int, dim: int) -> np.ndarray:
    z = rng.normal(size=(n, dim))
    z /= np.linalg.norm(z, axis=1)[:, None]
    r = rng.random(n) ** (1.0 / dim)
    return z * r[:, None]


def gaussian_isotropic(seed: int) -> ArrayPair:
    rng = _rng(seed)
    X = rng.normal(size=(40, 5))
    return _standardize(X), {"family": "gaussian", "n": 40, "p": 5}


def gaussian_clusters(seed: int) -> ArrayPair:
    rng = _rng(seed)
    centers = np.array([
        [-2.0, 0.0, 0.5, 0.0, 1.0, -1.0],
        [1.5, 1.8, -0.5, 1.0, -1.0, 0.0],
        [0.3, -2.2, 1.5, -1.0, 0.5, 1.3],
    ])
    sizes = [18, 22, 20]
    scales = [0.25, 0.45, 0.7]
    X = np.vstack([
        centers[i] + rng.normal(scale=scales[i], size=(sizes[i], centers.shape[1]))
        for i in range(3)
    ])
    return _standardize(X), {"family": "gaussian_clusters", "n": 60, "p": 6}


def noisy_circle_uniform(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 64
    theta = np.sort(rng.uniform(0, 2 * np.pi, size=n))
    r = 1.0 + rng.normal(scale=0.045, size=n)
    X = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    return X, {"family": "noisy_circle", "sampling": "uniform", "n": n}


def noisy_circle_nonuniform(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 72
    mix = rng.choice(3, size=n, p=[0.45, 0.35, 0.20])
    theta = np.empty(n)
    theta[mix == 0] = rng.normal(loc=0.35 * np.pi, scale=0.18, size=np.sum(mix == 0))
    theta[mix == 1] = rng.normal(loc=1.25 * np.pi, scale=0.22, size=np.sum(mix == 1))
    theta[mix == 2] = rng.uniform(0, 2 * np.pi, size=np.sum(mix == 2))
    theta %= 2 * np.pi
    r = 1.0 + rng.normal(scale=0.055, size=n)
    X = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    return X, {"family": "noisy_circle", "sampling": "nonuniform", "n": n}


def noisy_circle_outliers(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 72
    n_out = 7
    theta = rng.uniform(0, 2 * np.pi, size=n - n_out)
    r = 1.0 + rng.normal(scale=0.045, size=n - n_out)
    circle = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    outliers = rng.uniform(low=-1.8, high=1.8, size=(n_out, 2))
    X = np.vstack([circle, outliers])
    return X, {"family": "noisy_circle", "sampling": "outlier_contaminated", "n": n}


def quadform_graph(seed: int, dim: int, index: int, n: int) -> ArrayPair:
    rng = _rng(seed)
    U = _sample_ball(rng, n=n, dim=dim)
    y = np.sum(U[:, :index] ** 2, axis=1) - np.sum(U[:, index:] ** 2, axis=1)
    X = np.column_stack([U, y])
    return X, {"family": "quadratic_form_graph", "latent_dim": dim, "index": index, "n": n}


def crescent_open_arc(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 70
    theta = np.sort(rng.uniform(0.15 * np.pi, 1.65 * np.pi, size=n))
    r = 1.0 + rng.normal(scale=0.045, size=n)
    X = np.column_stack([r * np.cos(theta), 0.75 * r * np.sin(theta)])
    X += rng.normal(scale=0.01, size=X.shape)
    return X, {"family": "open_curve", "shape": "crescent", "n": n}


def asymmetric_y_tree(seed: int) -> ArrayPair:
    rng = _rng(seed)
    trunk_n, b1_n, b2_n = 28, 22, 25
    trunk_t = np.linspace(0, 1.1, trunk_n)
    b1_t = np.linspace(0.05, 1.0, b1_n)
    b2_t = np.linspace(0.05, 1.25, b2_n)
    trunk = np.column_stack([np.zeros(trunk_n), trunk_t])
    b1 = np.column_stack([0.78 * b1_t, 1.1 + 0.55 * b1_t])
    b2 = np.column_stack([-0.58 * b2_t, 1.1 + 0.88 * b2_t])
    X = np.vstack([trunk, b1, b2])
    X += rng.normal(scale=0.025, size=X.shape)
    return X, {"family": "branching", "shape": "asymmetric_y", "n": X.shape[0]}


def phate_dla_tree(seed: int) -> ArrayPair:
    X, labels = phate.tree.gen_dla(
        n_dim=12, n_branch=4, branch_length=20, rand_multiplier=1.6, seed=seed, sigma=1.2
    )
    return _standardize(X), {
        "family": "phate_paper_style_tree",
        "shape": "dla",
        "n": int(X.shape[0]),
        "labels": labels.astype(int).tolist(),
    }


def noisy_circle_gap(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 90
    theta = []
    while len(theta) < n:
        cand = rng.uniform(0, 2 * np.pi)
        if not (0.82 * np.pi < cand < 1.08 * np.pi):
            theta.append(cand)
    theta = np.array(theta)
    r = 1.0 + rng.normal(scale=0.05, size=n)
    return np.column_stack([r * np.cos(theta), r * np.sin(theta)]), {
        "family": "noisy_circle",
        "sampling": "gap",
        "n": n,
    }


def thick_annulus(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 100
    theta = rng.uniform(0, 2 * np.pi, size=n)
    r = rng.uniform(0.78, 1.22, size=n)
    X = np.column_stack([r * np.cos(theta), r * np.sin(theta)])
    return X, {"family": "annulus", "n": n}


def highdim_noisy_circle(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n, p = 90, 30
    theta = rng.uniform(0, 2 * np.pi, size=n)
    base = np.column_stack([np.cos(theta), np.sin(theta)])
    Q, _ = np.linalg.qr(rng.normal(size=(p, 2)))
    X = base @ Q.T
    X += rng.normal(scale=0.12, size=(n, p))
    X += rng.normal(scale=0.4, size=(n, 2)) @ rng.normal(size=(2, p)) * 0.15
    return _standardize(X), {"family": "highdim_noisy_circle", "n": n, "p": p}


def spiral(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 90
    t = np.linspace(0.4, 4.5 * np.pi, n)
    r = 0.12 + 0.055 * t
    X = np.column_stack([r * np.cos(t), r * np.sin(t)])
    X += rng.normal(scale=0.015, size=X.shape)
    return X, {"family": "open_curve", "shape": "spiral", "n": n}


def s_curve(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 90
    t = np.linspace(-1.4 * np.pi, 1.4 * np.pi, n)
    X = np.column_stack([np.sin(t), np.sign(t) * (np.cos(t) - 1), t / np.pi])
    X += rng.normal(scale=0.025, size=X.shape)
    return _standardize(X), {"family": "folded_manifold", "shape": "s_curve", "n": n}


def swiss_roll(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 100
    t = 1.5 * np.pi * (1 + 2 * rng.random(n))
    h = rng.uniform(-1.0, 1.0, size=n)
    X = np.column_stack([t * np.cos(t), h, t * np.sin(t)])
    X += rng.normal(scale=0.08, size=X.shape)
    return _standardize(X), {"family": "folded_manifold", "shape": "swiss_roll", "n": n}


def figure_eight(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n = 100
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    X = np.column_stack([np.sin(t), np.sin(t) * np.cos(t)])
    X += rng.normal(scale=0.018, size=X.shape)
    return X, {"family": "self_approaching_loop", "shape": "figure_eight", "n": n}


def near_crossing_strands(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n_each = 48
    x = np.linspace(-1, 1, n_each)
    strand1 = np.column_stack([x, 0.22 * x + 0.04])
    strand2 = np.column_stack([x, -0.22 * x - 0.04])
    X = np.vstack([strand1, strand2])
    X += rng.normal(scale=0.012, size=X.shape)
    return X, {"family": "near_crossing", "shape": "separate_strands", "n": X.shape[0]}


def cluster_bridge_cluster(seed: int) -> ArrayPair:
    rng = _rng(seed)
    left = rng.normal(loc=[-1.4, 0], scale=[0.18, 0.18], size=(36, 2))
    right = rng.normal(loc=[1.4, 0], scale=[0.18, 0.18], size=(36, 2))
    t = np.linspace(-1.1, 1.1, 24)
    bridge = np.column_stack([t, 0.20 * np.sin(2 * np.pi * (t + 1.1) / 2.2)])
    bridge += rng.normal(scale=0.035, size=bridge.shape)
    X = np.vstack([left, bridge, right])
    return X, {"family": "cluster_bridge_cluster", "n": X.shape[0]}


def compositional_cyclic_gradient(seed: int) -> ArrayPair:
    rng = _rng(seed)
    n, p = 90, 25
    theta = rng.uniform(0, 2 * np.pi, size=n)
    optima = np.linspace(0, 2 * np.pi, p, endpoint=False)
    logits = np.empty((n, p))
    for j, opt in enumerate(optima):
        circular_delta = np.angle(np.exp(1j * (theta - opt)))
        logits[:, j] = 1.4 * np.cos(circular_delta) + 0.25 * rng.normal(size=n)
    logits += rng.normal(scale=0.25, size=(n, 1))
    probs = np.exp(logits - logits.max(axis=1, keepdims=True))
    probs /= probs.sum(axis=1, keepdims=True)
    depths = rng.integers(2500, 7000, size=n)
    counts = np.vstack([rng.multinomial(int(depths[i]), probs[i]) for i in range(n)])
    rel = (counts + 0.5) / (counts.sum(axis=1, keepdims=True) + 0.5 * p)
    clr = np.log(rel) - np.log(rel).mean(axis=1, keepdims=True)
    return _standardize(clr), {"family": "compositional_cyclic_gradient", "n": n, "p": p}


CORE_CASES: Dict[str, Tuple[int, Callable[[int], ArrayPair]]] = {
    "gaussian_isotropic_n40_p5": (6101, gaussian_isotropic),
    "gaussian_clusters_n60_p6": (6102, gaussian_clusters),
    "noisy_circle_uniform_n64": (6103, noisy_circle_uniform),
    "noisy_circle_nonuniform_n72": (6104, noisy_circle_nonuniform),
    "noisy_circle_outliers_n72": (6105, noisy_circle_outliers),
    "quadform_graph_n2_k0": (6110, lambda seed: quadform_graph(seed, dim=2, index=0, n=60)),
    "quadform_graph_n2_k1": (6111, lambda seed: quadform_graph(seed, dim=2, index=1, n=60)),
    "quadform_graph_n2_k2": (6112, lambda seed: quadform_graph(seed, dim=2, index=2, n=60)),
    "quadform_graph_n3_k0": (6120, lambda seed: quadform_graph(seed, dim=3, index=0, n=64)),
    "quadform_graph_n3_k1": (6121, lambda seed: quadform_graph(seed, dim=3, index=1, n=64)),
    "quadform_graph_n3_k2": (6122, lambda seed: quadform_graph(seed, dim=3, index=2, n=64)),
    "quadform_graph_n3_k3": (6123, lambda seed: quadform_graph(seed, dim=3, index=3, n=64)),
    "crescent_open_arc_n70": (6130, crescent_open_arc),
    "asymmetric_y_tree_n75": (6140, asymmetric_y_tree),
    "phate_dla_tree_n80": (6150, phate_dla_tree),
}


STRESS_CASES: Dict[str, Tuple[int, Callable[[int], ArrayPair]]] = {
    "noisy_circle_gap_n90": (6201, noisy_circle_gap),
    "thick_annulus_n100": (6202, thick_annulus),
    "highdim_noisy_circle_n90_p30": (6203, highdim_noisy_circle),
    "spiral_n90": (6204, spiral),
    "s_curve_n90": (6205, s_curve),
    "swiss_roll_n100": (6206, swiss_roll),
    "figure_eight_n100": (6207, figure_eight),
    "near_crossing_strands_n96": (6208, near_crossing_strands),
    "cluster_bridge_cluster_n96": (6209, cluster_bridge_cluster),
    "compositional_cyclic_gradient_n90_p25": (6210, compositional_cyclic_gradient),
}


def _as_dense(x):
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def _write_csv(path: Path, x: np.ndarray) -> None:
    np.savetxt(path, np.asarray(x, dtype=float), delimiter=",", fmt="%.17g")


def generate_case(case_id: str, seed: int, factory: Callable[[int], ArrayPair], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    X, extra = factory(seed)
    X = np.asarray(X, dtype=float)

    params = {
        "knn": 5,
        "decay": 40,
        "t": "auto",
        "t_max": 30,
        "gamma": 1,
        "n_pca": None,
        "n_components": 2,
        "mds_solver": "smacof",
        "random_state": seed,
        "n_jobs": 1,
    }

    operator = phate.PHATE(
        knn=params["knn"],
        decay=params["decay"],
        t=params["t"],
        gamma=params["gamma"],
        n_pca=params["n_pca"],
        n_components=params["n_components"],
        mds="classic",
        mds_solver=params["mds_solver"],
        n_jobs=params["n_jobs"],
        random_state=params["random_state"],
        verbose=0,
    )
    operator.fit_transform(X)
    K = _as_dense(operator.graph.K)
    P = _as_dense(operator.graph.P)
    vne = phate.vne.compute_von_neumann_entropy(P, t_max=params["t_max"])
    t_grid = np.arange(params["t_max"])
    t_auto = int(phate.vne.find_knee_point(vne, x=t_grid))
    Pt = np.linalg.matrix_power(P, t_auto)
    U = -np.log(Pt + 1e-7)
    D_pot = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(U, "euclidean"))
    emb_classic = phate.mds.embed_MDS(
        U, ndim=2, how="classic", distance_metric="euclidean",
        solver="smacof", n_jobs=1, seed=seed, verbose=0
    )
    emb_metric = phate.mds.embed_MDS(
        U, ndim=2, how="metric", distance_metric="euclidean",
        solver="smacof", n_jobs=1, seed=seed, verbose=0
    )
    emb_nonmetric = phate.mds.embed_MDS(
        U, ndim=2, how="nonmetric", distance_metric="euclidean",
        solver="smacof", n_jobs=1, seed=seed, verbose=0
    )
    emb3_metric = phate.mds.embed_MDS(
        U, ndim=3, how="metric", distance_metric="euclidean",
        solver="smacof", n_jobs=1, seed=seed, verbose=0
    )

    _write_csv(out_dir / "X.csv", X)
    _write_csv(out_dir / "K.csv", K)
    _write_csv(out_dir / "P.csv", P)
    _write_csv(out_dir / "Pt.csv", Pt)
    _write_csv(out_dir / "U.csv", U)
    _write_csv(out_dir / "D.csv", D_pot)
    _write_csv(out_dir / "vne.csv", vne.reshape(1, -1))
    _write_csv(out_dir / "e_c.csv", emb_classic)
    _write_csv(out_dir / "e_m.csv", emb_metric)
    _write_csv(out_dir / "e_nm.csv", emb_nonmetric)
    _write_csv(out_dir / "e3_m.csv", emb3_metric)

    metadata = {
        "case_id": case_id,
        "seed": seed,
        "params": params,
        "t_auto": t_auto,
        "shape": {"n": int(X.shape[0]), "p": int(X.shape[1])},
        "example": extra,
        "versions": {
            "phate": getattr(phate, "__version__", None),
            "scipy": scipy.__version__,
            "numpy": np.__version__,
        },
    }
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)


def generate(cases: Dict[str, Tuple[int, Callable[[int], ArrayPair]]], out_root: Path) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    manifest = []
    for case_id, (seed, factory) in cases.items():
        print(f"Generating {case_id}")
        generate_case(case_id, seed, factory, out_root / case_id)
        manifest.append({"case_id": case_id, "seed": seed})
    with open(out_root / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--suite", choices=["core", "stress", "all"], default="all")
    parser.add_argument("--core-out", default="tests/testthat/fixtures/ph6/core")
    parser.add_argument("--stress-out", default="dev/phate-phase6-validation/cache/stress")
    args = parser.parse_args(argv)

    if args.suite in {"core", "all"}:
        generate(CORE_CASES, Path(args.core_out))
    if args.suite in {"stress", "all"}:
        generate(STRESS_CASES, Path(args.stress_out))


if __name__ == "__main__":
    main()
