# TIRPClo: efficient and complete mining of time intervals‑related patterns

## Introduction
Welcome!

This repository is related to the paper "TIRPClo: efficient and complete mining of time intervals‑related patterns", published in *Data Mining and Knowledge Discovery, 2023*.

Please cite as:

```bibtex
@article{harel2023tirpclo,
  title={TIRPClo: Efficient and complete mining of time intervals-related patterns},
  author={Harel, Omer and Moskovitch, Robert},
  journal={Data Mining and Knowledge Discovery},
  volume={37},
  number={5},
  pages={1806--1857},
  year={2023},
  publisher={Springer}
}
```

## Dependencies

- Python 3.8.10
- Packages:
    ```
    pytest==8.3.1
    pre-commit==3.5.0
    ```

## Run

```commandline
python -m tirpclo.run -n {num-entities} -s {min-support} -g {maximal-gap} -f 'datasets/{path-to-dataset}'
```
