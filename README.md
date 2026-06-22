# WASCO

This repository proposes the WASCO algorithm in the following paper submitted to VLDB 2027.

## Source code info

Programming Language: C++17

The implementation is built with CMake (>= 3.12) and a C++17-capable compiler (e.g. `g++` or `clang++`).

## Build

```bash
cd code
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces the `wasco` executable at `code/build/wasco`.

## Dataset

The actor dataset was compressed into a ZIP file due to its large size.

Each dataset is stored in a separate directory, and the corresponding network file is named `network.dat`.

Each line in `network.dat` represents a weighted edge and consists of three values separated by a single whitespace: `(node_u, node_v, weight)`.

Here, `node_u` and `node_v` denote the two endpoint node IDs, and `weight` specifies the edge weight between them.

## Command Line Options

| Option               | Type    | Default                              | Description |
|----------------------|---------|--------------------------------------|-------------|
| `--s`                | `int`   | `15`                                  | s threshold parameter |
| `--b`                | `int`   | `60`                                 | b threshold parameter |
| `--algorithm`        | `str`   | `exp`                                | The algorithm to run. Options: `exp`, `exact`, `compare` |
| `--network`          | `str`   | `../dataset/test/new_network.dat`    | Path to the input dataset file |
| `--tactics`          | `str`   | `TTT`                                | Optimization tactics for the `exp` algorithm (Format: "XYZ"). X (T1): Candidate pruning, Y (T2): Upper bound pruning, Z (T3): Reuse Strategy (e.g., `TTT` enables all (GPA), `FFF` disables all (GIA)) |
| `--output_path`      | `str`   | `../output/results.csv`              | Path to save the CSV result file |
| `--delta_tactic`     | `str`   | `compute`                            | How to calculate delta. Options: `compute` (mee), a percentage integer (must set tactics to `TFF`), else (naive) |
| `--compare_tactic`   | `str`   | `random`                             | Heuristic for the `compare` algorithm. Options: `degree`, `high_degree`, `weight_sum`, `high_weight_sum` |
| `--calculating_iter` | `str`   | `F`                                  | `T` to track detailed iteration logs, `F` for a standard run |

## Usage

Run the binary from the `code` directory (so the default relative dataset/output paths resolve correctly):

```bash
./build/wasco --network ../dataset/real/final/moreno_names/network.dat --algorithm exp --s 15 --b 60 --tactics TTT --output_path ../output/GPA.csv
```
```bash
./build/wasco --network ../dataset/real/final/moreno_names/network.dat --algorithm compare --s 15 --b 60 --compare_tactic degree --output_path ../output/compare.csv
```
```bash
./build/wasco --network ../dataset/real/exact/karate/network.dat --algorithm exact --s 5 --b 4 --output_path ../output/exact.csv
```
