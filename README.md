# Joint SeqStruct Search

**Joint SeqStruct Search** is a sequence- and structure-based retrieval workflow for viral protein analysis.

The project compares protein relationships identified by **sequence similarity** and **structural similarity**, using MMseqs2 and Foldseek on the same query/reference protein sets. The resulting retrieval results are integrated for downstream comparison and analysis.

## Data

The original Viro3D dataset and predicted protein structures can be downloaded from Zenodo:

[https://zenodo.org/records/15622906](https://zenodo.org/records/15622906)

Required files:

- `viro3d_metadata.tar.gz`
- `colabfold_pdb.tar.gz`

Direct download links:

- Metadata:  
  [https://zenodo.org/records/15622906/files/viro3d_metadata.tar.gz?download=1](https://zenodo.org/records/15622906/files/viro3d_metadata.tar.gz?download=1)

- ColabFold structures:  
  [https://zenodo.org/records/15622906/files/colabfold_pdb.tar.gz?download=1](https://zenodo.org/records/15622906/files/colabfold_pdb.tar.gz?download=1)

## Workflow

A concise description of the complete analysis workflow is provided in:

```text
notebooks/workflow.ipynb
```

The workflow contains four main stages:

1. **Dataset preparation**
   - Filter the original Viro3D dataset.
   - Extract Baculoviridae proteins.
   - Match protein sequences with predicted structures.

2. **Query and reference selection**
   - Normalize protein annotations.
   - Select representative query proteins.
   - Construct the reference protein set.

3. **Sequence and structure retrieval**
   - Sequence-based retrieval with MMseqs2.
   - Structure-based retrieval with Foldseek.

4. **Result comparison and analysis**
   - Integrate MMseqs2 and Foldseek retrieval results.
   - Compare sequence- and structure-derived protein relationships.
   - Generate downstream statistics and visualizations.

## Software

The retrieval experiments were performed using:

- **MMseqs2**: `v18.8cc5c`
- **Foldseek**: `v10.941cd33`

The preprocessing and comparison programs are implemented mainly in **C++**, while shell scripts are used to run MMseqs2 and Foldseek.

