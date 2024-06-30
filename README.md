# Dataset
The SignalP6.0 dataset focuses on 6 types of signaling protein.

This project will help crawl 3D data from AlphaFold with _protein_ID_ from SignalP6.0 and transform to Graph 3D representation.

# Instruction
1. Crawling 3D data from AlphaFold 
- If you are using Kaggle, import notebook `protein3dcrawler.ipynb` into the Kaggle coding environment and run it.
- If you are using local, just run `crawler.py`.

2. Transforming 3D data to Graph 3D.
- If you are using Kaggle, import notebook `protein3dbuilder.ipynb` into the Kaggle coding environment and run it.
- If you are using local, just run `utils.py`.

# Note
1. If you are using Kaggle for crawling, remember import SignalP6.0 as input dataset.
2. If you are using Kaggle for transforming/building, remember import [`3d_dataset.json`](https://drive.google.com/file/d/1qDjdufOXU4lAZNs-dpYqcLMofqmMLIoj/view?usp=sharing) and [`smiles_string_aa.csv`](https://drive.google.com/file/d/14vYaDhmD6m1DOHKpc1H2cXm-XhXeZHmF/view?usp=sharing) as input datasets.
