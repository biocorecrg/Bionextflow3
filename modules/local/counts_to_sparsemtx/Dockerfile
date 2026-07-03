FROM python:3.12-slim
LABEL description="Convert STAR ReadsPerGene to Cell Ranger sparse matrix format"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    pandas \
    scipy \
    anndata



