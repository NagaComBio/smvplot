FROM mambaorg/micromamba

LABEL smvplot="An IGV-like screenshot of small variants from SAM or BAM file"

ENV C_INCLUDE_PATH=/usr/local/include
ENV LIBRARY_PATH=/usr/local/lib
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV PATH=$PATH:/opt/conda/bin/

COPY smvplot.yaml /opt/

RUN micromamba install -y -n base -f /opt/smvplot.yaml \
    && micromamba clean --all --yes

