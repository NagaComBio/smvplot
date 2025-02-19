FROM mambaorg/micromamba

LABEL smvplot="An IGV-like screenshot of small variants from SAM or BAM file"

USER root
RUN useradd smvplot
RUN mkdir /build && chown smvplot:smvplot /build
USER smvplot

ENV C_INCLUDE_PATH=/usr/local/include
ENV LIBRARY_PATH=/usr/local/lib
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV PATH=$PATH:/opt/conda/bin/

WORKDIR /build
COPY --chown=smvplot:smvplot . /build

RUN micromamba install -y -n base -f /build/containers/smvplot.yaml \
    && micromamba clean --all --yes

# Install build dependencies
RUN micromamba install -y -n base pip
RUN pip install build

# Build and install the Python package
RUN python -m build && pip install dist/*tar.gz
