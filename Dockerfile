FROM --platform=linux/amd64 python@sha256:ed86c82274b3c69b52fb5820f358f0bd7df0b603332063cb5c6e32bd220c3e6e
ENV PIP_DISABLE_PIP_VERSION_CHECK=1 PIP_NO_CACHE_DIR=1 PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg MPLCONFIGDIR=/tmp/matplotlib IPYTHONDIR=/tmp/ipython \
    XDG_CACHE_HOME=/tmp/cache PYTHONPATH=/opt/mfobs/src:/opt/mfobs
RUN apt-get update && apt-get install -y --no-install-recommends \
      libstdc++6 libgomp1 libxrender1 libx11-6 libxext6 libexpat1 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/mfobs
COPY requirements.txt .
RUN python -m pip install --extra-index-url https://py.molfoundry.de/dl/ -r requirements.txt \
    && python -m ipykernel install --name molfoundry --display-name MolFoundry
COPY . .
ARG SOURCE_COMMIT=unversioned-source
RUN python -m pip install --no-deps -e . && python -m pip check \
    && python scripts/freeze_source.py "$SOURCE_COMMIT"
ENTRYPOINT ["python", "-m", "tutorials"]
