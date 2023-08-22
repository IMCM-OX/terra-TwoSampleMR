FROM mrcieu/twosamplemr:0.5.7

LABEL org.opencontainers.image.source=https://github.com/anand-imcm/terra-TwoSampleMR-wf1
LABEL org.opencontainers.image.description="An image containing tools and packages for the two sample MR analysis."

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y --no-install-recommends \
    libudunits2-dev \
    libgdal-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --deps TRUE --error --skipinstalled --ncpus -1 \
    optparse \
    qpdf \
    berryFunctions

RUN installGithub.r --deps TRUE --repos https://cloud.r-project.org \
    NightingaleHealth/ggforestplot

COPY scripts /scripts

ENV PATH "/scripts/:$PATH"
