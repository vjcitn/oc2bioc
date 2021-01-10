# container should be built in root folder of oc2bioc R package source
# THIS DOES NOT WORK -- something wrong with devtools::install command
# ANYWAY THE DESCRIPTION.md file can be used with a 'community cloud environment' and
# no customization is really needed

FROM us.gcr.io/anvil-gcr-public/anvil-rstudio-bioconductor:0.0.9
# intended for AnVIL with Rstudio

# This is to avoid the error
# 'debconf: unable to initialize frontend: Dialog'
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update ; \
    apt-get -y install libharfbuzz-dev libfribidi-dev 


RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories(), upgrade='never')"


# Add back other env vars
RUN echo "TERRA_R_PLATFORM='anvil-rstudio-bioconductor'" >> /usr/local/lib/R/etc/Renviron.site \
    && echo "TERRA_R_PLATFORM_BINARY_VERSION='0.99.1'" >> /usr/local/lib/R/etc/Renviron.site

USER root

# Init command for s6-overlay
CMD ["/init"]

