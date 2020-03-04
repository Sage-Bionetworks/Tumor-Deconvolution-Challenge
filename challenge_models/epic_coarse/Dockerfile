## Start from this Docker image
FROM rocker/tidyverse:3.6.1

## Install R packages in Docker image
RUN Rscript -e "library(devtools);devtools::install_github('GfellerLab/EPIC')"

## Copy files into Docker image
COPY docker_files/* ./

## Make script executable
RUN chmod a+x epic.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "epic.R"]

