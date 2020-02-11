## Start from this Docker image
FROM rocker/tidyverse:3.6.0

## Install R packages in Docker image
RUN Rscript -e "library(devtools);devtools::install_github('"grst/immunedeconv"')"

## Copy files into Docker image
COPY docker_files/* ./

## Make script executable
RUN chmod a+x timer.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "timer.R"]

