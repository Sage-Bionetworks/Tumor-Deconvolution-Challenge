## Start from this Docker image
FROM rocker/tidyverse:3.6.0

RUN Rscript -e "install.packages('e1071')"
RUN Rscript -e "install.packages('BiocManager')"

RUN Rscript -e "BiocManager::install('preprocessCore')"

## Copy files into Docker image
COPY docker_files/* ./

## Make script executable
RUN chmod a+x do_cibersort.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "do_cibersort.R"]


