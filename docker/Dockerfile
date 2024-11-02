FROM r-base

RUN apt-get update && \
  apt-get install -y r-cran-devtools libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev && \
  R -e "install.packages(c('devtools', 'testthat', 'roxygen2'))"

WORKDIR /home