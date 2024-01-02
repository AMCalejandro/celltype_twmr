FROM bioconductor/bioconductor_docker:devel

RUN apt-get update && \
    apt-get install -y \
    git-core \
    libcurl4-openssl-dev \
    libgit2-dev \
    libicu-dev \
    libssl-dev \
    make pandoc \
    pandoc-citeproc \
    zlib1g-dev \
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libsbml5-dev \
	qpdf \
	cmake \
	&& apt-get clean \
    && rm -rf /var/lib/apt/lists/*


ADD . /tmp/
WORKDIR /tmp/

RUN chmod +x R_WORKENV_CONFIG/*

RUN bash R_WORKENV_CONFIG/get_pandoc.sh
RUN bash R_WORKENV_CONFIG/install_r_deps.sh NCPUS=30
RUN bash R_WORKENV_CONFIG/get_plink_gcta.sh

ENV PATH "$PATH:/home/rstudio/software"

RUN cd ../tmp
RUN rm -rf /tmp/*
