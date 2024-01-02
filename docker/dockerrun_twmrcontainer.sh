#!/bin/env bash

docker run --env PASSWORD=introduce_your_pwd \
--publish 8080:8787 \
--detach \
--name amc-celltyping \
-v /mnt/acarrasco/RDS/acarrasco/ANALYSES_WORKSPACE/CELLTYPING_PROJECT/TWMR_DATA/:/home/rstudio/DATA \
-v /mnt/acarrasco/RDS/acarrasco/ANALYSES_WORKSPACE/CELLTYPING_PROJECT/CELLTYPING_V2:/home/rstudio/celltyping_wd \
amc-celltyping:latest
