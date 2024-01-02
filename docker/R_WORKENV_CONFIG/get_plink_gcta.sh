#!/bin/bash

set -e


# Create directory where we are going to keep the binaries we need
mkdir -p /home/rstudio/software

#CGTA
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
mv gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 /home/rstudio/software/
echo -e "Check gcta was installed...\n"
ls /home/rstudio/software/gcta-1.94.1

#PLINK
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20221210.zip
unzip plink_linux_x86_64_20221210.zip
mv plink /home/rstudio/software/
echo -e "Check plink was installed...\n"
ls /home/rstudio/software/plink
