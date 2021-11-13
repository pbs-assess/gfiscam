# Docker script for ADMB containers
#
# This script was based loosely on this Dockerfile:
# https://github.com/rocker-org/rocker-versioned/blob/master/verse/3.6.3.Dockerfile

# Base rocker/tidyverse image https://hub.docker.com/u/rocker/
FROM rocker/tidyverse:latest
LABEL maintainer = "chris.grandin@gmail.com"

RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install flex && \
    apt-get install bash-completion && \
    apt-get install htop

ENV INST=/usr/bin
ENV ADMB_HOME=$INST/admb/build/admb
ENV ADMB_AD2CSV=$INST/admb/contrib/ad2csv
ENV ISCAM_BASE=$INST/gfiscam
ENV ISCAM_HOME=$ISCAM_BASE/build/dist/bin
ENV PATH=$ISCAM_BASE:$ISCAM_HOME:$ADMB_AD2CSV:$ADMB_HOME/bin:$HOME:$PATH

COPY install_admb.sh $INST/install_admb.sh
COPY install_iscam.sh $INST/install_iscam.sh

RUN chmod +x $INST/install_admb.sh && \
    chmod +x $INST/install_iscam.sh && \
    cd $INST && ./install_admb.sh -d && \
    cd $INST && ./install_iscam.sh

COPY bash_rc.sh /root/.bashrc
COPY i $ISCAM_BASE/i
COPY clean_admb $ISCAM_BASE/clean_admb

ENV USER=rstudio PASSWORD=qwerty

# To run the container, go into the directory your document is in/will be in and run
# this (set project_name in --name=project_name to unique for each project):
#
# docker run -d -p 8787:8787 --name=ADMB --mount type=bind,source="$(pwd)",target=/home/rstudio cgrandin/admb
#
# If prompted for the username and password, they are rstudio and qwerty