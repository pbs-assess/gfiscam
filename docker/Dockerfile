# Docker script for ADMB containers

FROM rocker/tidyverse:latest
LABEL maintainer = "chris.grandin@gmail.com"

SHELL ["/bin/bash", "-c"]

RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install flex

ENV INST=/usr/bin
ENV ADMB_HOME=$INST/admb/build/admb
ENV ADMB_AD2CSV=$INST/admb/contrib/ad2csv
ENV ISCAM_BASE=$INST/gfiscam
ENV ISCAM_HOME=$ISCAM_BASE/build/dist/bin
ENV PATH=$ISCAM_BASE:$ISCAM_HOME:$ADMB_AD2CSV:$ADMB_HOME/bin:$HOME:$PATH

COPY install_admb.sh \
     install_iscam.sh $INST/

WORKDIR /usr/bin

RUN ./install_admb.sh -d && \
    ./install_iscam.sh

COPY bash_rc.sh /root/.bashrc
COPY i \
     clean_admb $ISCAM_BASE

ENV USER=rstudio PASSWORD=qwerty

WORKDIR /home/rstudio

# To run the container, go into the directory your document is in/will be in and run
# this (set project_name in --name=project_name to unique for each project):
#
# docker run -d -p 8787:8787 --name=ADMB --mount type=bind,source="$(pwd)",target=/home/rstudio cgrandin/admb
#
# If prompted for the username and password, they are rstudio and qwerty