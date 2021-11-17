# Dockerfile for building an image containing admb, ad2csv, and iscam executables with
# all linkages

FROM cgrandin/ubuntu-essentials:latest
LABEL maintainer = "chris.grandin@gmail.com"

SHELL ["/bin/bash", "-c"]

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

WORKDIR /home/container
