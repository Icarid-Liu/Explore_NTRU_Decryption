FROM sagemath/sagemath:10.1

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        clang \
        git \
        libssl-dev \
        make \
        valgrind \
    && rm -rf /var/lib/apt/lists/*

USER sage
ENV PATH="/home/sage/.local/bin:${PATH}"
WORKDIR /home/sage/artifact

COPY --chown=sage:sage . .

RUN sage -pip install --user --no-cache-dir \
        tqdm==4.67.1 \
        jupyter==1.0.0 \
        nbconvert==7.16.4 \
    && bash artifact/setup_estimators.sh

CMD ["bash"]
