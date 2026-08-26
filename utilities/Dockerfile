FROM tacc/tacc-ubuntu18-mvapich2.3-ib

RUN apt-get -q update && apt-get install -y --no-install-recommends \
    wget libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev \
    libgdbm-dev libc6-dev libbz2-dev build-essential libffi-dev \
    python3-pip python3-dev python3-setuptools python3-wheel git \
    && cd /opt && wget https://www.python.org/ftp/python/3.12.9/Python-3.12.9.tgz \
    && tar xzf Python-3.12.9.tgz \
    && rm -rf Python-3.12.9.tgz \
    && cd Python-3.12.9 && ./configure && make install

RUN pip3 install --upgrade pip \
    && pip3 install --no-cache-dir --no-binary=mpi4py mpi4py \
    && pip3 install ucla_plha h5py

ADD ./ucla_plha_mpi.py /usr/local/bin/ucla_plha_mpi.py
ADD ./Process_H5_MPI.py /usr/local/bin/Process_H5.py

RUN chmod a+rx /usr/local/bin/ucla_plha_mpi.py \
    && chmod a+rx /usr/local/bin/Process_H5.py

ENV MV2_ENABLE_AFFINITY=0
ENV MV2_SMP_USE_CMA=0
