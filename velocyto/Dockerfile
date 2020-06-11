FROM conda/miniconda3


RUN conda install numpy scipy cython numba matplotlib scikit-learn h5py click
RUN pip install pysam

RUN apt-get update && apt-get install -y gcc
RUN pip install velocyto

