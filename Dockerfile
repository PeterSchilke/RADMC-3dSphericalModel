FROM alpine:latest
RUN apk add --no-cache python3 py3-pip gfortran make musl-dev python3-dev 
RUN pip3  install astropy --break-system-packages
RUN pip3  install scipy --break-system-packages
RUN pip3  install matplotlib --break-system-packages
RUN pip3  install spectral_cube --break-system-packages
COPY radmc3d-2.0 /usr/local/radmc3d-2.0
VOLUME /data
WORKDIR /usr/local/radmc3d-2.0/src
RUN make
WORKDIR /usr/local/radmc3d-2.0/python/radmc3dPy
RUN python3 setup.py install 
WORKDIR /data
ENV PATH="$PATH:/usr/local/radmc3d-2.0/src"
