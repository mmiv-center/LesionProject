FROM ubuntu:18.04

RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget boost

# install itk
RUN cd /tmp/ \
    && wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.0.0/InsightToolkit-5.0.0.tar.gz \
    && cd /opt/ && tar xzvf /tmp/InsightToolkit-5.0.0.tar.gz && cd /opt/InsightToolkit-5.0.0 \
    && mkdir bin && cd bin && cmake .. && make

RUN mkdir /ConnectedComponents && cd /ConnectedComponents/ \
    && git clone https://github.com/mmiv-center/LesionProject.git . \
    && cmake .

ENTRYPOINT [ "/ConnectedComponents/ConnectedComponents" ]
