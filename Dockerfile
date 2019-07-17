FROM ubuntu:18.04

RUN RUN apt-get update -qq && apt-get install -yq build-essential \
    cmake git wget

