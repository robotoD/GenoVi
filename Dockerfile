# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

WORKDIR /app

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

ADD ./circos-0.69-9.tgz /app/
RUN tar -xzvf circos-0.69-9.tgz
ENV PATH="$PATH:/app/circos-0.69-9/bin/"

RUN pip3 install -e .

RUN apt update && apt install -y cmake \
    cpanminus \
    pkg-config \
    libgd-graph-perl

RUN cpanm Clone \
      Font::TTF::Font \
      Regexp::Common \
      SVG \
      Math::VecStat \
      Math::Bezier \
      Text::Format \
      Math::Round \
      Config::General \
      Font::TTF::Font \
      GD \
      GD::Polyline \
      Math::Bezier \
      Math::Round \
      Math::VecStat \
      Regexp::Common \
      SVG \
      Set::IntSpan \
      Statistics::Basic \
      Text::Format

