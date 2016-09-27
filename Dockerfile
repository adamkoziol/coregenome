# Dockerfile for core genome calculation
FROM ubuntu:16.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	build-essential \
	git \
	python-dev \
	python-pip \
	python-setuptools \
	wget \
	curl \
	ruby-full \	
	libdatetime-perl \ 
	libxml-simple-perl \
	libdigest-md5-perl \
	bioperl \
	roary && \
	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/

# Install QSimScan
RUN git clone https://github.com/abadona/qsimscan.git
RUN cd /qsimscan && make
ENV PATH /qsimscan:$PATH

# Install Linuxbrew
RUN useradd -m -s /bin/bash linuxbrew
RUN echo 'linuxbrew ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers

USER linuxbrew
WORKDIR /home/linuxbrew
ENV PATH /home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:$PATH
ENV SHELL /bin/bash
RUN yes |ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
RUN brew doctor || true

# Install prokka
# Enable homebrew-science tap
RUN brew tap homebrew/science
RUN brew update
# Install Prokka and all its dependencies:
RUN brew install prokka

#
ENV PATH /qsimscan/nsimscan:$PATH
ENV PATH /qsimscan/psimscan:$PATH
