# Use the official Ubuntu 20.04 as a base image
FROM ubuntu:20.04

# Set environment variable to force IPv4 to avoid issues with IPv6 DNS resolution
RUN echo 'Acquire::ForceIPv4 "true";' > /etc/apt/apt.conf.d/99force-ipv4

# Set environment variables for non-interactive apt-get installs and timezone configuration
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Lisbon

# Install tzdata package to avoid interactive prompt
RUN apt-get update --fix-missing && apt-get install -y tzdata

# Remove existing apt lists and update package lists
RUN rm -rf /var/lib/apt/lists/* \
    && apt-get update --fix-missing \
    && apt-get install -y \
        python3 \
        python3-pip \
        wget \
        curl \
        pkg-config \
        libhdf5-dev

# Set environment variables for non-interactive apt-get installs
RUN apt-get clean && rm -rf /var/lib/apt/lists/* \
    && sed -i 's/http:\/\/archive.ubuntu.com/http:\/\/ports.ubuntu.com/g' /etc/apt/sources.list \
    && apt-get update --fix-missing

# Install necessary Python packages
RUN pip3 install --no-cache-dir planemo

# Copy the local project directory to the Docker container
COPY . /usr/src/app

# Set the working directory in the Docker container
WORKDIR /usr/src/app

# Install any additional Python dependencies specified in requirements.txt
RUN if [ -f requirements.txt ]; then pip3 install --no-cache-dir -r requirements.txt; fi

# Expose the port that Planemo serves on (adjust as needed)
EXPOSE 8000

# Define the default command to run when starting the container
CMD ["planemo", "test"]























