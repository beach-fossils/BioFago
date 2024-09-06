#!/bin/bash

# Function to compare version numbers
version_compare() {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # Fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # Fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}

# Check if BLAST is already installed
if command -v blastn &> /dev/null
then
    installed_version=$(blastn -version | awk 'NR==1{print $2}' | cut -d'+' -f1)
    required_version="2.15.0"
    echo "BLAST version $installed_version is installed"

    version_compare $installed_version $required_version
    case $? in
        0) echo "Installed BLAST version is up to date" ; exit 0 ;;
        1) echo "Installed BLAST version is newer than required" ; exit 0 ;;
        2) echo "Installed BLAST version is older than required. Proceeding with update." ;;
    esac
else
    echo "BLAST is not installed. Proceeding with installation."
fi

# Detect OS and install/update BLAST
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    sudo apt-get update
    sudo apt-get install -y ncbi-blast+
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    if command -v brew &> /dev/null
    then
        brew install blast
    else
        echo "Homebrew is not installed. Please install Homebrew first."
        exit 1
    fi
else
    echo "Unsupported OS. Please install BLAST manually."
    exit 1
fi

# Verify installation
if command -v blastn &> /dev/null
then
    new_version=$(blastn -version | awk 'NR==1{print $2}' | cut -d'+' -f1)
    echo "BLAST has been successfully installed/updated to version $new_version"
else
    echo "BLAST installation/update failed"
    exit 1
fi