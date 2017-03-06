#!/usr/bin/env bash

## Initialization script for setting up directory structure, downloading
## HMMER if necessary, and configuring system architecture for use with
## Meta-MARC.

shopt -s extglob


#############
## Methods ##
#############
config_metamarc() {
    local OS=${2}
    local ARCH=${3}
    local RELPATH=${1}
    local HMMER="n"
    
    echo -n "Download HMMER binaries to the meta-marc directory? (yes/no) [no]: "
    read HMMER
    
    if [ "$HMMER" == "yes" ]; then
        if [ "$OS" == "Darwin" ]; then
            wget -P "$RELPATH" http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-macosx-intel.tar.gz
            tar -xzvf "${RELPATH}/hmmer-3.1b2-macosx-intel.tar.gz" -C "${RELPATH}"
            mv "${RELPATH}/hmmer-3.1b2-macosx-intel" "${RELPATH}/hmmer-3.1b2"
            rm "${RELPATH}/hmmer-3.1b2-macosx-intel.tar.gz"
        elif [ "$OS" == "Linux" ] && [ "$ARCH" == "64" ]; then
            wget -P "$RELPATH" http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
            tar -xzvf "${RELPATH}/hmmer-3.1b2-linux-intel-x86_64.tar.gz" -C "${RELPATH}"
            mv "${RELPATH}/hmmer-3.1b2-linux-intel-x86_64" "${RELPATH}/hmmer-3.1b2"
            rm "${RELPATH}/hmmer-3.1b2-linux-intel-x86_64.tar.gz"
        elif [ "$OS" == "Linux" ] && [ "$ARCH" == "32" ]; then
            wget -P "$RELPATH" http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-ia32.tar.gz
            tar -xzvf "${RELPATH}/hmmer-3.1b2-linux-intel-ia32.tar.gz" -C "${RELPATH}"
            mv "${RELPATH}/hmmer-3.1b2-linux-intel-ia32" "${RELPATH}/hmmer-3.1b2"
            rm "${RELPATH}/hmmer-3.1b2-linux-intel-ia32.tar.gz"
        fi
        cd bin
        ln -s ../hmmer-3.1b2/binaries/* .
        cd ..
    fi
}


##########
## Main ##
##########
OS=""
ARCH=""
RELPATH="${BASH_SOURCE%/*}"

chmod +x ${RELPATH}/bin/*

if [ "$( uname )" == "Darwin" ]; then
    OS="Darwin"
elif [ "$( expr substr $( uname -s ) 1 5 )" == "Linux" ] && [ $( uname -a | grep -o -m 1 "x86_64" | head -1 ) == "x86_64" ]; then
    OS="Linux"
    ARCH="64"
elif [ "$( expr substr $( uname -s ) 1 5 )" == "Linux" ] && [ $( uname -a | grep -o -m 1 "i386" | head -1 ) == "i386" ]; then
    OS="Linux"
    ARCH="32"
else
    echo -e "\nError: Not a valid system architecture for Meta-MARC.  Linux 32/64-bit or Mac OS required."
    exit 1
fi

config_metamarc ${RELPATH} ${OS} ${ARCH}

gunzip ${RELPATH}/src/HMMs/*

exit 0

