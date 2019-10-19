#!/bin/bash

set -xe

export PATH="$HOME/miniconda/bin:$PATH"
source $HOME/miniconda/bin/activate
#
cd Usecases/test && sh run-all.sh && cd -

