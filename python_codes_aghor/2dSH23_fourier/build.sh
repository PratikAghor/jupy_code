#! /bin/bash
bash clean.sh
#
mkdir data
mkdir snapshots
#
python 2dSH23_fourier_PratikAghor.py
python post_process.py
rm -rf *.pyc
