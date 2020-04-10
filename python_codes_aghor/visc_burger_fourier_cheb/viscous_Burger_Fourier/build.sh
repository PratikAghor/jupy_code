#! /bin/bash
rm -rf data
mkdir data
python burger_fourier_PratikAghor.py
python post_process.py
rm -rf *.pyc
