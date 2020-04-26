#! /bin/bash
rm -rf data
mkdir data
python upperbound.py
python post_process.py
rm -rf *.pyc
