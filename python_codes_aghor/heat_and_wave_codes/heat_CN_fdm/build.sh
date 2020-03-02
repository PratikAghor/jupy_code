#! /bin/bash
rm -rf data
mkdir data
python heat_CN_main.py
python post_process.py
rm -rf heat_CN_functions.pyc
