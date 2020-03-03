#! /bin/bash
rm -rf data
mkdir data
python Burger_AB2_main.py
python Burger_AB4_main.py
python post_process.py
rm -rf Burger_AB_functions.pyc
rm -rf params.pyc
