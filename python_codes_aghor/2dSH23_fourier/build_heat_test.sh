#! /bin/bash
rm -rf heat_data
mkdir heat_data
python test_2d_heat.py
python test_post_process_heat.py
rm -rf *.pyc
