#!/usr/bin/bash
#SBATCH -p defq
#SBATCH -t 24:00:00
#SBATCH -N 3

date
python3 sri_eigenvalue_3d.py --re=5182 --eta=0.6 --G=50 --epsilon=0 --Pr=4.35 --m=6 --k=3 --mu=0 --ar=3 --alpha_z=3.13

date
