#!/bin/bash
#$ -cwd
#$ -l m_mem_free=2G
#$ -pe threads 16

python analyse.py

