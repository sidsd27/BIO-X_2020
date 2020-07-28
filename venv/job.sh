#!/bin/bash
cd /Users/sidsd27/PycharmProjects/untitled/venv
while read LINE;
do
  python main.py $LINE;
  conda activate azimuth
  python /Users/sidsd27/PycharmProjects/untitled/venv/demo_doench.py
  conda activate crispr
  python new_scoring.py
done < repeat_input.txt

