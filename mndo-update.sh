#!/bin/bash
url=https://mndo.blusteve.com
wget "$url/MNDOParam.jar" -O MNDOParam.jar
wget "$url/atom-properties.json" -O atom-properties.json
wget "$url/molecules.txt" -O molecules.txt
wget "$url/params.csv" -O params.csv
wget "$url/param-numbers.csv" -O param-numbers.csv
