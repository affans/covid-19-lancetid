#!/bin/bash
echo "fetching full cases data..."
curl https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv > data/csv/COVID_ontario_data.csv
echo "done."
echo "fetching google distancing data..."
curl https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv?cachebust=fbb340e43a0602e1 | { head -1; grep ",CA-ON,"; } > data/csv/distancing_data.csv
echo "done."