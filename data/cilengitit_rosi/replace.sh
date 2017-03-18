awk 'FNR==NR { a[$1]=$2; next } $3 in a { $3=a[$3] } $4 in a { $4=a[$4] }1' replace.csv rosi_tvol.csv | sed 's/ /      /g' > rosi_tvol_reid.csv
