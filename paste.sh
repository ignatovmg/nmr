paste <(sort -k1 distance_filtered) <(sort -k1 spring_min.log | cut -d' ' -f2) | LC_ALL=C sort -k2 -n > combo
