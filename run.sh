# Extract all the Nitrogen atoms inside the polygen
python3.4  extract_elements_inside_polygon.py  < input_extract_elements_inside_polygon
xmgrace -block ts_N_990.dat -bxy 3:4 | xmgrace pol_N_990.dat
