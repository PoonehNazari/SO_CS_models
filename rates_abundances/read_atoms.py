import os

##########################################################################################
##########################################################################################
def read_extract_specs():
    myfile = open('rate12_full_atomic_v2.specs', 'r')
    lines = myfile.readlines()
    myfile.close()
    specs = []
    for line in lines:
        parts = line.split()
        second_column_value = parts[1]
        specs.append(second_column_value)

    specs_string = ' '.join(specs)
    os.system('./time_evolution.pl ../outputs/output_dark_cloud.dat final_dark_cloud_abundances_pre_final_format 0 0 '+specs_string)
##########################################################################################
##########################################################################################
read_extract_specs()
