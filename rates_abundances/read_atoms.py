import os

##########################################################################################
##########################################################################################
def read_extract_specs(C_O):
    myfile = open('rate12_C_O_'+str(C_O)+'.specs', 'r')
    lines = myfile.readlines()
    myfile.close()
    specs = []
    for line in lines:
        parts = line.split()
        second_column_value = parts[1]
        specs.append(second_column_value)

    specs_string = ' '.join(specs)
    os.system('./time_evolution.pl ../outputs/output_dark_cloud_C_O_'+str(C_O)+'.dat final_dark_cloud_abundances_pre_final_format_C_O_'+str(C_O)+' 0 0 '+specs_string)
##########################################################################################
##########################################################################################
C_O = [0.2,0.44,0.9,1.2,1.5]
for i in range(len(C_O)):
    read_extract_specs(C_O[i])
