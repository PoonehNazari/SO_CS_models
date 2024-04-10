import os

##########################################################################################
##########################################################################################
def read_extract_specs(C_O):
    """
    A function that reads the output*.dat from the Fortran models of the dark cloud runs and then writes a new
    input rate12_*specs file for th epost-dark cloud chemistry
    Takes:
        C_O : the assumed C/O ratios
    Returns:
        None
    """
    #Read the input .specs file of the dark cloud run to find the elements and molecules that we need to input
    homedir = '/Users/pooneh/Library/Mobile Documents/com~apple~CloudDocs/Academic/ESO/projects/SO_CS_models/'
    myfile = open('rate12_C_O_'+str(C_O)+'.specs', 'r')
    lines_input = myfile.readlines()
    myfile.close()
    specs = []
    for line in lines_input:
        parts = line.split()
        second_column_value = parts[1]
        specs.append(second_column_value)

    #Read the output file and find the abundances of the species we picked from previous loop at the final time
    output = open(homedir+'outputs/output_dark_cloud_C_O_'+str(C_O)+'.dat')
    lines_output = output.readlines()[16:]
    output.close()
    dict_output = {}
    for species in specs:
        for line in lines_output:
            if line.startswith(' '+species+' '):
                parts = line.split()
                value = [float(value) for value in parts[1:]][-1]
                dict_output[species] = {'abund':value}

    #write the final file of *.specs to input in the protostellar chemistry :)
    init_abund_file = open('rate12_init_abunds_C_O_'+str(C_O)+'.specs', 'w')
    for i_line in range(len(lines_input)):
        parts = lines_input[i_line].split()
        parts[2] = dict_output[parts[1]]['abund']
        init_abund_file.write(str(parts[0])+'\t'+parts[1]+'\t'+f'{parts[2]:.3e}'+'\t'+str(parts[3])+'\n')

    init_abund_file.close()

    #specs_string = ' '.join(specs)
    #os.system('./time_evolution.pl ../outputs/output_dark_cloud_C_O_'+str(C_O)+'.dat final_dark_cloud_abundances_pre_final_format_C_O_'+str(C_O)+' 0 0 '+specs_string)
##########################################################################################
##########################################################################################
C_O = [0.2,0.44,0.9,1.2,1.5]
for i in range(len(C_O)):
    read_extract_specs(C_O[i])
