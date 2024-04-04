import numpy as np
###########################################################################
###########################################################################
def write_specs_file():
    myfile = open('rate12_full_atomic_v2.specs', 'r')
    lines = myfile.readlines()
    myfile.close()
    n = 0
    for i_line in range(len(lines)):
        parts = lines[i_line].split()
        if parts[1]=='H2*': #This one for some reason was not transferred with time*.pl script
            index = n
        n = n+1

    print(index)

    data = np.loadtxt('final_dark_cloud_abundances_pre_final_format.dat', skiprows=8)
    time = data[-1,0] #just interesetd in the final abundnaces nit the time evolution that's why have -1
    #1,2 are H2 and He but I already have them in my list so we ignore here
    abunds = []
    for i in range(len(data[-1])):
        if i > 2:
            abunds.append(data[-1,i])
    abunds = np.array(abunds)
    abunds = np.insert(abunds, index,0.000e+00)
    init_abund_file = open('rate12_init_abunds.specs', 'w')
    if len(lines) == len(abunds):
        print('All good with lengths :)')
        for i_line in range(len(lines)):
            parts = lines[i_line].split()
            parts[2] = abunds[i_line]
            init_abund_file.write(str(parts[0])+'\t'+parts[1]+'\t'+f'{parts[2]:.3e}'+'\t'+str(parts[3])+'\n')

    init_abund_file.close()
###########################################################################
###########################################################################
write_specs_file()
