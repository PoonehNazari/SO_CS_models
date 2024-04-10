import numpy as np
###########################################################################
###########################################################################
def write_specs_file(C_O):
    myfile = open('rate12_C_O_'+str(C_O)+'.specs', 'r')
    lines = myfile.readlines()
    myfile.close()
    n = 0
    for i_line in range(len(lines)):
        parts = lines[i_line].split()
        if parts[1]=='H2*': #This one for some reason was not transferred with time*.pl script (did the check but every time you add a new C/O check if H2* is teh one and if it is zero in in the outputs dir)
            index = n
        n = n+1

    print(index)

    data = np.loadtxt('final_dark_cloud_abundances_pre_final_format_C_O_'+str(C_O)+'.dat', skiprows=8)
    time = data[-1,0] #just interesetd in the final abundnaces nit the time evolution that's why have -1
    #1,2 are H2 and He but I already have them in my list so we ignore here
    abunds = []
    for i in range(len(data[-1])):
        if i > 2:
            abunds.append(data[-1,i])
    abunds = np.array(abunds)
    abunds = np.insert(abunds, index,0.000e+00)
    init_abund_file = open('rate12_init_abunds_C_O_'+str(C_O)+'.specs', 'w')
    if len(lines) == len(abunds):
        print('All good with lengths :)')
        for i_line in range(len(lines)):
            parts = lines[i_line].split()
            parts[2] = abunds[i_line]
            init_abund_file.write(str(parts[0])+'\t'+parts[1]+'\t'+f'{parts[2]:.3e}'+'\t'+str(parts[3])+'\n')

    init_abund_file.close()
###########################################################################
###########################################################################
C_O = [0.2,0.44,0.9,1.2,1.5]
for i in range(len(C_O)):
    write_specs_file(C_O[i])
