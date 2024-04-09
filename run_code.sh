# Save the starting directory
start_dir=$(pwd)
# Define directories
#M_dir=('M_star_0.1' 'M_star_0.5' 'M_star_1.0')
M_dir=('M_star_0.5')
#C_O=('C_O_0.2' 'C_O_0.44' 'C_O_0.9' 'C_O_1.2' 'C_O_1.5')
C_O=('C_O_0.9' 'C_O_1.2' 'C_O_1.5')
sub_grids=('grid_0' 'grid_1' 'grid_2' 'grid_3' 'grid_4')

for i_star in "${!M_dir[@]}"; do
    if [[ -d "${M_dir[i_star]}" ]]; then
        cd "${M_dir[i_star]}" || exit
        for i_CO in "${!C_O[@]}"; do
            if [[ -d "${C_O[i_CO]}" ]]; then
                cd "${C_O[i_CO]}" || exit
                for x in "${sub_grids[@]}"; do
                    if [[ -d "$x/core_chemical_model" ]]; then
                        (
                            cd "$x/core_chemical_model" && {
                                echo "Running in directory: $PWD in background"
                                # Assuming you meant to run a Makefile here. Adjust as needed.
                                ./makefile && ./disk_model file_parameters.txt
                            }
                        ) &
                    else
                        echo "Directory $x/core_chemical_model does not exist."
                    fi
                done
                wait
                echo "All disk_model scripts in ${C_O[i_CO]} have completed."
                # Return to the M_dir directory
                cd "$start_dir/${M_dir[i_star]}" || exit
            else
                echo "Directory ${C_O[i_CO]} does not exist."
            fi
        done
        # Return to the starting directory
        cd "$start_dir" || exit
    else
        echo "Directory ${M_dir[i_star]} does not exist."
    fi
done
