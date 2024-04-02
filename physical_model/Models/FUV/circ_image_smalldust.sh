M_dir=('M_star_0.1' 'M_star_0.5' 'M_star_1.0')

for x in "${M_dir[@]}";
do
    (
    cd "$x" && {
        echo "running in directory: $PWD/ in background"
        python3 setup_model.py
}
    ) &
done
wait
echo "All setup_model.py scripts in ${M_dir[*]} have completed."
