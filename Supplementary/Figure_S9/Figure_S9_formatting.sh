while IFS=$'\t' read -r col1 col2 col3;
do
    sbatch --output=Simulations_${col1}_${col2}_${col3}.out --error=Simulations_${col1}_${col2}_${col3}.err --job-name=Simulations_${col1}_${col2}_${col3} Simulation_test.sh ${col1} ${col2} ${col3}
done < combinations.txt



touch fraction_plots.txt
while read file_name;
do
    amount=$(awk 'NR > 1 && $7 >= 0.45 && $8 <= 0.13' ${file_name} | wc -l)
    total=$(awk 'NR > 1' ${file_name} | wc -l)
    fraction=$(echo "scale=5; ${amount} / ${total}" | bc)
    echo ${fraction} >> fraction_plots.txt
done < grouped_files.txt




