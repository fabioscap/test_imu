seq_len=(100 200 300 400 500 600 700 800)
folders=("full" "slim" "ukf" "ukf_slim")

do_gtsam=false
if [ "$do_gtsam" = true ]; then
  echo "running GTSAM"
  mkdir -p "eval/gtsam"
  rosrun test_imu nclt_gtsam > eval/gtsam/log.txt 
  for seq in "${seq_len[@]}"; do
      mkdir -p eval/gtsam/$seq
      rosrun test_imu compute_rpe $seq ./nclt_gtsam.csv ./data/nclt/gt_float.csv
      mv t_err_$seq.txt eval/gtsam/$seq/t_err.txt
      mv r_err_$seq.txt eval/gtsam/$seq/r_err.txt  
  done
fi

for folder in "${folders[@]}"; do
    mkdir -p "eval/$folder"
done

for folder in "${folders[@]}"; do
  echo "running $folder"
  rosrun test_imu nclt $folder > eval/$folder/log.txt # folder is also the preintegration type
  for seq in "${seq_len[@]}"; do
      mkdir -p eval/$folder/$seq
      rosrun test_imu compute_rpe $seq ./nclt.csv ./data/nclt/gt_float.csv
      mv t_err_$seq.txt eval/$folder/$seq/t_err.txt
      mv r_err_$seq.txt eval/$folder/$seq/r_err.txt
      gnuplot -e "set term png;
      set title 'Translation RPE';
      set xlabel 'poses';
      set ylabel 'error [%]';
      set output 'eval/${folder}/${seq}/t_err.png';
      plot 'eval/${folder}/${seq}/t_err.txt' with lines, 'eval/gtsam/${seq}/t_err.txt' with lines;"
      gnuplot -e "set term png;
      set title 'Orientation RPE';
      set output 'eval/${folder}/${seq}/r_err.png';
      set xlabel 'poses';
      set ylabel 'error [deg]';
      plot 'eval/${folder}/${seq}/r_err.txt' with lines, 'eval/gtsam/${seq}/r_err.txt' with lines;";
                
  done
done
