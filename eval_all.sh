seq_len=(100 200 300 400 500 600 700 800)
folders=("full" "slim" "ukf" "ukfslim")

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
      set ylabel 'error [deg/m]';
      plot 'eval/${folder}/${seq}/r_err.txt' with lines, 'eval/gtsam/${seq}/r_err.txt' with lines;";
                
  done
done

for seq in "${seq_len[@]}"; do
    gnuplot -e "set term png;
    set title 'Combined Translation RPE for ${seq}m';
    set output 'eval/t_err_${seq}_all.png';
    set xlabel 'poses';
    set ylabel 'error [%]';
    plot 'eval/full/${seq}/t_err.txt' with lines title 'full', \
         'eval/slim/${seq}/t_err.txt' with lines title 'slim', \
         'eval/ukf/${seq}/t_err.txt' with lines title 'ukf', \
         'eval/ukfslim/${seq}/t_err.txt' with lines title 'ukfslim', \
         'eval/gtsam/${seq}/t_err.txt' with lines title 'gtsam';"
done

for seq in "${seq_len[@]}"; do
    gnuplot -e "set term png;
    set title 'Combined Rotation RPE for ${seq}m';
    set output 'eval/r_err_${seq}_all.png';
    set xlabel 'poses';
    set ylabel 'error [deg/m]';
    plot 'eval/full/${seq}/r_err.txt' with lines title 'full', \
         'eval/slim/${seq}/r_err.txt' with lines title 'slim', \
         'eval/ukf/${seq}/r_err.txt' with lines title 'ukf', \
         'eval/ukfslim/${seq}/r_err.txt' with lines title 'ukfslim', \
         'eval/gtsam/${seq}/r_err.txt' with lines title 'gtsam';"
done


for seq in "${seq_len[@]}"; do
    echo "Method, Translation RPE [\%], Rotation RPE [deg/m]" > eval/rpe_averages_$seq.csv

    for folder in "${folders[@]}"; do
            avg_t_err=$(awk '{sum += $1} END {print sum/NR}' eval/$folder/$seq/t_err.txt)
            avg_r_err=$(awk '{sum += $1} END {print sum/NR}' eval/$folder/$seq/r_err.txt) 
            echo "$folder, $avg_t_err, $avg_r_err" >> eval/rpe_averages_$seq.csv
    done
    avg_t_err=$(awk '{sum += $1} END {print sum/NR}' eval/gtsam/$seq/t_err.txt)
    avg_r_err=$(awk '{sum += $1} END {print sum/NR}' eval/gtsam/$seq/r_err.txt) 
    echo "gtsam, $avg_t_err, $avg_r_err" >> eval/rpe_averages_$seq.csv
done

echo "\\noindent" > eval/tables.tex
for ((i=0; i<${#seq_len[@]}; i+=2)); do
    seq1=${seq_len[$i]}
    seq2=${seq_len[$i+1]}
    if [ -z "$seq2" ]; then
        break
    fi
    awk -v seqlen="$seq1" 'BEGIN {
        FS=","; 
        print "\\begin{minipage}{.5\\textwidth}";
        print "\\resizebox{\\textwidth}{!}{";
        print "\\begin{tabular}{|c|c|c|}";
        print "\\hline";
        print "\\multicolumn{3}{|c|}{\\textbf{Average RPE@" seqlen " [m]}} \\\\";
        print "\\hline";
    }
    {
        printf "%s & %s & %s \\\\\n", $1, $2, $3;
        print "\\hline";
    } 
    END {
        print "\\end{tabular}}";
        print "\\end{minipage}";
    }' eval/rpe_averages_$seq1.csv >> eval/tables.tex
    awk -v seqlen="$seq2" 'BEGIN {
        FS=","; 
        print "\\begin{minipage}{.5\\textwidth}";
        print "\\resizebox{\\textwidth}{!}{";
        print "\\begin{tabular}{|c|c|c|}";
        print "\\hline";
        print "\\multicolumn{3}{|c|}{\\textbf{Average RPE@" seqlen " [m]}} \\\\";
        print "\\hline";
    }
    {
        printf "%s & %s & %s \\\\\n", $1, $2, $3;
        print "\\hline";
    } 
    END {
        print "\\end{tabular}}";
        print "\\end{minipage} \\\\ \\\\";
    }' eval/rpe_averages_$seq2.csv >> eval/tables.tex
done