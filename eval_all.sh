seq_len=(100 200 300 400 500 600 700 800)
folders=("full" "slim" "ukf" "ukfslim")

do_gtsam=false


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