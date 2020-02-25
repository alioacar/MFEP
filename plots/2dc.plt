set size square 
set contour
set encoding iso_8859_1
unset title
set pm3d map
set cont base
set cntrparam levels incremental -4.5,0.35,0
set cntrlabel onecolor format '%8.3g' font 'Times,13' start 2 interval 20
set xlabel "x ({\305})" font "Courier,13"
set ylabel "z ({\305})" font "Courier,13"
set cblabel "PMF (kcal/mol)" font "Courier,13"
set label 1 "MFEP" font "Courier,14" tc rgb "white" front
set label 1 at graph 0.41,0.56
#set label 2 "B" font "Times,12" tc rgb "black" front
#set label 2 at graph 0.51,0.26
#set label 3 "C" font "Times,12" tc rgb "black" front
#set label 3 at graph 0.62,0.42
#set label 4 "D" font "Times,12" tc rgb "black" front
#set label 4 at graph 0.36,0.46
#set label 5 "E_{1}" font "Times,12" tc rgb "yellow" front
#set label 5 at graph 0.33,0.40
#set label 6 "E_{2}" font "Times,12" tc rgb "yellow" front
#set label 6 at graph 0.57,0.20
#set label 7 "G" font "Times,12" tc rgb "yellow" front
#set label 7 at graph 0.95,0.74
#set label 8 "H" font "Times,12" tc rgb "yellow" front
#set label 8 at graph 0.62,0.74
#set label 9 "F" font "Times,12" tc rgb "yellow" front
#set label 9 at graph 0.08,0.30
set mxtics 5
set mytics 5
set cbtics 1
set cbrange [-4.3:0]
set xtics offset 0,0.8,1 font "Courier,12"
set ytics offset 0,0,1 font "Courier,12"
set cbtics font "Courier,11"
set xlabel offset 0.9,1
set ylabel offset -1.2,1
set cblabel offset 2,1
set style line 100 lt 0.5 lc rgb "black" lw 0.6
set style line 101 lt 0.5 lc rgb "black" lw 0.6
#set yrange [25.8:26.4] 
#set xrange [2:2.8] 
set yrange [20.0:26.5] 
set xrange [-5.5:2.0] 
set palette rgbformulae 23,3,28
set parametric
set style line 2 lc rgb 'black' pt 1 
set style line 3 lc rgb 'white' pt 7 
set style line 4 lc rgb 'yellow' pt 5 
set style line 5 lc rgb 'blue' pt 4 
set style line 6 lc rgb 'violet' pt 6 
splot "fort.17" with pm3d,\
"fort.53" using 1:2:(0) with lines ls 3 notitle
#"backward1" using 1:2:(0) with lines ls 4 notitle,"backward2" using 1:2:3 with points ls 4
#"forward1" using 1:2:(0) with lines ls 4 notitle,"forward2" using 1:2:3 with points ls 4
#splot "fort.17" with pm3d,"fort.32" using 1:2:3 with points ls 2,"fort.33" using 1:2:3 with points ls 3,\
#splot "fort.17" with pm3d,"GtoHstate-pathway" using 1:2:(0) with lines ls 3 notitle,\
#"HtoE1state-pathway" using 1:2:(0) with lines ls 3 notitle,\
#"fort.53" using 1:2:(0) with lines ls 3 notitle,\
#"E1toFstate-pathway" using 1:2:(0) with lines ls 3 notitle,\
