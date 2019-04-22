 set offset 1,1,1,1
 set terminal png
 set output 'plot.png'
 plot 'ports_xy' using 2:3:(sprintf("(%d)", $1)) with labels point  pt 7 offset char 1,1 notitle 

