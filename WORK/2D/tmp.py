#!/usr/bin/python
import os

for i in range(0,50):
   print "data/Mach_%2.2i.d" % (i)
   os.system("""
gnuplot <<EOF
set grid
set xlabel "x"
set ylabel "Mach Number"
set xrange [0:1]
set terminal png
set key below
set output "png/Mach_%2.2i.png"
p "data/Mach_%2.2i.d" w lp
EOF""" % (i,i))
   print "data/density_%2.2i.d" % (i)
   os.system("""
gnuplot <<EOF
set grid
set xlabel "x"
set ylabel "Density"
set xrange [0:1]
set key below
set terminal png
set output "png/density_%2.2i.png"
p "data/density_%2.2i.d" w lp
EOF""" % (i,i))

   print "data/pressure_%2.2i.d" % (i)
   os.system("""
gnuplot <<EOF
set grid
set xlabel "x"
set ylabel "Pressure"
set xrange [0:1]
set terminal png
set key below
set output "png/pressure_%2.2i.png"
p "data/pressure_%2.2i.d" w lp
EOF""" % (i,i))
   print "data/U_Velocity_%2.2i.d" % (i)
   os.system("""
gnuplot <<EOF
set grid
set xlabel "x"
set ylabel "U_Velocity"
set xrange [0:1]
set terminal png
set key below
set output "png/U_Velocity_%2.2i.png"
p "data/U_Velocity_%2.2i.d" w lp
EOF""" % (i,i))
   print "data/V_Velocity_%2.2i.d" % (i)
   os.system("""
gnuplot <<EOF
set grid
set xlabel "x"
set ylabel "V_Velocity"
set xrange [0:1]
set terminal png
set key below
set output "png/V_Velocity_%2.2i.png"
p "data/V_Velocity_%2.2i.d" w lp
EOF""" % (i,i))
#os.system("animate -delay 100 png/Mach*.png")
#os.system("animate -delay 100 png/Pressure*.png")
#os.system("animate -delay 100 png/Density*.png")
#os.system("animate -delay 100 png/Velocity*.png")
