# Collisions

Recommended use case detailed below. Visualising the animations required gnuplot.

Compile as preferred, I use g++ compiler.
```
g++ collision.cc -o collision
```

Run the program, pointing the output to data.txt
```
./collision > data.txt
```

Open gnuplot, and load the required scripts.
```
gnuplot

load('plotting/plot16.txt')
```

Other plotting scripts are available depending on the selected output from the code. Instructive methods to try are increasing the total number of collisions and plotting the average energy using 'energy.txt'
