echo 
echo compiling...
echo 
g++ conduction.cpp -o conduction-cpp
echo 
echo running...
echo 
./conduction-cpp
echo plotting
python3 plot-script.py

echo Do you want to keep the renders ? [Y/n]  
read rep

if rep == "N" or "n":
	rm renders/*.png
else:
	echo ok !
