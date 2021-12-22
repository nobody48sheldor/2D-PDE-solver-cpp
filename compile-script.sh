echo 
echo compiling...
echo 
g++ solver.cpp -o solver-cpp
echo 
echo running...
echo 
./solver-cpp

echo Do you want to keep the renders ? [Y/n]  
read rep

if $rep == "N" or "n":
	rm renders/*.png
else:
	echo ok !
