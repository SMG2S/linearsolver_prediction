make compile
for i in {1..10000}; do /usr/bin/mpiexec -n 6 ./main.exe; done
	
