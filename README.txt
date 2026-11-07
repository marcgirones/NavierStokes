FFT library
sudo apt install libfftw3-dev


Build project:
cmake -S . -B build
 - To execute case A:
	cmake --build build --target caseA
 - To execute case B: 
 	cmake --build build --target caseB
 - To execute case C:
 	cmake --build build --target caseC
 - To execute all 3 cases:
 	cmake --build build --target caseAll
