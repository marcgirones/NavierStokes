Valoration execution command:
python convergence_check.py --sizes 10 50 100 150 200

Main compilation command withh fft:
g++ -std=c++17 -O2 *.cpp -lfftw3 

FFT library
sudo apt install libfftw3-dev

Execution cpp
./a.out N tFinal nu dt

Execution time convergence python (random values)
python3 convergence_time_check.py --exe ./a.out --sizes 16 32 64 --time 0.01 --dt 0.0001 --nu 1.51e-5