Run the following commands:
$ mkdir build && cd build && cmake .. && make && mpirun -n 2 ./parallel_systems

Clion editor configuration for MPI (linux):
![image](https://user-images.githubusercontent.com/51324239/188848801-9861b1e1-24ec-486d-a8ae-3c19d819cfab.png)



Sequential program (last commit):

Iterations= 50 Elapsed MPI Wall time is 0.235326
Time taken 0 seconds 440 milliseconds
Residual 5.42574e-09
The error of the iterative solution is 0.000634165
[argo354@argo seqJacobi]$ 


OpenMP + MPI(last commit on 21/9):

Δοκιμασμενο σε 840*840 (μονο τοπικα)


Με το openMP+ MPI καταφέρνω να ριξω το mpi time κατα λίγο, αλλα ανέβηκε ο συνολικος χρόνος..πως γινεται αυτο; δεν εχω καταλαβει ακομα.
Υπάρχουν καποια σχολια για τη δικη μας κατανόηση στο κώδικα του #pragma omp
