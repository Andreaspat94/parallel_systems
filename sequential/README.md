# sequential

### Initial version

Average execution in argo:
```
-> 1680, 1680, 0.8, 1, 1e-13, 50
Iterations= 50 Elapsed MPI Wall time is 0.612990
Time taken 0 seconds 769 milliseconds
Residual 6.75761e-10
The error of the iterative solution is 0.000317239
```

### Precalculate "f"

Average execution in argo:
```
-> 1680, 1680, 0.8, 1, 1e-13, 50
Iterations= 50 Elapsed MPI Wall time is 0.280155
Time taken 0 seconds 437 milliseconds
Residual 6.7576e-10
The error of the iterative solution is 0.000317239
```

The problem with pre-calculating the whole ***f*** domain is that it
requires a `n*m`-sized array to store its values.

Your solution which was pre-calculating only the ***fX*** and ***fY***
and then calculate ***f*** inside each loop iteration, yielded about
the same improvement in execution time but with less memory usage
(i.e. `n+m`).
