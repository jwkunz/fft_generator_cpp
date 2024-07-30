This is a python script that will generate a C++ function to compute a power-of-2 FFT using the split radix architecture.

The point of this is that the C++ function is completely flattend with no for loops or recursion.  
It just describes the mathmatical operations and order they need to occur.  The compiler can be free to
optimize the mathmatical operations in the function as it sees fit.

Warning:  The resulting file tends to be huge and takes significant compilation time.
