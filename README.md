

# QuadMatrix
## A simple library for matrix calculations with different levels of precision.


**QuadMatrix** is a simple library for matrix calculations on square matrices of real numbers, 
with different levels of precision. It includes the following public classes:

-   **Matrix**  -- a generic abstract class which defines the set of operations that all 
its subclasses must implement;
-   **DoubleMatrix**  -- a subclass of the Matrix that stores the matrix data and implements 
operations on them using standard Java primitive  `double`  type;
-   **BigDecimalMatrix**  -- a subclass of the Matrix that stores the matrix data and implements 
operations on them using the standard `BigDecimal`  class, which allows to gain (theoretically) unlimited precision;
-   **QuadrupleMatrix**  -- a subclass of the Matrix that stores the matrix data and implements 
operations on them using [`Quadruple`](https://m-vokhm.github.io/Quadruple/src/main/javadoc/com/mvohm/quadruple/Quadruple.html "class in com.mvohm.quadruple") class (see  the [Github repository](https://github.com/m-vokhm/Quadruple)), which allows to gain an accuracy of about 36 decimal places *(may vary depending on the matrix size, it's data, and the nature of performed operation*) while performing calculations much faster than  `BigDecimalMatrix`.

***The set of the operations*** which can be performed by the concrete public classes is defined by their common parent,
abstract class `Matrix`. They include:
- creating new instances with data provided as constructor arguments of different types; 
- solving system of linear equations of form  ***Ax = b***; 
- solving matrix equations of form  ***AX = B***; 
- inversion,
- transposition,
- multiplication by a matrix, by a vector and by a scalar,
- matrix addition and subtraction, 
- computing determinant, norm and condition number.
 
For solving equations and inversion, there are additional methods with higher accuracy which use
iterative refinement of the solution.

The current version is ***limited to square matrices*** and ***does not provide any optimization for sparse matrices***.
Although the classes are immutable in the sense that the user can not change their inner states after they are created,
in the current version they are ***not thread-safe***.

For more details, see the [documentation](https://m-vokhm.github.io/QuadMatrix/doc/com/mvohm/quadmatrix/package-summary.html).

#### Motivation

The main goal of the project was to provide an opportunity to solve matrix-dependent tasks with precision higher than can be provided by matrix libraries using `double` arithmetic and to do it not so slow as it is done by libraries using `BigDeciamal`. 

Average times of execution and values of computational errors (`sqrt(MSE)`) for the most common operations on different types of matrices of size 400 x 400 are shown in the diagrams below. The abbreviation `BDMatrix (40)` refers to `BigDecimalMatrix` with precision set to 40, `BDMatrix (80)` means `BigDecimalMatrix` with precision 80. The Y-axes are given in logarithmic scale, since the range of the values is too large to be depicted in linear scale.  The results obtained with the classic JAMA matrix library are added for reference. The times are measured with Bellsoft Liberica JDK 17 on a computer with Intel(R) Core(TM) i5-4460 3200 processor. 
![errors](https://github.com/m-vokhm/QuadMatrix/blob/master/doc/readme-images/errors.png)
![times](https://github.com/m-vokhm/QuadMatrix/blob/master/doc/readme-images/times.png)

More precise values of the data shown above are given in the following table:

    ------------------------------------------------------------------------------------------------
    Times, ms                  Jama.Matrix  DoubleMatrix   QuadrMatrix   BDMatrix_40   BDMatrix_80  
    ------------------------------------------------------------------------------------------------
    Simple vector solution            21.4           9.3       2,293.2      19,342.4      27,833.5 
    Accurate vector solution                        19.7       2,772.1      22,161.8      32,501.9 
    Simple SPD solution               10.6          10.5       1,578.6       9,063.8      13,771.7 
    Accurate SPD solution                           19.5       1,758.9      11,354.9      17,057.5 
    Simple matrix solution            66.7          55.5       9,217.0      75,628.2     108,988.2 
    Accurate matrix solution                     5,135.4     293,045.4   1,322,195.2   1,649,843.5 
    Simple inversion                  67.0          55.5       7,618.2      59,322.9      84,461.2 
    Accurate inversion                           4,676.5     272,901.8   1,162,496.2   1,527,991.2 
    Multiplication                    61.7         370.6      12,819.6      60,059.1      65,332.5 
    ------------------------------------------------------------------------------------------------

    ------------------------------------------------------------------------------------------------
    sqrt(MSE)                   Jama.Matrix  DoubleMatrix   QuadrMatrix   BDMatrix_40   BDMatrix_80                                                                   
    ------------------------------------------------------------------------------------------------
    Simple vector solution         1.1E-13       9.5E-14       1.2E-36       1.1E-36       2.0E-77
    Accurate vector solution                     3.3E-15       1.2E-38       8.9E-40       5.4E-80
    Simple SPD solution            1.2E-13       1.1E-13       1.6E-36       2.6E-37       2.8E-77
    Accurate SPD solution                        1.8E-14       1.7E-37       1.5E-38       1.6E-78
    Simple matrix solution         9.1E-14       4.1E-14       2.3E-36       1.8E-37       3.3E-77
    Accurate matrix solution                     1.4E-15       1.3E-38       6.7E-40       1.1E-79
    Simple inversion               4.4E-14       4.0E-14       4.0E-37       5.3E-38       6.2E-78
    Accurate inversion                           2.8E-16       5.8E-38       8.0E-40       1.8E-79
    Multiplication                 1.2E-18       1.6E-19       1.0E-82       7.9E-42       7.9E-82
    ------------------------------------------------------------------------------------------------
    
More detailed information on performance and errors on different matrix sizes along with a tool to measure them can be found in a separate project, [QuadMatrixMeasurements](https://github.com/m-vokhm/QuadMatrixMeasurements)  

#### Testing
For each non-abstract `Matrix` class, the tests are divided into four groups:
- tests for constructors and data getters (envelope tests);
- tests for vector solutions;
- tests for matrix solutions;
- tests for other operations (inversions, transposition, multiplications and others)  

The names of test classes are formed from the name of the tested class and the group of operations tested with this class, so that `DoubleMatrixEnvelopeTests.java` tests constructors and data getters of `DoubleMatrix`, and `QuadrupleMatrixMatrixSolutionsTests.java` tests matrix solutions performed by `QuadrupleMatrix.java`. The test classes are located in `/src/test/java/com/mvohm/quadmatrix/`. The complete execution of all tests may take up to a few tens of minutes, depending on the computer performance.
