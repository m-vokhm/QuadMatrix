

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
operations on them using 
[`Quadruple`](https://m-vokhm.github.io/Quadruple/src/main/javadoc/com/mvohm/quadruple/Quadruple.html "class in com.mvohm.quadruple")
class (see  the [Github repository](https://github.com/m-vokhm/Quadruple)), which allows to gain an accuracy of about 36 decimal 
places (may vary depending on the matrix size, it's data, and the nature of performed operation) and performs calculations
much faster than  `BigDecimalMatrix`.

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

### Motivation

The main goal of the project was to provide an opportunity to solve matrix-dependent tasks with precision 
higher than can be provided by matrix libraries using `double` arithmetic and do it not so slow
as it is done by libraries using `BigDeciamal`. 

Approximate characteristics of different types of matrices for different types of operations and different
sizes are shown in the diagrams below (***coming soon***).