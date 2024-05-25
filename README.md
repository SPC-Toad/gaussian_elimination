# Gaussian Elimination Calculator

This project is a Python3 implementation of a Gaussian Elimination calculator using NumPy. The calculator converts a given matrix into its reduced echelon form using Gaussian Elimination, a method used in linear algebra for solving systems of linear equations.

## Features

- **Leftmost Nonzero Index**: Finds the leftmost nonzero entry in a matrix after a given row.
- **Inconsistent Row Check**: Checks if a row represents an inconsistent equation.
- **Leading Entry Index**: Finds the index of the leading entry of a row.
- **Row Swapping**: Swaps the rows of a matrix.
- **Zeroing in Pivot Column**: Subtracts a row with a pivot position from another row to zero out a position at the pivot column index.
- **Scaling to One in Pivot Column**: Scales a pivot position to 1.
- **Elimination Phase**: Converts a matrix into echelon form.
- **Back Substitution Phase**: Converts a matrix in echelon form into reduced echelon form.
- **Gaussian Elimination**: Combines elimination and back substitution phases to convert a matrix into reduced echelon form.

## Requirements

- Python 3.x
- NumPy library

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/gaussian-elimination.git
    ```

2. Navigate to the project directory:
    ```bash
    cd gaussian-elimination
    ```

3. Install the required dependencies:
    ```bash
    pip install numpy
    ```

## Usage

1. To run the Gaussian Elimination calculator, use the following command:
    ```bash
    make run
    ```

2. The script will perform Gaussian Elimination on the predefined test cases and display the results before and after the elimination process.

## Example

Here is an example of how to use the Gaussian Elimination calculator:

```python
import numpy as np

a = np.array([[1., 2, 3, 4, 5], 
              [5., 6, 7, 8, 5], 
              [9., 10, 11, 12, 5], 
              [13., 14, 15, 16, 5]])
print("before elimination:")
print(a)
gaussian_elimination(a)
print()
print("after elimination:")
print(a)

b = np.array([[3., -4, 3, -9], [6, 7, -3, 0], [1, 0, 10, -21]])
print("before elimination:")
print(b)
gaussian_elimination(b)
print()
print("after elimination:")
print(b)
```

## Functions

- `leftmost_nonzero_index(a, row_index)`: Finds the leftmost nonzero entry in a matrix after a given row.
- `is_inconsistent_row(row)`: Checks if a row represents an inconsistent equation.
- `leading_entry_index(row)`: Finds the index of the leading entry of a row.
- `num_of_rows(a)`: Returns the number of rows in a matrix.
- `swap_rows(a, row_index_1, row_index_2)`: Swaps the rows of a matrix.
- `zero_in_pivot_column(a, target_row_index, pivot_row_index, pivot_col_index)`: Zeroes out a position at the pivot column index.
- `scale_to_one_in_pivot_column(a, pivot_row_index, pivot_col_index)`: Scales a pivot position to 1.
- `elimination_phase(a)`: Converts a matrix into echelon form.
- `back_substitution_phase(a)`: Converts a matrix in echelon form into reduced echelon form.
- `gaussian_elimination(a)`: Converts a matrix into reduced echelon form.

## Makefile

The provided `Makefile` contains a simple command to run the Gaussian Elimination calculator:
```sh
    make
```

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue to discuss any improvements or bugs.
