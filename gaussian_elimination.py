import numpy as np

def leftmost_nonzero_index(a, row_index):
    """finds the leftmost nonzero entry after a given row

    parameters:

    a: 2D numpy array
    row_index: index for a row of A

    returns:

    (i, j) such that A[i, j] is the leftmost nonzero entry in A at or
    below the row ROW_INDEX, or None if no such value exists

    """
    nonzeros = np.transpose(np.nonzero(a[row_index:]))
    nonzeros.sort(key=lambda x: x[1])
    if len(nonzeros) != 0:
        l_r, l_c = nonzeros[0]
        return l_r + row_index, l_c
    return None

def is_inconsistent_row(row):
    """checks if a row represents an inconsistent equation

    paramaters:

    row: 1D numpy array

    returns:

    True if ROW is of the form [0 0 0 ... a] where a != 0, and False
    otherwise.
    """
    # TODO: FILL IN THIS FUNCTION AND CHANGE THE RETURN VALUE
    if (np.all(row[:-1] == 0) and row[-1] != 0) : 
        return True; # check if any values of the row are zeros except last index.
    else: 
        return False

def leading_entry_index(row):
    """finds the index of the leading entry of a row

    parameters:

    row: 1D numpy array

    returns:

    an index i such that ROW[i] is the leading entry (the first
    nonzero value) of ROW

    """
    nonzeros = np.nonzero(row)[0]
    if len(nonzeros) != 0:
        return nonzeros[0]
    return None

def num_of_rows(a):
    """the number of rows if a matrix

    parameters:

    a: 2D numpy array

    returns:

    the number of rows of A

    """
    return a.shape[0]

def swap_rows(a, row_index_1, row_index_2):
    """swaps the rows of a matrix

    paramters:

    a: 2D numpy array
    row_index_1: index for a row of A
    row_index_2: index for a row of A

    returns:

    None

    A is mutated so that row ROW_INDEX_1 and row ROW_INDEX_2 are
    swapped

    """
    a[[row_index_1, row_index_2]] = a[[row_index_2, row_index_1]]

def zero_in_pivot_column(a, target_row_index, pivot_row_index, pivot_col_index):
    """subtracts a row with a pivot position from another row to
    zero-out a position at the pivot column index

    paramters:

    a: 2D numpy array
    target_row_index: index of a row of A
    pivot_row_index: row index of a pivot position of A
    pivot_col_index: column index of a pivot position of A

    returns:

    True if A[TARGET_ROW_INDEX] becomes inconsistent after mutation,
    and false otherwise

    A is mutated by

        A[TARGET_ROW_INDEX] -= c * A[PIVOT_ROW_INDEX]

    where c is chosen so that A[TARGET_ROW_INDEX, PIVOT_COL_INDEX] == 0.0

    sources:

    based on an implementation from the Jupyter Guide to Linear Algebra
    https://bvanderlei.github.io/jupyter-guide-to-linear-algebra/intro.html

    """
    tamp_row = np.vectorize(lambda x: 0.0 if np.isclose(x, 0.0) else x)
    pivot_val = a[pivot_row_index, pivot_col_index]
    elim_val = a[target_row_index, pivot_col_index]
    a[target_row_index] -= a[pivot_row_index] * elim_val / pivot_val
    a[target_row_index] = tamp_row(a[target_row_index])
    a[target_row_index, pivot_col_index] = 0.0
    return is_inconsistent_row(a[target_row_index])

def scale_to_one_in_pivot_column(a, pivot_row_index, pivot_col_index):
    """scales a pivot position to 1

    parameters:

    a: 2D numpy array
    pivot_row_index: row index of a pivot position
    pivot_col_index: column index of a pivot position

    returns:
    None

    A is mutated by

        A[PIVOT_ROW_INDEX] /= c

    where c is the entry at the pivot position of A

    """
    #print(a, pivot_row_index, pivot_col_index)
    pivot_value = a[pivot_row_index, pivot_col_index]
    
    a[pivot_row_index] /= pivot_value

def elimination_phase(a):
    """converts a matrix into echelon form

    paramters:

    a: 2D numpy array

    returns:

    True if A represents an CONSISTENT system, and False otherwise

    A is mutated so that if A is consistent, it is converted into
    echelon form

    NOTE: this should short circuit and return False if a row
    operation ever creates an inconsistent row

    """
    # TODO: FILL IN THIS FUNCTION AND CHANGE THE RETURN VALUE
    # n is the row length.
    n = num_of_rows(a)
    x = n - 1
    
    if (leading_entry_index(a[0])) : 
        for i in range(1, n) :
            if leading_entry_index(a[i]) is not None :
                swap_rows(a, 0, i)
                break
    
    for i in range(n) : 
        for row in range(x) :
            if leading_entry_index(a[row]) is not None and leading_entry_index(a[row + 1]) is not None : # check if there is pivot column, if not then go next.
                if (leading_entry_index(a[row]) == leading_entry_index(a[row + 1])) : 
                    zero_in_pivot_column(a, row + 1, row, leading_entry_index(a[row]))
                    scale_to_one_in_pivot_column(a, row, leading_entry_index(a[row]))
                elif (leading_entry_index(a[row]) > leading_entry_index(a[row + 1])) :
                    scale_to_one_in_pivot_column(a, row, leading_entry_index(a[row]))
                    scale_to_one_in_pivot_column(a, row + 1, leading_entry_index(a[row + 1]))
                    swap_rows(a, row + 1, row)
    if (leading_entry_index(a[n-1]) is not None) : 
        scale_to_one_in_pivot_column(a, n - 1, leading_entry_index(a[n - 1]));
        
    return True

def back_substitution_phase(a):
    """converts a matrix in echlon form into one in reduced echelon
    form

    parameters:
    a: 2D numpy array

    returns:

    None

    A is mutated so that if A is in echelon form, it is converted into
    reduced echelon form

    """
    # TODO: FILL IN THIS FUNCTION
    n = num_of_rows(a)
    # You do not want to start from LAST row becasuse that is constant value column.
    for i in range(n-1, 0, -1) : 
        if (leading_entry_index(a[i]) is not None) :
        # you want to remove all the other values on the same pivot column
            for r in range(i) :
                if (leading_entry_index(a[r]) is not None) :
                    zero_in_pivot_column(a, r, i, leading_entry_index(a[i]))
    

def gaussian_elimination(a):
    """converts a matrix into reduced echelon form

    parameters:

    a: 2D numpy array

    returns:

    None

    A is mutated so that is converted to reduced echelon form

    Note: back substitution should only be performed if A is
    consistent

    """
    if elimination_phase(a):
        back_substitution_phase(a)

# # small test case
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

# # # larger test case
# # matrix in reduced echelon form
# Y = np.array(
#     [[1., 0, 0, 0, 2, 3,  0, 0, -9.1],
#      [0,  1, 0, 0, 0, 5.32,  0, 0, 24.3],
#      [0,  0, 1, 0, 1, 11.01, 0, 0, 6],
#      [0,  0, 0, 1, 8, 12, 0, 0, 1],
#      [0,  0, 0, 0, 0, 0,  1, 0, 16],
#      [0,  0, 0, 0, 0, 0,  0, 1, 12]])

# print("before elimination:")
# print(Y)
# gaussian_elimination(Y)
# print()
# print("after elimination:")
# print(Y)

# X = np.copy(Y)

# # a bunch of random row operations
# X[0] += 3 * X[1] + 2 * X[2] + 8 * X[4] + -3 * X[5]
# X[1] += -6 * X[0] + X[2] + 3 * X[3] + 4 * X[5]
# X[2] += X[0] + X[1] + X[3] + X[4]
# X[3] += -X[0] + -X[1] + X[2]
# X[4] += X[5] + 3 * X[2]
# X[5] += X[1] + X[3]
# X[0] *= 3
# X[1] /= 3
# X[2] -= X[4] + 0.5 * X[3]

# print("before elimination:")
# np.round(X, decimals=2)
# print(X)
# gaussian_elimination(X)
# np.round(X, decimals=2)
# print()
# print("after elimination:")
# print(X)
# print()
# if np.isclose(X, Y).all():
#     print("got back the same matrix you started with!")
# else:
#     print("something happened, got the wrong matrix...")
