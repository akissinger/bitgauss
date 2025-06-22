"""
MyPy type stubs for bitgauss module.

This file provides type annotations for the BitMatrix Python bindings
created with PyO3.
"""

from typing import Callable, List, Optional, Tuple

class BitMatrix:
    """
    A binary matrix implementation optimized for bit operations.
    
    BitMatrix provides efficient storage and operations for matrices containing
    only boolean values (0 or 1), using bit-level operations for performance.
    """
    
    def __init__(self, rows: int, cols: int) -> None:
        """
        Create a new BitMatrix of size rows x cols with all bits set to 0.
        
        Args:
            rows: Number of rows in the matrix
            cols: Number of columns in the matrix
        """
        ...
    
    def bit(self, i: int, j: int) -> bool:
        """
        Get the bit at position (i, j).
        
        Args:
            i: Row index
            j: Column index
            
        Returns:
            The boolean value at position (i, j)
            
        Raises:
            ValueError: If indices are out of bounds
        """
        ...
    
    def set_bit(self, i: int, j: int, b: bool) -> None:
        """
        Set the bit at position (i, j) to the given value.
        
        Args:
            i: Row index
            j: Column index
            b: Boolean value to set
            
        Raises:
            ValueError: If indices are out of bounds
        """
        ...
    
    @staticmethod
    def build(rows: int, cols: int, func: Callable[[int, int], bool]) -> 'BitMatrix':
        """
        Build a BitMatrix from a function that determines the value of each bit.
        
        Args:
            rows: Number of rows in the matrix
            cols: Number of columns in the matrix
            func: Function that takes (row, col) and returns a boolean value
            
        Returns:
            New BitMatrix with values determined by the function
        """
        ...
    
    @staticmethod
    def zeros(rows: int, cols: int) -> 'BitMatrix':
        """
        Create a new BitMatrix of size rows x cols with all bits set to 0.
        
        Args:
            rows: Number of rows in the matrix
            cols: Number of columns in the matrix
            
        Returns:
            New zero matrix
        """
        ...
    
    def is_zero(self) -> bool:
        """
        Check if the matrix consists of all zero bits.
        
        Returns:
            True if all bits are zero, False otherwise
        """
        ...
    
    @staticmethod
    def identity(size: int) -> 'BitMatrix':
        """
        Create a new identity BitMatrix of size size x size.
        
        Args:
            size: Size of the square identity matrix
            
        Returns:
            New identity matrix
        """
        ...
    
    @staticmethod
    def random(rows: int, cols: int, seed: Optional[int] = None) -> 'BitMatrix':
        """
        Create a new random BitMatrix of size rows x cols.
        
        Args:
            rows: Number of rows in the matrix
            cols: Number of columns in the matrix
            seed: Optional random seed for reproducible results
            
        Returns:
            New random matrix
        """
        ...
    
    @staticmethod
    def random_invertible(size: int, seed: Optional[int] = None) -> 'BitMatrix':
        """
        Create a new random invertible BitMatrix of size size x size.
        
        Args:
            size: Size of the square matrix
            seed: Optional random seed for reproducible results
            
        Returns:
            New random invertible matrix
        """
        ...
    
    @property
    def rows(self) -> int:
        """Number of logical rows in the matrix."""
        ...
    
    @property
    def cols(self) -> int:
        """Number of logical columns in the matrix."""
        ...
    
    @property
    def shape(self) -> Tuple[int, int]:
        """Shape of the matrix as a tuple (rows, cols)."""
        ...
    
    def transposed(self) -> 'BitMatrix':
        """
        Return a transposed copy of the matrix.
        
        Returns:
            New transposed matrix
        """
        ...
    
    def transpose_inplace(self) -> None:
        """Transpose the matrix in place."""
        ...
    
    def gauss(self, full: bool = False) -> None:
        """
        Perform gaussian elimination.
        
        Args:
            full: If True, perform full Gauss-Jordan to produce reduced 
                  echelon form, otherwise just echelon form
        """
        ...
    
    def rank(self) -> int:
        """
        Compute the rank of the matrix using gaussian elimination.
        
        Returns:
            The rank of the matrix
        """
        ...
    
    def inverse(self) -> 'BitMatrix':
        """
        Compute the inverse of an invertible matrix.
        
        Returns:
            The inverse matrix
            
        Raises:
            ValueError: If matrix is not square or not invertible
        """
        ...
    
    def vstack(self, other: 'BitMatrix') -> 'BitMatrix':
        """
        Vertically stack this matrix with another one.
        
        Args:
            other: Matrix to stack below this one
            
        Returns:
            New vertically stacked matrix
            
        Raises:
            ValueError: If matrices have different number of columns
        """
        ...
    
    def hstack(self, other: 'BitMatrix') -> 'BitMatrix':
        """
        Horizontally stack this matrix with another one.
        
        Args:
            other: Matrix to stack to the right of this one
            
        Returns:
            New horizontally stacked matrix
            
        Raises:
            ValueError: If matrices have different number of rows
        """
        ...
    
    def nullspace(self) -> List['BitMatrix']:
        """
        Compute a basis for the nullspace of the matrix.
        
        Returns:
            List of BitMatrix instances representing basis vectors
        """
        ...
    
    def copy(self) -> 'BitMatrix':
        """
        Return a copy of the matrix.
        
        Returns:
            New matrix that is a copy of this one
        """
        ...
    
    def __str__(self) -> str:
        """String representation of the matrix."""
        ...
    
    def __repr__(self) -> str:
        """Python representation of the matrix."""
        ...
    
    def __len__(self) -> int:
        """Return the number of rows (for len() function)."""
        ...
    
    def __getitem__(self, key: Tuple[int, int]) -> bool:
        """
        Support for indexing with [i, j]
        
        Args:
            key: A (row, col) tuple for single bit access
                 
        Returns:
            Boolean value at the specified position
            
        Raises:
            ValueError: If indices are out of bounds or invalid type
        """
        ...
    
    def __setitem__(self, key: Tuple[int, int], value: bool) -> None:
        """
        Support for item assignment with [i, j] = value.
        
        Args:
            key: (row, col) tuple for bit position
            value: Boolean value to set
            
        Raises:
            ValueError: If indices are out of bounds or invalid type
        """
        ...
    
    def to_list(self) -> List[List[bool]]:
        """
        Convert matrix to a list of lists.
        
        Returns:
            List of lists representation of the matrix
        """
        ...
    
    @staticmethod
    def from_list(data: List[List[bool]]) -> 'BitMatrix':
        """
        Create matrix from a list of lists.
        
        Args:
            data: List of lists containing boolean values
            
        Returns:
            New BitMatrix created from the data
            
        Raises:
            ValueError: If rows have inconsistent lengths
        """
        ...
    
    def __matmul__(self, other: 'BitMatrix') -> 'BitMatrix':
        """
        Matrix multiplication using the @ operator.
        
        Args:
            other: Right-hand side matrix
            
        Returns:
            Result of matrix multiplication
            
        Raises:
            ValueError: If matrix dimensions are incompatible
        """
        ...
    
    def __mul__(self, other: 'BitMatrix') -> 'BitMatrix':
        """
        Matrix multiplication using the * operator.
        
        Args:
            other: Right-hand side matrix
            
        Returns:
            Result of matrix multiplication
            
        Raises:
            ValueError: If matrix dimensions are incompatible
        """
        ...
    
    def __rmul__(self, other: 'BitMatrix') -> 'BitMatrix':
        """
        Right-hand matrix multiplication.
        
        Args:
            other: Left-hand side matrix
            
        Returns:
            Result of matrix multiplication
            
        Raises:
            ValueError: If matrix dimensions are incompatible
        """
        ...
    
    def __imatmul__(self, other: 'BitMatrix') -> None:
        """
        In-place matrix multiplication using @=.
        
        Args:
            other: Right-hand side matrix
            
        Raises:
            ValueError: If matrix dimensions are incompatible
        """
        ...
    
    def matmul(self, other: 'BitMatrix') -> 'BitMatrix':
        """
        Matrix multiplication method (alternative to operators).
        
        Args:
            other: Right-hand side matrix
            
        Returns:
            Result of matrix multiplication
            
        Raises:
            ValueError: If matrix dimensions are incompatible
        """
        ...
    
    def __pow__(self, exponent: int) -> 'BitMatrix':
        """
        Matrix power (repeated multiplication).
        
        Args:
            exponent: Non-negative integer exponent
            
        Returns:
            Matrix raised to the given power
            
        Raises:
            ValueError: If matrix is not square or exponent is invalid
        """
        ...
    