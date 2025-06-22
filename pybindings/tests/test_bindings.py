import pytest
import bitgauss
from bitgauss import BitMatrix

def test_available():
    # Test if the bitgauss module can be imported and used
    assert bitgauss is not None, "bitgauss module should be importable"

class TestBitMatrixConstruction:
    """Test BitMatrix construction methods."""
    
    def test_new_basic(self):
        """Test basic constructor."""
        matrix = BitMatrix(3, 4)
        assert matrix.rows == 3
        assert matrix.cols == 4
        assert matrix.is_zero()
    
    def test_zeros(self):
        """Test zeros constructor."""
        matrix = BitMatrix.zeros(2, 3)
        assert matrix.rows == 2
        assert matrix.cols == 3
        assert matrix.is_zero()
    
    def test_identity(self):
        """Test identity matrix construction."""
        matrix = BitMatrix.identity(3)
        assert matrix.rows == 3
        assert matrix.cols == 3
        
        # Check diagonal elements are 1
        for i in range(3):
            assert matrix.bit(i, i) == True
            
        # Check off-diagonal elements are 0
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert matrix.bit(i, j) == False
    
    def test_random(self):
        """Test random matrix construction."""
        matrix = BitMatrix.random(3, 4)
        assert matrix.rows == 3
        assert matrix.cols == 4
        
        # Test with seed for reproducibility
        matrix1 = BitMatrix.random(3, 3, seed=42)
        matrix2 = BitMatrix.random(3, 3, seed=42)
        assert matrix1 == matrix2
        
        # Different seeds should produce different matrices (with high probability)
        matrix3 = BitMatrix.random(3, 3, seed=43)
        assert matrix1 != matrix3  # This could theoretically fail but very unlikely
    
    def test_random_invertible(self):
        """Test random invertible matrix construction."""
        matrix = BitMatrix.random_invertible(3)
        assert matrix.rows == 3
        assert matrix.cols == 3
        
        # Should not raise an error when computing inverse
        inverse = matrix.inverse()
        assert inverse.rows == 3
        assert inverse.cols == 3
        
        # Test with seed
        matrix1 = BitMatrix.random_invertible(3, seed=42)
        matrix2 = BitMatrix.random_invertible(3, seed=42)
        assert matrix1 == matrix2
    
    def test_build_with_function(self):
        """Test build method with function."""
        # Create a checker pattern
        def checker(i, j):
            return (i + j) % 2 == 0
        
        matrix = BitMatrix.build(3, 3, checker)
        assert matrix.rows == 3
        assert matrix.cols == 3
        
        # Verify the pattern
        for i in range(3):
            for j in range(3):
                expected = (i + j) % 2 == 0
                assert matrix.bit(i, j) == expected
    
    def test_from_list(self):
        """Test construction from list of lists."""
        data = [
            [True, False, True],
            [False, True, False],
            [True, True, False]
        ]
        matrix = BitMatrix.from_list(data)
        assert matrix.rows == 3
        assert matrix.cols == 3
        
        for i in range(3):
            for j in range(3):
                assert matrix.bit(i, j) == data[i][j]
    
    def test_from_list_empty(self):
        """Test construction from empty list."""
        matrix = BitMatrix.from_list([])
        assert matrix.rows == 0
        assert matrix.cols == 0
    
    def test_from_list_inconsistent_rows(self):
        """Test error handling for inconsistent row lengths."""
        data = [
            [True, False],
            [True, False, True]  # Different length
        ]
        with pytest.raises(ValueError):
            BitMatrix.from_list(data)


class TestBitMatrixAccess:
    """Test bit access and modification."""
    
    def test_bit_access(self):
        """Test getting bits."""
        matrix = BitMatrix.identity(3)
        assert matrix.bit(0, 0) == True
        assert matrix.bit(0, 1) == False
        assert matrix.bit(1, 1) == True
    
    def test_bit_access_out_of_bounds(self):
        """Test bounds checking for bit access."""
        matrix = BitMatrix(2, 3)
        
        with pytest.raises(ValueError, match="Index out of bounds"):
            matrix.bit(2, 0)
        
        with pytest.raises(ValueError, match="Index out of bounds"):
            matrix.bit(0, 3)
    
    def test_set_bit(self):
        """Test setting bits."""
        matrix = BitMatrix.zeros(2, 2)
        matrix.set_bit(0, 1, True)
        matrix.set_bit(1, 0, True)
        
        assert matrix.bit(0, 0) == False
        assert matrix.bit(0, 1) == True
        assert matrix.bit(1, 0) == True
        assert matrix.bit(1, 1) == False
    
    def test_set_bit_out_of_bounds(self):
        """Test bounds checking for bit setting."""
        matrix = BitMatrix(2, 3)
        
        with pytest.raises(ValueError, match="Index out of bounds"):
            matrix.set_bit(2, 0, True)
        
        with pytest.raises(ValueError, match="Index out of bounds"):
            matrix.set_bit(0, 3, True)
    
    def test_indexing_getitem(self):
        """Test indexing with []."""
        matrix = BitMatrix.identity(3)
        assert matrix[0, 0] == True
        assert matrix[0, 1] == False
        assert matrix[1, 1] == True
    
    def test_indexing_getitem_out_of_bounds(self):
        """Test bounds checking for indexing."""
        matrix = BitMatrix(2, 3)
        
        with pytest.raises(ValueError, match="Index out of bounds"):
            _ = matrix[2, 0]
    
    def test_indexing_setitem(self):
        """Test item assignment with []."""
        matrix = BitMatrix.zeros(2, 2)
        matrix[0, 1] = True
        matrix[1, 0] = 1  # Should be truthy
        matrix[1, 1] = 0  # Should be falsy
        
        assert matrix[0, 0] == False
        assert matrix[0, 1] == True
        assert matrix[1, 0] == True
        assert matrix[1, 1] == False
    
    def test_indexing_setitem_out_of_bounds(self):
        """Test bounds checking for item assignment."""
        matrix = BitMatrix(2, 3)
        
        with pytest.raises(ValueError, match="Index out of bounds"):
            matrix[2, 0] = True


class TestBitMatrixOperations:
    """Test matrix operations."""
    
    def test_transpose(self):
        """Test matrix transpose."""
        data = [
            [True, False, True],
            [False, True, False]
        ]
        matrix = BitMatrix.from_list(data)
        transposed = matrix.transposed()
        
        assert transposed.rows == 3
        assert transposed.cols == 2
        
        for i in range(2):
            for j in range(3):
                assert transposed[j, i] == matrix[i, j]
    
    def test_transpose_inplace(self):
        """Test in-place transpose."""
        data = [
            [True, False],
            [False, True],
            [True, True]
        ]
        matrix = BitMatrix.from_list(data)
        original_data = matrix.to_list()
        
        matrix.transpose_inplace()
        
        assert matrix.rows == 2
        assert matrix.cols == 3
        
        for i in range(3):
            for j in range(2):
                assert matrix[j, i] == original_data[i][j]
    
    def test_matrix_multiplication(self):
        """Test matrix multiplication."""
        # 2x3 @ 3x2 = 2x2
        a_data = [
            [True, False, True],
            [False, True, True]
        ]
        b_data = [
            [True, False],
            [True, True],
            [False, True]
        ]
        
        a = BitMatrix.from_list(a_data)
        b = BitMatrix.from_list(b_data)
        
        # Test @ operator
        c = a @ b
        assert c.rows == 2
        assert c.cols == 2
        
        # Test * operator
        c2 = a * b
        assert c == c2
        
        # Test matmul method
        c3 = a.matmul(b)
        assert c == c3
    
    def test_matrix_multiplication_dimension_mismatch(self):
        """Test matrix multiplication with incompatible dimensions."""
        a = BitMatrix(2, 3)
        b = BitMatrix(2, 2)  # Should be 3x2 for compatibility
        
        with pytest.raises(ValueError, match="Matrix multiplication failed"):
            _ = a @ b
    
    def test_inplace_matrix_multiplication(self):
        """Test in-place matrix multiplication."""
        a = BitMatrix.identity(3)
        b = BitMatrix.random(3, 3, seed=42)
        b_copy = b.copy()
        
        a @= b
        assert a == b_copy
    
    def test_matrix_power(self):
        """Test matrix exponentiation."""
        matrix = BitMatrix.identity(3)
        
        # Identity to any power should be identity
        powered = matrix ** 5
        assert powered == matrix
        
        # Any matrix to power 0 should be identity
        random_matrix = BitMatrix.random(3, 3, seed=42)
        power_zero = random_matrix ** 0
        identity = BitMatrix.identity(3)
        assert power_zero == identity
    
    def test_matrix_power_non_square(self):
        """Test matrix power with non-square matrix."""
        matrix = BitMatrix(2, 3)
        
        with pytest.raises(ValueError):
            _ = matrix ** 2
    
    def test_inverse(self):
        """Test matrix inversion."""
        # Create a known invertible matrix
        matrix = BitMatrix.identity(3)
        matrix[0, 1] = True
        matrix[1, 2] = True
        
        inverse = matrix.inverse()
        product = matrix @ inverse
        identity = BitMatrix.identity(3)
        
        assert product == identity
    
    def test_inverse_non_invertible(self):
        """Test inverse of non-invertible matrix."""
        # Zero matrix is not invertible
        matrix = BitMatrix.zeros(3, 3)
        
        with pytest.raises(ValueError):
            matrix.inverse()
    
    def test_gauss_elimination(self):
        """Test Gaussian elimination."""
        matrix = BitMatrix.random_invertible(4, seed=42)
        original = matrix.copy()
        
        matrix.gauss(full=False)
        
        # After partial Gaussian elimination, should be in row echelon form
        # The rank should be preserved
        assert matrix.rank() == original.rank()
        
        # Test full Gaussian elimination
        matrix2 = original.copy()
        matrix2.gauss(full=True)
        assert matrix2.rank() == original.rank()
    
    def test_rank(self):
        """Test rank computation."""
        # Identity matrix should have full rank
        identity = BitMatrix.identity(3)
        assert identity.rank() == 3
        
        # Zero matrix should have rank 0
        zero = BitMatrix.zeros(3, 3)
        assert zero.rank() == 0
    
    def test_stacking(self):
        """Test vertical and horizontal stacking."""
        a = BitMatrix.identity(2)
        b = BitMatrix.zeros(2, 2)
        
        # Vertical stack
        v_stacked = a.vstack(b)
        assert v_stacked.rows == 4
        assert v_stacked.cols == 2
        
        # Horizontal stack
        h_stacked = a.hstack(b)
        assert h_stacked.rows == 2
        assert h_stacked.cols == 4
    
    def test_stacking_dimension_mismatch(self):
        """Test stacking with incompatible dimensions."""
        a = BitMatrix(2, 3)
        b = BitMatrix(2, 2)
        
        with pytest.raises(ValueError):
            a.vstack(b)  # Different column counts
    
    def test_nullspace(self):
        """Test nullspace computation."""
        # For identity matrix, nullspace should be empty
        identity = BitMatrix.identity(3)
        nullspace = identity.nullspace()
        assert len(nullspace) == 0
        
        # Create a matrix with known nullspace
        matrix = BitMatrix.zeros(2, 3)
        matrix[0, 0] = True
        matrix[1, 1] = True
        # This matrix has rank 2, so nullspace dimension should be 1
        
        nullspace = matrix.nullspace()
        assert len(nullspace) > 0


class TestBitMatrixUtilities:
    """Test utility methods."""
    
    def test_copy(self):
        """Test matrix copying."""
        original = BitMatrix.random(3, 3, seed=42)
        copy = original.copy()
        
        assert copy == original
        assert copy is not original
        
        # Modifying copy shouldn't affect original
        copy[0, 0] = not copy[0, 0]
        assert copy != original
    
    def test_is_zero(self):
        """Test zero matrix checking."""
        zero = BitMatrix.zeros(3, 3)
        assert zero.is_zero()
        
        non_zero = BitMatrix.identity(3)
        assert not non_zero.is_zero()
    
    def test_to_list(self):
        """Test conversion to list of lists."""
        data = [
            [True, False, True],
            [False, True, False]
        ]
        matrix = BitMatrix.from_list(data)
        result = matrix.to_list()
        
        assert result == data
    
    def test_equality(self):
        """Test equality and inequality."""
        matrix1 = BitMatrix.identity(3)
        matrix2 = BitMatrix.identity(3)
        matrix3 = BitMatrix.zeros(3, 3)
        
        assert matrix1 == matrix2
        assert matrix1 != matrix3
        assert not (matrix1 == matrix3)
        assert matrix1.__ne__(matrix3)


class TestBitMatrixErrorHandling:
    """Test error handling and edge cases."""
    
    def test_build_with_invalid_function(self):
        """Test build with function that raises exceptions."""
        def bad_function(i, j):
            if i == 1 and j == 1:
                raise ValueError("Test error")
            return True
        
        # Should handle the exception and default to False
        matrix = BitMatrix.build(3, 3, bad_function)
        assert matrix[1, 1] == False  # Should be False due to exception
        assert matrix[0, 0] == True   # Should work normally
    
    def test_indexing_with_invalid_types(self):
        """Test indexing with invalid key types."""
        matrix = BitMatrix(3, 3)
        
        with pytest.raises(ValueError, match="Invalid index type"):
            _ = matrix["invalid"]
        
        with pytest.raises(ValueError, match="Invalid index type"):
            matrix["invalid"] = True


class TestBitMatrixIntegration:
    """Integration tests combining multiple operations."""
    
    def test_matrix_operations_chain(self):
        """Test chaining multiple operations."""
        # Create a random invertible matrix
        matrix = BitMatrix.random_invertible(4, seed=42)
        
        # Perform a series of operations
        inverse = matrix.inverse()
        product = matrix @ inverse
        identity = BitMatrix.identity(4)
        
        assert product == identity
        
        # Test with transpose
        transposed = matrix.transposed()
        double_transpose = transposed.transposed()
        assert double_transpose == matrix


if __name__ == "__main__":
    pytest.main([__file__])