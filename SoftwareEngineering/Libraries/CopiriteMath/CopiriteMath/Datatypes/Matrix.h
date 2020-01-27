#pragma once
#include "Transform.h"

#include <cassert>

#ifdef INCLUDE_DIRECTX_MATH
#define _XM_NO_INTRSECTS
#include <DirectXMath.h>
#endif


// A 2D arrangement of numbers in rows and columns.
// @template Columns - The amount of values along the X axis of this matrix.
// @template Rows - The amount of values along the Y axis of this matrix.
// @template Type - The datatype this matrix uses (default is float).
template <uint Columns, uint Rows, typename Type = float>
struct SMatrix
{
protected:
	/// Properties

	// The stored data in the matrix.
	STVector<Columns, Type> Data[Rows];


public:
	/// Constructors

	// Constructor, Default.
	SMatrix()
	{}

	// Constructor, Floods the entire matrix with a single floating point value.
	// @param Value - The float value used to flood the matrix.
	INLINE SMatrix(Type Value);

	// Constructor, Initiates the matrix using an array of vectors.
	// @note - The vector size must be the same size as the matrix's columns.
	// @note - The array size must be the same as the matrix's rows.
	// @param V - The array of vectors to assign to this matrix.
	INLINE SMatrix(STVector<Columns, Type> V[Rows]);



	/// Operators

	// Operator, Adds each component of this matrix by the same component on another matrix.
	INLINE SMatrix<Columns, Rows, Type> operator+(const SMatrix<Columns, Rows, Type>& Other) const;

	// Operator, Adds each component of this matrix by an inputted value.
	INLINE SMatrix<Columns, Rows, Type> operator+(const Type& Value) const;

	// Operator, Sets this matrix to the result of each component of this matrix added by the same component on another matrix.
	INLINE SMatrix<Columns, Rows, Type> operator+=(const SMatrix<Columns, Rows, Type>& Other);

	// Operator, Sets this matrix to the result of an inputted value being added to each component of this matrix.
	INLINE SMatrix<Columns, Rows, Type> operator+=(const Type& Value);

	// Operator, Subtracts each component of this matrix by the same component on another matrix.
	INLINE SMatrix<Columns, Rows, Type> operator-(const SMatrix<Columns, Rows, Type>& Other) const;

	// Operator, Subtracts each component of this matrix by an inputted value.
	INLINE SMatrix<Columns, Rows, Type> operator-(const Type& Value) const;

	// Operator, Sets this matrix to the result of each component of this matrix subtracted by the same component on another matrix.
	INLINE SMatrix<Columns, Rows, Type> operator-=(const SMatrix<Columns, Rows, Type>& Other);

	// Operator, Sets this matrix to the result of an inputted value being subtracted to each component of this matrix.
	INLINE SMatrix<Columns, Rows, Type> operator-=(const Type& Value);

	// Operator, Returns the result of this matrix being multiplied by another matrix.
	// Matrix multiplication requires that the row of this matrix is equal to the columns of the other matrix.
	// The return will be a matarix with the columns of the other matrix and the rows of this matrix.
	template <uint Columns2, uint Rows2, typename Type2>
	INLINE SMatrix<Columns2, Rows, Type> operator*(const SMatrix<Columns2, Rows2, Type2>& Other) const;

	// Operator, Multiplies each component of this matrix by an inputted value.
	INLINE SMatrix<Columns, Rows, Type> operator*(const Type& Value) const;

	// Operator, Multiplies a vector of the same size as the rows of this matrix, against this matrix.
	INLINE STVector<Rows, Type> operator*(const STVector<Rows, Type>& Value) const;

	// Operator, Sets this matrix to the result of this matrix being multiplied by another matrix.
	// Matrix multiplication requires that the row of this matrix is equal to the columns of the other matrix.
	// The return will be a matrix with the columns of the other matrix and the rows of this matrix.
	template <uint Columns2, uint Rows2, typename Type2>
	INLINE SMatrix<Columns2, Rows, Type> operator*=(const SMatrix<Columns2, Rows2, Type2>& Other);

	// Operator, Sets this matrix to the result of an inputted value being multiplied to each component of this matrix.
	INLINE SMatrix<Columns, Rows, Type> operator*=(const Type& Value);

	// Operator, Retrieves a row of the matrix based on an index.
	INLINE STVector<Columns, Type>& operator[](const uint& Index);

	// Operator, Retrieves a row of the matrix based on an index.
	INLINE STVector<Columns, Type> operator[](const uint& Index) const;

	// Operator, Retrieves a row of the matrix based on an axis.
	INLINE STVector<Columns, Type>& operator[](const EAxis& Axis);

	// Operator, Retrieves a row of the matrix based on an axis.
	INLINE STVector<Columns, Type> operator[](const EAxis& Axis) const;

	// Operator, Compares if this matrix has the same values as another matrix.
	// All values must be the same to return true.
	INLINE bool operator==(const SMatrix<Columns, Rows, Type>& Other) const;

	// Operator, Compares if this matrix has the same values as an inputted value.
	// All values must be the same to return true.
	INLINE bool operator==(const Type& Value) const;

	// Operator, Compares if this matrix does not have the same values as another matrix.
	// All values must be different to return true.
	INLINE bool operator!=(const SMatrix<Columns, Rows, Type>& Other) const;

	// Operator, Compares if this matrix does not have the same values as an inputted value.
	// All values must be different to return true.
	INLINE bool operator!=(const Type& Value) const;



	/// Conversions

	// Adds up all components of every column then puts it in a vector.
	INLINE STVector<Rows, Type> ToVector() const;



	/// Debug

	// Debug diagnostics handle for when a matrix contains NaN.
	INLINE void CheckNaN() const;

	// Checks if this matrix contains NaN.
	INLINE bool ContainsNaN() const;

	// Prints out all values in this matrix.
	INLINE void Print() const;



	/// Functions

	// Sets this matrix to the identity matrix.
	// @note - The last number in the matrix will always be 1.0f.
	// @param Value - The number that will be assigned along the diagonal.
	INLINE void Identity(Type Value = 1.0f);

	// Calculates the inverse of this matrix.
	// @note - This does not apply the inverse to this matrix.
	// @note - If the result is an empty matrix, then the inputted matrix could not be inverted.
	// @param Inverted - A pass-by-reference to return if this matrix was successfully inverted.
	// @return - The inverse matrix.
	INLINE SMatrix<Columns, Rows, Type> Inverse(bool& Inverted) const;

	// Calculates the inverse of this matrix.
	// @note - This does not apply the inverse to this matrix.
	// @note - If the result is an empty matrix, then the inputted matrix could not be inverted.
	// @return - The inverse matrix.
	INLINE SMatrix<Columns, Rows, Type> Inverse() const;

	// Calculates the transpose of this matrix.
	// @note - This does not apply the transpose to this matrix.
	// @return - The transposed matrix.
	INLINE SMatrix<Columns, Rows, Type> Transpose() const;

	// Calculates the determinant of this matrix.
	// @note - This function can only run if the matrix it is being called on is a squared matrix.
	// @param Size - The X and Y size of the matrix.
	// @return - The calculated determinant of this matrix.
	INLINE float Determinant(const uint& Size = Columns) const;

	// Creates a matrix by taking the transpose of the cofactor matrix.
	INLINE SMatrix<Columns, Rows, Type> Adjoint() const;



	/// Setters

	// Sets a vertical row of vectors in a specified column.
	// @param Column - The column to be assigned.
	// @param Vec - The vector to assign at the given column.
	INLINE void SetColumn(const uint& Column, const STVector<Rows, Type>& Vec);

	// Sets a vertical row of vectors in a specified column.
	// @param Axis - The column to be assigned.
	// @param Vec - The vector to assign at the given column.
	INLINE void SetColumn(const EAxis& Axis, const STVector<Rows, Type>& Vec);



	/// Getters

	// Returns the number of rows this matrix has.
	INLINE uint GetRowCount() const { return Rows; }

	// Returns the number of columns this matrix has.
	INLINE uint GetColumnCount() const { return Columns; }

	// Returns a value at a specific element of the matrix.
	// @param X - The column to be retrieved.
	// @param Y - The row to be retrieved.
	// @return - The value at the specified row and column.
	INLINE Type& GetValue(const uint& X, const uint& Y);

	// Returns a value at a specific element of the matrix.
	// @param X - The column to be retrieved.
	// @param Y - The row to be retrieved.
	// @return - The value at the specified row and column.
	INLINE Type GetValue(const uint& X, const uint& Y) const;

	// Returns a value at a specific element of the matrix.
	// @param X - The column to be retrieved.
	// @param Y - The row to be retrieved.
	// @return - The value at the specified row and column.
	INLINE Type& GetValue(const EAxis& Column, const EAxis& Row);

	// Returns a value at a specific element of the matrix.
	// @param X - The column to be retrieved.
	// @param Y - The row to be retrieved.
	// @return - The value at the specified row and column.
	INLINE Type GetValue(const EAxis& Column, const EAxis& Row) const;

	// Gets a matrix by removing a row and column from this matrix.
	// @param Row - The row to be removed.
	// @param Column - The column to be removed.
	// @param Size - The column and row size of the new matrix.
	// @return - Returns a copy of this matrix after being confactored.
	INLINE SMatrix<Columns, Rows, Type> GetCofactor(const uint& Row, const uint& Column, const uint& Size) const;

	// Returns a row of the matrix.
	// @param Row - The row to be retrieved.
	// @return - The row at the inputted index.
	INLINE STVector<Rows, Type>& GetRow(const uint& Row);

	// Returns a row of the matrix.
	// @param Row - The row to be retrieved.
	// @return - The row at the inputted index.
	INLINE STVector<Rows, Type> GetRow(const uint& Row) const;

	// Returns a row of the matrix.
	// @param Row - The row to be retrieved.
	// @return - The row at the inputted index.
	INLINE STVector<Rows, Type>& GetRow(const EAxis& Row);

	// Returns a row of the matrix.
	// @param Row - The row to be retrieved.
	// @return - The row at the inputted index.
	INLINE STVector<Rows, Type> GetRow(const EAxis& Row) const;

	// Returns a column of the matrix.
	// @param Column - The column to be retrieved.
	// @return - The column at the inputted index.
	INLINE STVector<Columns, Type>& GetColumn(const uint& Column);

	// Returns a column of the matrix.
	// @param Column - The column to be retrieved.
	// @return - The column at the inputted index.
	INLINE STVector<Columns, Type> GetColumn(const uint& Column) const;

	// Returns a column of the matrix.
	// @param Column - The column to be retrieved.
	// @return - The column at the inputted index.
	INLINE STVector<Columns, Type>& GetColumn(const EAxis& Column);

	// Returns a column of the matrix.
	// @param Column - The column to be retrieved.
	// @return - The column at the inputted index.
	INLINE STVector<Columns, Type> GetColumn(const EAxis& Column) const;

	// Returns a vector at an inputted axis of the matrix.
	// @param Axis - The axis to retrieve.
	// @return - The vector at the given axis.
	INLINE STVector<Columns, Type> GetScaledAxis(EAxis Axis) const;



	/// Statics

	// Converts an inputted vector to a matrix.
	// @param Vec - The vector to convert.
	// @return - The resulting matrix from the vector.
	static INLINE SMatrix<1, Rows, Type> ToMatrix(STVector<Rows, float> Vec)
	{
		SMatrix<1, Rows, Type> Result;
		for (uint i = 0; i < Rows; ++i)
		{
			Result[i][0] = Vec[i];
		}
		return Result;
	}

	// Returns an identity matrix.
	static INLINE SMatrix<Columns, Rows, Type> MatrixIdentity()
	{
		SMatrix<Columns, Rows, Type> Matrix;
		Matrix.Identity();
		return Matrix;
	}
};


// A floating point matrix with 1 column and 2 rows.
typedef SMatrix<1, 2, float> SMatrix12;

// A floating point matrix with 2 columns and 1 row.
typedef SMatrix<2, 1, float> SMatrix21;

// A floating point matrix with 1 column and 3 rows.
typedef SMatrix<1, 3, float> SMatrix13;

// A floating point matrix with 3 columns and row.
typedef SMatrix<3, 1, float> SMatrix31;

// A floating point matrix with 1 column and 4 rows.
typedef SMatrix<1, 4, float> SMatrix14;

// A floating point matrix with 4 columns and 1 row.
typedef SMatrix<4, 1, float> SMatrix41;

// A floating point matrix with 2 columns and 2 rows.
typedef SMatrix<2, 2, float> SMatrix2;

// A floating point matrix with 2 columns and 3 rows.
typedef SMatrix<2, 3, float> SMatrix23;

// A floating point matrix with 3 columns and 2 rows.
typedef SMatrix<3, 2, float> SMatrix32;

// A floating point matrix with 2 columns and 4 rows.
typedef SMatrix<2, 4, float> SMatrix24;

// A floating point matrix with 4 columns and 2 rows.
typedef SMatrix<4, 2, float> SMatrix42;

// A floating point matrix with 3 columns and 4 rows.
typedef SMatrix<3, 4, float> SMatrix34;

// A floating point matrix with 4 columns and 3 rows.
typedef SMatrix<4, 3, float> SMatrix43;

// A double type matrix with 1 column and 2 rows.
typedef SMatrix<1, 2, double> SMatrix12d;

// A double type matrix with 2 columns and 1 row.
typedef SMatrix<2, 1, double> SMatrix21d;

// A double type matrix with 1 column and 3 rows.
typedef SMatrix<1, 3, double> SMatrix13d;

// A double type matrix with 3 columns and row.
typedef SMatrix<3, 1, double> SMatrix31d;

// A double type matrix with 1 column and 4 rows.
typedef SMatrix<1, 4, double> SMatrix14d;

// A double type matrix with 4 columns and 1 row.
typedef SMatrix<4, 1, double> SMatrix41d;

// A double type matrix with 2 columns and 2 rows.
typedef SMatrix<2, 2, double> SMatrix2d;

// A double type matrix with 2 columns and 3 rows.
typedef SMatrix<2, 3, double> SMatrix23d;

// A double type matrix with 3 columns and 2 rows.
typedef SMatrix<3, 2, double> SMatrix32d;

// A double type matrix with 2 columns and 4 rows.
typedef SMatrix<2, 4, double> SMatrix24d;

// A double type matrix with 4 columns and 2 rows.
typedef SMatrix<4, 2, double> SMatrix42d;

// A double type matrix with 3 columns and 4 rows.
typedef SMatrix<3, 4, double> SMatrix34d;

// A double type matrix with 4 columns and 3 rows.
typedef SMatrix<4, 3, double> SMatrix43d;

// An integer matrix with 1 column and 2 rows.
typedef SMatrix<1, 2, int> SMatrix12i;

// An integer matrix with 2 columns and 1 row.
typedef SMatrix<2, 1, int> SMatrix21i;

// An integer matrix with 1 column and 3 rows.
typedef SMatrix<1, 3, int> SMatrix13i;

// An integer matrix with 3 columns and row.
typedef SMatrix<3, 1, int> SMatrix31i;

// An integer matrix with 1 column and 4 rows.
typedef SMatrix<1, 4, int> SMatrix14i;

// An integer matrix with 4 columns and 1 row.
typedef SMatrix<4, 1, int> SMatrix41i;

// An integer matrix with 2 columns and 2 rows.
typedef SMatrix<2, 2, int> SMatrix2i;

// An integer matrix with 2 columns and 3 rows.
typedef SMatrix<2, 3, int> SMatrix23i;

// An integer matrix with 3 columns and 2 rows.
typedef SMatrix<3, 2, int> SMatrix32i;

// An integer matrix with 2 columns and 4 rows.
typedef SMatrix<2, 4, int> SMatrix24i;

// An integer matrix with 4 columns and 2 rows.
typedef SMatrix<4, 2, int> SMatrix42i;

// An integer matrix with 3 columns and 4 rows.
typedef SMatrix<3, 4, int> SMatrix34i;

// An integer matrix with 4 columns and 3 rows.
typedef SMatrix<4, 3, int> SMatrix43i;


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type>::SMatrix(Type Value)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i] = Value;
	}
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type>::SMatrix(STVector<Columns, Type> V[Rows])
	:Data{ V }
{}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator+(const SMatrix<Columns, Rows, Type>& M) const
{
	SMatrix<Columns, Rows, Type> Result;
	for (uint i = 0; i < Rows; ++i)
	{
		Result[i] = Data[i] + M[i];
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator+(const Type& Value) const
{
	SMatrix<Columns, Rows, Type> Result;
	for (uint i = 0; i < Rows; ++i)
	{
		Result[i] = Data[i] + Value;
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator+=(const SMatrix<Columns, Rows, Type>& M)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i] += M[i];
	}
	return *this;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator+=(const Type& Value)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i] += Value;
	}
	return *this;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator-(const SMatrix<Columns, Rows, Type>& M) const
{
	SMatrix<Columns, Rows, Type> Result;
	for (uint i = 0; i < Rows; ++i)
	{
		Result[i] = Data[i] - M[i];
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator-(const Type& Value) const
{
	SMatrix<Columns, Rows, Type> Result;
	for (uint i = 0; i < Rows; ++i)
	{
		Result[i] = Data[i] - Value;
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator-=(const SMatrix<Columns, Rows, Type>& M)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i] -= M[i];
	}
	return *this;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator-=(const Type& Value)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i] -= Value;
	}
	return *this;
}


template <uint Columns, uint Rows, typename Type>
template <uint Columns2, uint Rows2, typename Type2>
INLINE SMatrix<Columns2, Rows, Type> SMatrix<Columns, Rows, Type>::operator*(const SMatrix<Columns2, Rows2, Type2>& M) const
{
	ASSERT(Columns == Rows2, "These matrices can not be multiplied together.");
	SMatrix<Columns2, Rows> Result{ 0.0f };
	for (uint y = 0; y < Rows; ++y)
	{
		for (uint x = 0; x < Columns2; ++x)
		{
			for (uint i = 0; i < Columns; ++i)
			{
				Result[y][x] += (Data[y][i] * M[i][x]);
			}
		}
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator*(const Type& Value) const
{
	SMatrix<Columns, Rows, Type> Result;
	for (uint i = 0; i < Rows; ++i)
	{
		Result[i] = Data[i] * Value;
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Rows, Type> SMatrix<Columns, Rows, Type>::operator*(const STVector<Rows, Type>& V) const
{
	return (*this * SMatrix<1, Rows>::ToMatrix(V)).ToVector();
}


template <uint Columns, uint Rows, typename Type>
template <uint Columns2, uint Rows2, typename Type2>
INLINE SMatrix<Columns2, Rows, Type> SMatrix<Columns, Rows, Type>::operator*=(const SMatrix<Columns2, Rows2, Type2>& M)
{
	*this = *this * M;
	return *this;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::operator*=(const Type& Value)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i] *= Value;
	}
	return *this;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type>& SMatrix<Columns, Rows, Type>::operator[](const uint& Index)
{
	return Data[Index];
}

template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type> SMatrix<Columns, Rows, Type>::operator[](const uint& Index) const
{
	return Data[Index];
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type>& SMatrix<Columns, Rows, Type>::operator[](const EAxis& Index)
{
	return Data[Index];
}

template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type> SMatrix<Columns, Rows, Type>::operator[](const EAxis& Index) const
{
	return Data[Index];
}


template <uint Columns, uint Rows, typename Type>
INLINE bool SMatrix<Columns, Rows, Type>::operator==(const SMatrix<Columns, Rows, Type>& Other) const
{
	for (uint i = 0; i < Rows; ++i)
	{
		if (Data[i] != Other[i]) return false;
	}
	return true;
}


template <uint Columns, uint Rows, typename Type>
INLINE bool SMatrix<Columns, Rows, Type>::operator==(const Type& Value) const
{
	for (uint i = 0; i < Rows; ++i)
	{
		if (Data[i] != Value) return false;
	}
	return true;
}


template <uint Columns, uint Rows, typename Type>
INLINE bool SMatrix<Columns, Rows, Type>::operator!=(const SMatrix<Columns, Rows, Type>& Other) const
{
	for (uint i = 0; i < Rows; ++i)
	{
		if (Data[i] == Other[i]) return false;
	}
	return true;
}


template <uint Columns, uint Rows, typename Type>
INLINE bool SMatrix<Columns, Rows, Type>::operator!=(const Type& Value) const
{
	for (uint i = 0; i < Rows; ++i)
	{
		if (Data[i] == Value) return false;
	}
	return true;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Rows, Type> SMatrix<Columns, Rows, Type>::ToVector() const
{
	STVector<Rows, Type> Result{ 0.0f };
	for (uint y = 0; y < Rows; ++y)
	{
		for (uint x = 0; x < Columns; ++x)
		{
			Result[y] += Data[y][x];
		}
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE void SMatrix<Columns, Rows, Type>::CheckNaN() const
{
	if (ContainsNaN())
	{
		printf("Matrix contains NaN\n");
		*const_cast<SMatrix<Columns, Rows, Type>*>(this) = SMatrix<Columns, Rows, Type>::Identity();
	}
}


template <uint Columns, uint Rows, typename Type>
INLINE bool SMatrix<Columns, Rows, Type>::ContainsNaN() const
{
	for (uint y = 0; y < Rows; ++y)
	{
		for (uint x = 0; x < Columns; ++x)
		{
			if (TMath::IsFinite(Data[y][x])) return false;
		}
	}
	return true;
}


template <uint Columns, uint Rows, typename Type>
INLINE void SMatrix<Columns, Rows, Type>::Print() const
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i].Print();
	}
	printf("\n");
}


template <uint Columns, uint Rows, typename Type>
INLINE void SMatrix<Columns, Rows, Type>::Identity(Type Value)
{
	ASSERT(Columns == Rows, "The matrix must have equal rows and columns to be set to an indentity matrix.");
	for (uint y = 0; y < Rows; ++y)
	{
		for (uint x = 0; x < Columns; ++x)
		{
			Data[y][x] = ((x == y) * Value);
		}
	}
	Data[Columns - 1][Rows - 1] = (Type)1.0f;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::Inverse(bool& Inverted) const
{
	ASSERT(Columns == Rows, "This matrix size must be the same in both dimensions to calculate the inverse.");
	float Det{ Determinant() };
	if (Det = 0.0f)
	{
		Inverted = false;
		return 0.0f;
	}

	SMatrix<Columns, Rows, Type> Result;
	SMatrix<Columns, Rows, Type> Adj{ Adjoint() };
	for (uint x = 0; x < Columns; ++x)
	{
		for (uint y = 0; y < Rows; ++y)
		{
			Result[x][y] = Adj[x][y] / Det;
		}
	}
	Inverted = true;
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::Inverse() const
{
	ASSERT(Columns == Rows, "This matrix size must be the same in both dimensions to calculate the inverse.");
	float Det{ Determinant() };
	if (Det = 0.0f)
	{
		return 0.0f;
	}

	SMatrix<Columns, Rows, Type> Result;
	SMatrix<Columns, Rows, Type> Adj{ Adjoint() };
	for (uint x = 0; x < Columns; ++x)
	{
		for (uint y = 0; y < Rows; ++y)
		{
			Result[x][y] = Adj[x][y] / Det;
		}
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::Transpose() const
{
	// TODO
}


template <uint Columns, uint Rows, typename Type>
INLINE float SMatrix<Columns, Rows, Type>::Determinant(const uint& Size) const
{
	ASSERT(Columns == Rows, "This matrix must be the same in both dimensions to calculate the determinate.");
	float Result{ 0.0f };

	if (Size == 1) return Data[0][0];

	SMatrix<Columns, Rows, Type> Matrix;
	int8 Sign{ 1 };

	for (uint i = 0; i < Size; ++i)
	{
		Matrix = GetCofactor(0, i, Size);
		Result += Sign * Data[0][i] * Matrix.Determinant(Size - 1);
		Sign = -Sign;
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::Adjoint() const
{
	ASSERT(Columns == Rows, "This matrix size must be the same in both dimensions to calculate the adjoint.");
	SMatrix<Columns, Rows> Result;

	if (Columns == 1)
	{
		Result[0][0] = 1.0f;
		return Result;
	}

	int8 Sign{ 1 };
	SMatrix<Columns, Rows> Matrix;
	for (uint x = 0; x < Columns; ++x)
	{
		for (uint y = 0; y < Rows; ++y)
		{
			Matrix = GetCofactor(x, y, Columns);
			Sign = ((x + y) % 2 == 0) ? 1 : -1;
			Result[y][x] = Sign * Matrix.Determinant(Columns - 1);
		}
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE void SMatrix<Columns, Rows, Type>::SetColumn(const uint& Column, const STVector<Rows, Type>& Vec)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i][Column] = Vec[i];
	}
}


template <uint Columns, uint Rows, typename Type>
INLINE void SMatrix<Columns, Rows, Type>::SetColumn(const EAxis& Axis, const STVector<Rows, Type>& Vec)
{
	for (uint i = 0; i < Rows; ++i)
	{
		Data[i][Axis] = Vec[i];
	}
}


template <uint Columns, uint Rows, typename Type>
INLINE Type& SMatrix<Columns, Rows, Type>::GetValue(const uint& X, const uint& Y)
{
	return Data[Y][X];
}


template <uint Columns, uint Rows, typename Type>
INLINE Type SMatrix<Columns, Rows, Type>::GetValue(const uint& X, const uint& Y) const
{
	return Data[Y][X];
}


template <uint Columns, uint Rows, typename Type>
INLINE Type& SMatrix<Columns, Rows, Type>::GetValue(const EAxis& Column, const EAxis& Row)
{
	return Data[Row][Column];
}


template <uint Columns, uint Rows, typename Type>
INLINE Type SMatrix<Columns, Rows, Type>::GetValue(const EAxis& Column, const EAxis& Row) const
{
	return Data[Row][Column];
}


template <uint Columns, uint Rows, typename Type>
INLINE SMatrix<Columns, Rows, Type> SMatrix<Columns, Rows, Type>::GetCofactor(const uint& Row, const uint& Column, const uint& Size) const
{
	int8 i{ 0 };
	int8 j{ 0 };

	SMatrix<Columns, Rows> Result;

	for (uint y = 0; y < Size; ++y)
	{
		for (uint x = 0; x < Size; ++x)
		{
			if (y != Row && x != Column)
			{
				Result[i][j++] = Data[y][x];
				if (j == Size - 1)
				{
					j = 0;
					++i;
				}
			}
		}
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Rows, Type>& SMatrix<Columns, Rows, Type>::GetRow(const uint& Row)
{
	return Data[Row];
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Rows, Type> SMatrix<Columns, Rows, Type>::GetRow(const uint& Row) const
{
	return Data[Row];
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Rows, Type>& SMatrix<Columns, Rows, Type>::GetRow(const EAxis& Row)
{
	return Data[Row];
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Rows, Type> SMatrix<Columns, Rows, Type>::GetRow(const EAxis& Row) const
{
	return Data[Row];
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type>& SMatrix<Columns, Rows, Type>::GetColumn(const uint& Column)
{
	STVector<Columns, Type> Result;
	for (uint i = 0; i < Columns; ++i)
	{
		Result[i] = Data[i][Column];
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type> SMatrix<Columns, Rows, Type>::GetColumn(const uint& Column) const
{
	STVector<Columns, Type> Result;
	for (uint i = 0; i < Columns; ++i)
	{
		Result[i] = Data[i][Column];
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type>& SMatrix<Columns, Rows, Type>::GetColumn(const EAxis& Column)
{
	STVector<Columns, Type> Result;
	for (uint i = 0; i < Columns; ++i)
	{
		Result[i] = Data[i][Column];
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type> SMatrix<Columns, Rows, Type>::GetColumn(const EAxis& Column) const
{
	STVector<Columns, Type> Result;
	for (uint i = 0; i < Columns; ++i)
	{
		Result[i] = Data[i][Column];
	}
	return Result;
}


template <uint Columns, uint Rows, typename Type>
INLINE STVector<Columns, Type> SMatrix<Columns, Rows, Type>::GetScaledAxis(EAxis Axis) const
{
	if (Axis == EAxis::W) return STVector<Columns, Type>{ 0.0f };
	return Data[Axis];
}



// A 2D arrangement of numbers in rows and columns.
// This matrix type contains 3 rows and 3 columns.
// This matrix is mostly used for 2D math and has a lot more functions independent from the other matrix types.
struct SMatrix3 :public SMatrix<3, 3, float>
{
public:
	/// Constructors

	// Constructor, Default.
	SMatrix3()
	{}

	// Constructor, Copies the values of a 3x3 matrix into this matrix.
	// @param Other - The original matrix to copy the values from.
	INLINE SMatrix3(SMatrix<3, 3, float> Other);

	// Constructor, Initiates the matrix with a location, rotation and scale.
	// @note - This constructor multiplies the matrices: Scale * Rotation * Translation.
	// @param Transform - The location, rotation and scale this matrix should be created at.
	INLINE SMatrix3(STransform Transform);

	// Constructor, Initiates the matrix with a location, rotation and scale.
	// @note - This constructor multiplies the matrices: Scale * Rotation * Translation.
	// @param Location - The location this matrix should be placed at.
	// @param Rotation - The rotation this matrix should be at.
	// @param Scale - The scale this matrix should be at.
	INLINE SMatrix3(SVector Location, SQuaternion Rotation, SVector Scale);



	/// Conversions

	// Converts this matrix to a transformation.
	INLINE STransform ToTransform() const;



	/// Functions

	// Moves this matrix by an inputted amount.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @return - The resulting matrix from the translation.
	INLINE SMatrix3 Translate(const SVector& Translation) const;

	// Moves this matrix by an inputted amount.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @param W - The W component.
	// @return - The resulting matrix from the translation.
	INLINE SMatrix3 Translate(const SVector2& Translation, const float& W = 1.0f) const;

	// Moves this matrix by an inputted amount.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param W - The W component.
	// @return - The resulting matrix from the translationg.
	INLINE SMatrix3 Translate(const float& X, const float& Y, const float& W = 1.0f) const;

	// Moves this matrix by an inputted amount and sets this matrix to the returning result.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @return - The resulting matrix from the translation.
	INLINE SMatrix3 SetTranslate(const SVector& Translation);

	// Moves this matrix by an inputted amount and sets this matrix to the returning result.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @param W - The W component.
	// @return - The resulting matrix from the translation.
	INLINE SMatrix3 SetTranslate(const SVector2& Translation, const float& W = 1.0f);

	// Moves this matrix by an inputted amount and sets this matrix to the returning result.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param W - The W component.
	// @return - The resulting matrix from the translationg.
	INLINE SMatrix3 SetTranslate(const float& X, const float& Y, const float& W);

	// Makes this matrix a translation matrix and represents the inputted location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @return - The translation matrix.
	INLINE SMatrix3 ToTranslation(const SVector& Translation);

	// Makes this matrix a translation matrix and represents the inputted location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @param W - The W component.
	// @return - The translation matrix.
	INLINE SMatrix3 ToTranslation(const SVector2& Translation, const float& W = 1.0f);

	// Makes this matrix a translation matrix and represents the inputted location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param W - The W component.
	// @return - The translation matrix.
	INLINE SMatrix3 ToTranslation(const float& X, const float& Y, const float& W = 1.0f);

	// Rotates along the Z axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 RotateRoll(const float& Angle) const;

	// Rotates along the Z axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 SetRotateRoll(const float& Angle);

	// Makes this matrix a rotation matrix along the Z axis and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 ToRotationRoll(const float& Angle);

	// Rotates along the Z axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 RotateZ(const float& Angle) const;

	// Rotates along the Z axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 SetRotateZ(const float& Angle);

	// Makes this matrix a rotation matrix along the Z axis and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 ToRotationZ(const float& Angle);

	// Rotates this matrix along all axis and returns the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The amount to rotate in the Z axis.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 Rotate(const float& Rotation) const;

	// Rotates this matrix along all axis and returns the result and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The amount to rotate in the Z axis.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 SetRotate(const float& Rotation);

	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The amount to rotate in the Z axis.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix3 ToRotation(const float& Rotation);

	// Scales this matrix and returns the resulting matrix.
	// @param InScale - The scale this matrix will be multiplied by.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix3 Scale(const SVector2& InScale) const;

	// Scales this matrix and returns the resulting matrix.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix3 Scale(const float& X, const float& Y) const;

	// Scales this matrix and returns the resulting matrix and sets this matrix the resulting matrix.
	// @param InScale - The scale this matrix will be multiplied by.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix3 SetScale(const SVector2& InScale);

	// Scales this matrix and returns the resulting matrix and sets this matrix the resulting matrix.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix3 SetScale(const float& X, const float& Y);

	// Makes this matrix a scale matrix and represents the inputted scale.
	// @param InScale - The scale this matrix will be multiplied by.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix3 ToScale(const SVector2& NewScale);

	// Makes this matrix a scale matrix and represents the inputted scale.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix3 ToScale(const float& X, const float& Y);

	// Transforms this matrix using a location, rotation and scale.
	// @param InTransform - The location, rotation and scale to multiply against this matrix.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix3 Transform(const STransform& InTransform) const;

	// Transforms this matrix using a location, rotation and scale.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The scale to be scaled in.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix3 Transform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale) const;

	// Transforms this matrix using a location, rotationa and scale and set this matrix to the resulting matrix.
	// @param InTransform - The location, rotation and scale to multiply against this matrix.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix3 SetTransform(const STransform& InTransform);

	// Transforms this matrix using a location, rotationa and scale and set this matrix to the resulting matrix.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The scale to be scaled in.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix3 SetTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale);

	// Makes this matrix a transform matrix and represents the inputted transformation.
	// @param InTransform - The location, rotation and scale to multiply against this matrix.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix3 ToTransform(const STransform& NewTransform);

	// Makes this matrix a transform matrix and represents the inputted transformation.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The scale to be scaled in.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix3 ToTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale);



	/// Getters

	// Returns this matrix as a transformation.
	INLINE STransform GetTransform() const;

	// Returns this matrix's location.
	INLINE SVector2 GetLocation() const;

	// Returns this matrix's rotation.
	INLINE float GetRotation() const;

	// Returns this matrix's scale.
	INLINE SVector2 GetScale() const;



	/// Statics

	// Returns a matrix translated to a location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The translation to be moved by.
	// @return - The translation matrix.
	static INLINE SMatrix3 MatrixTranslate(const SVector& Translation)
	{
		SMatrix3 Result;
		return Result.ToTranslation(Translation);
	}


	// Returns a matrix translated to a location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param W - The W value. 
	// @return - The translation matrix.
	static INLINE SMatrix3 MatrixTranslate(const float& X, const float& Y, const float& W = 1.0f)
	{
		SMatrix3 Result;
		return Result.ToTranslation(X, Y, W);
	}


	// Returns a matrix rotated in the Z axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix3 MatrixRotateRoll(const float& Angle)
	{
		SMatrix3 Result;
		return Result.ToRotationRoll(Angle);
	}


	// Returns a matrix rotated in the Z axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix3 MatrixRotateZ(const float& Angle)
	{
		SMatrix3 Result;
		return Result.ToRotationRoll(Angle);
	}


	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix3 MatrixRotate(const float& Rotation)
	{
		SMatrix3 Result;
		return Result.ToRotation(Rotation);
	}


	// Returns a matrix scaled by a value in all dimensions.
	// @param InScale - The amount to be scaled.
	// @return - The resulting matrix from the scaling.
	static INLINE SMatrix3 MatrixScale(const SVector2& Scale)
	{
		SMatrix3 Result;
		return Result.ToScale(Scale);
	}


	// Returns a matrix scaled by a value in all dimensions.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @return - The resulting matrix from the scaling.
	static INLINE SMatrix3 MatrixScale(const float& X, const float& Y)
	{
		SMatrix3 Result;
		return Result.ToScale(X, Y);
	}


	// Returns a matrix scaled, rotated and translated by an inputted transformation.
	// @param Transform - The location, rotation and scale this matrix will be transformed.
	// @return - The resulting matrix after the transformation.
	static INLINE SMatrix3 MatrixTransform(const STransform& Transform)
	{
		SMatrix3 Result;
		return Result.ToTransform(Transform);
	}


	// Returns a matrix scaled, rotated and translated by an inputted transformation.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The scale to be scaled in.
	// @return - The resulting matrix after the transformation.
	static INLINE SMatrix3 MatrixTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& Scale)
	{
		SMatrix3 Result;
		return Result.ToTransform(Location, Rotation, Scale);
	}
};


INLINE SMatrix3::SMatrix3(SMatrix<3, 3, float> Other)
{
	*this = Other;
}


INLINE SMatrix3::SMatrix3(STransform Transform)
{
	ToTransform(Transform);
}


INLINE SMatrix3::SMatrix3(SVector Location, SQuaternion Rotation, SVector Scale)
{
	ToTransform(Location, Rotation, Scale);
}


INLINE SMatrix3 SMatrix3::Translate(const SVector& Translation) const
{
	SMatrix3 Result;
	Result.Identity();
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Result[i][2] = Translation[i];
	}
	return *this * Result;
}


INLINE SMatrix3 SMatrix3::Translate(const SVector2& Translation, const float& W) const
{
	return Translate(SVector{ Translation, W });
}


INLINE SMatrix3 SMatrix3::Translate(const float& X, const float& Y, const float& W) const
{
	return Translate(SVector{ X, Y, W });
}


INLINE SMatrix3 SMatrix3::SetTranslate(const SVector& Translation)
{
	*this = Translate(Translation);
	return *this;
}


INLINE SMatrix3 SMatrix3::SetTranslate(const SVector2& Translation, const float& W)
{
	*this = Translate(Translation, W);
	return *this;
}


INLINE SMatrix3 SMatrix3::SetTranslate(const float& X, const float& Y, const float& W)
{
	*this = Translate(X, Y, W);
	return *this;
}


INLINE SMatrix3 SMatrix3::ToTranslation(const SVector& Translation)
{
	Identity();
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Data[i][2] = Translation[i];
	}
	return *this;
}


INLINE SMatrix3 SMatrix3::ToTranslation(const SVector2& Translation, const float& W)
{
	return ToTranslation(SVector{ Translation, W });
}


INLINE SMatrix3 SMatrix3::ToTranslation(const float& X, const float& Y, const float& W)
{
	return ToTranslation(SVector{ X, Y, W });
}


INLINE SMatrix3 SMatrix3::RotateRoll(const float& Angle) const
{
	SMatrix3 Result;
	Result.Identity();

	Result[0][0] = TMath::Cos(Angle);
	Result[0][1] = TMath::Sin(Angle);
	Result[1][0] = -TMath::Sin(Angle);
	Result[1][1] = TMath::Cos(Angle);

	return Result;
}


INLINE SMatrix3 SMatrix3::SetRotateRoll(const float& Angle)
{
	*this = RotateRoll(Angle);
	return *this;
}


INLINE SMatrix3 SMatrix3::ToRotationRoll(const float& Angle)
{
	Identity();
	Data[0][0] = TMath::Cos(Angle);
	Data[0][1] = TMath::Sin(Angle);
	Data[1][0] = -TMath::Sin(Angle);
	Data[1][1] = TMath::Cos(Angle);

	return *this;
}


INLINE SMatrix3 SMatrix3::RotateZ(const float& Angle) const
{
	return RotateRoll(Angle);
}


INLINE SMatrix3 SMatrix3::SetRotateZ(const float& Angle)
{
	return SetRotateRoll(Angle);
}


INLINE SMatrix3 SMatrix3::ToRotationZ(const float& Angle)
{
	return ToRotationRoll(Angle);
}


INLINE SMatrix3 SMatrix3::Rotate(const float& Rotation) const
{
	SMatrix3 Result;
	Result = RotateRoll(Rotation);
	return *this * Result;
}


INLINE SMatrix3 SMatrix3::SetRotate(const float& Rotation)
{
	*this = Rotate(Rotation);
	return *this;
}


INLINE SMatrix3 SMatrix3::ToRotation(const float& Rotation)
{
	*this = RotateRoll(Rotation);
	return *this;
}


INLINE SMatrix3 SMatrix3::Scale(const SVector2& InScale) const
{
	SMatrix3 Result{ 0.0f };
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Result[i][i] = InScale[i];
	}
	Result[3][3] = 1.0f;
	return *this * Result;
}


INLINE SMatrix3 SMatrix3::Scale(const float& X, const float& Y) const
{
	return Scale(SVector{ X, Y, 1.0f });
}


INLINE SMatrix3 SMatrix3::SetScale(const SVector2& InScale)
{
	*this = Scale(InScale);
	return *this;
}


INLINE SMatrix3 SMatrix3::SetScale(const float& X, const float& Y)
{
	*this = Scale(X, Y);
	return *this;
}


INLINE SMatrix3 SMatrix3::ToScale(const SVector2& NewScale)
{
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Data[i][i] = NewScale[i];
	}
	Data[2][2] = 1.0f;
	return *this;
}


INLINE SMatrix3 SMatrix3::ToScale(const float& X, const float& Y)
{
	return ToScale(SVector2{ X, Y });
}


INLINE SMatrix3 SMatrix3::Transform(const STransform& InTransform) const
{
	SMatrix3 Result{ Scale(InTransform.Scale) * Rotate(InTransform.Rotation.Z) * Translate(InTransform.Location) };
	return *this * Result;
}


INLINE SMatrix3 SMatrix3::Transform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale) const
{
	return Transform(STransform{ Location, Rotation, InScale });
}


INLINE SMatrix3 SMatrix3::SetTransform(const STransform& InTransform)
{
	*this = Transform(InTransform);
	return *this;
}


INLINE SMatrix3 SMatrix3::SetTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale)
{
	*this = Transform(Location, Rotation, InScale);
	return *this;
}


INLINE SMatrix3 SMatrix3::ToTransform(const STransform& NewTransform)
{
	SMatrix3 TranslateMat;
	SMatrix3 RotationMat;
	SMatrix3 ScaleMat;

	TranslateMat.ToTranslation(NewTransform.Location);
	RotationMat.ToRotation(NewTransform.Rotation.Z);
	ScaleMat.ToScale(NewTransform.Scale);

	*this = ScaleMat * RotationMat * TranslateMat;
	return *this;
}


INLINE SMatrix3 SMatrix3::ToTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale)
{
	return ToTransform(STransform{ Location, Rotation, InScale });
}


INLINE STransform SMatrix3::GetTransform() const
{
	return STransform{ GetLocation(), GetRotation(), GetScale() };
}


INLINE SVector2 SMatrix3::GetLocation() const
{
	return GetColumn(2);
}


INLINE float SMatrix3::GetRotation() const
{
	float Result{ 0.0f };
	Result = TMath::ACos(Data[0][0]);
	Result *= TMath::ASin(Data[0][1]);
	Result *= -TMath::ASin(Data[1][0]);
	Result *= TMath::ACos(Data[1][1]);
	return Result;

	//float R22{ Data[2][2] };
	//float Sum10 = Data[1][1] + Data[0][0];
	//float Opr22 = 1.0f + R22;
	//if (Sum10 <= 0.0f)
	//{
	//	float FourZSqr = Opr22 - Sum10;
	//	float Inv4z = 0.5f / TMath::Sqrt(FourZSqr);
	//	return FourZSqr * Inv4z;
	//}
	//else
	//{
	//	float FourWSqr = Opr22 + Sum10;
	//	float Inv4w = 0.5f / TMath::Sqrt(FourWSqr);
	//	return (Data[0][1] - Data[1][0]) * Inv4w;
	//}
	
}


INLINE SVector2 SMatrix3::GetScale() const
{
	SVector2 Result;
	for (uint i = 0; i < Result.GetVectorSize(); ++i)
	{
		SVector2 Axis{ Data[i][0], Data[i][1] };
		Result[i] = SVector2::Distance(0.0f, Axis);
	}
	return Result;
}



// A 2D arrangement of numbers in rows and columns.
// This matrix type contains 4 rows and 4 columns.
// This matrix is mostly used for 3D math and has a lot more functions independent from the other matrix types.
struct SMatrix4 :public SMatrix<4, 4, float>
{
public:
	/// Constructors

	// Constructor, Default.
	SMatrix4()
	{}

	INLINE SMatrix4(SMatrix<4, 4, float> Other);

	// Constructor, Initiates the matrix with a location, rotation and scale.
	// @note - This constructor multiplies the matrices: Scale * Rotation * Translation.
	// @param Transform - The location, rotation and scale this matrix should be created by.
	INLINE SMatrix4(STransform Transform);

	// Constructor, Initiates the matrix with a location, rotation and scale.
	// @note - This constructor multiplies the matrices: Scale * Rotation * Translation.
	// @param Location - The location this matrix should be placed at.
	// @param Rotation - The rotation this matrix should be at.
	// @param Scale - The scale this matrix should be at.
	INLINE SMatrix4(SVector Location, SQuaternion Rotation, SVector Scale);



	/// Conversions

	// Converts this matrix to a transformation.
	INLINE STransform GetTransform() const;

#ifdef INCLUDE_DIRECTX_MATH
	// Converts this matrix to a DirectX::XMMATRIX
	INLINE DirectX::XMMATRIX ToXMMatrix() const;
#endif


	/// Functions

	// Moves this matrix by an inputted amount.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @return - The resulting matrix from the translation.
	INLINE SMatrix4 Translate(const SVector4& Translation) const;

	// Moves this matrix by an inputted amount.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @param W - The W value. 
	// @return - The resulting matrix from the translationg.
	INLINE SMatrix4 Translate(const SVector& Translation, const float& W = 1.0f) const;

	// Moves this matrix by an inputted amount.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param Z - The Z location to move in.
	// @param W - The W value. 
	// @return - The resulting matrix from the translationg.
	INLINE SMatrix4 Translate(const float& X, const float& Y, const float& Z, const float& W = 1.0f) const;

	// Moves this matrix by an inputted amount and sets this matrix to the returning result.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @return - The resulting matrix from the translation.
	INLINE SMatrix4 SetTranslate(const SVector4& Translation);

	// Moves this matrix by an inputted amount and sets this matrix to the returning result.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @param W - The W value. 
	// @return - The resulting matrix from the translationg.
	INLINE SMatrix4 SetTranslate(const SVector& Translation, const float& W = 1.0f);

	// Moves this matrix by an inputted amount and sets this matrix to the returning result.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param Z - The Z location to move in.
	// @param W - The W value. 
	// @return - The resulting matrix from the translationg.
	INLINE SMatrix4 SetTranslate(const float& X, const float& Y, const float& Z, const float& W = 1.0f);

	// Makes this matrix a translation matrix and represents the inputted location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @return - The translation matrix.
	INLINE SMatrix4 ToTranslation(const SVector4& Translation);

	// Makes this matrix a translation matrix and represents the inputted location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The location multiplied to this matrix.
	// @param W - The W value. 
	// @return - The translation matrix.
	INLINE SMatrix4 ToTranslation(const SVector& Translation, const float& W = 1.0f);

	// Makes this matrix a translation matrix and represents the inputted location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param Z - The Z location to move in.
	// @param W - The W value. 
	// @return - The translation matrix.
	INLINE SMatrix4 ToTranslation(const float& X, const float& Y, const float& Z, const float& W = 1.0f);

	// Rotates along the X axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 RotatePitch(const float& Angle) const;

	// Rotates along the Y axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 RotateYaw(const float& Angle) const;

	// Rotates along the Z axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 RotateRoll(const float& Angle) const;

	// Rotates along the X axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotatePitch(const float& Angle);

	// Rotates along the Y axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotateYaw(const float& Angle);

	// Rotates along the Z axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotateRoll(const float& Angle);

	// Makes this matrix a rotation matrix along the X axis and sets represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotationPitch(const float& Angle);

	// Makes this matrix a rotation matrix along the Y axis and sets represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotationYaw(const float& Angle);

	// Makes this matrix a rotation matrix along the Z axis and sets represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotationRoll(const float& Angle);

	// Rotates along the X axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 RotateX(const float& Angle) const;

	// Rotates along the Y axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 RotateY(const float& Angle) const;

	// Rotates along the Z axis by an inputted amount.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 RotateZ(const float& Angle) const;

	// Rotates along the X axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotateX(const float& Angle);

	// Rotates along the Y axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotateY(const float& Angle);

	// Rotates along the Z axis by an inputted amount and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotateZ(const float& Angle);

	// Makes this matrix a rotation matrix along the X axis and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotationX(const float& Angle);

	// Makes this matrix a rotation matrix along the X axis and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotationY(const float& Angle);

	// Makes this matrix a rotation matrix along the X axis and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotationZ(const float& Angle);

	// Rotates this matrix along all axis and returns the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 Rotate(const SQuaternion& Rotation) const;

	// Rotates this matrix along all axis and returns the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 Rotate(const SVector4& Rotation) const;

	// Rotates this matrix along all axis and returns the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @param W - The W component.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 Rotate(const SVector& Rotation) const;

	// Rotates this matrix along all axis and returns the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param X - The angle to rotate in the X axis.
	// @param Y - The angle to rotate in the Y axis.
	// @param Z - The angle to rotate in the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 Rotate(const float& Pitch, const float& Yaw, const float& Roll) const;

	// Rotates this matrix along all axis and returns the result and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotate(const SQuaternion& Rotation);

	// Rotates this matrix along all axis and returns the result and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotate(const SVector4& Rotation);

	// Rotates this matrix along all axis and returns the result and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @param W - The W component.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotate(const SVector& Rotation);

	// Rotates this matrix along all axis and returns the result and sets this matrix to the resulting matrix.
	// @note - This assumes the input is in Radians.
	// @param X - The angle to rotate in the X axis.
	// @param Y - The angle to rotate in the Y axis.
	// @param Z - The angle to rotate in the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 SetRotate(const float& Pitch, const float& Yaw, const float& Roll);

	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotation(const SQuaternion& Rotation);

	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotation(const SVector4& Rotation);

	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @param W - The W component.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotation(const SVector& Rotation);

	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param X - The angle to rotate in the X axis.
	// @param Y - The angle to rotate in the Y axis.
	// @param Z - The angle to rotate in the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the rotation.
	INLINE SMatrix4 ToRotation(const float& Pitch, const float& Yaw, const float& Roll);

	// Scales this matrix and returns the resulting matrix.
	// @param InScale - The scale this matrix will be multiplied by.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 Scale(const SVector4& InScale) const;

	// Scales this matrix and returns the resulting matrix.
	// @param InScale - The scale this matrix will be multiplied by.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 Scale(const SVector& InScale) const;

	// Scales this matrix and returns the resulting matrix.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @param Z - The scaling amount along the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 Scale(const float& X, const float& Y, const float& Z) const;

	// Scales this matrix and returns the resulting matrix and sets this matrix the the resulting matrix.
	// @param InScale - The scale this matrix will be multiplied by.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 SetScale(const SVector4& InScale);

	// Scales this matrix and returns the resulting matrix and sets this matrix the the resulting matrix.
	// @param InScale - The scale this matrix will be multiplied by.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 SetScale(const SVector& InScale);

	// Scales this matrix and returns the resulting matrix and sets this matrix the the resulting matrix.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @param Z - The scaling amount along the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 SetScale(const float& X, const float& Y, const float& Z);

	// Makes this matrix a scale matrix and represents the inputted scale.
	// @param InScale - The scale this matrix will be multiplied by.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 ToScale(const SVector4& NewScale);

	// Makes this matrix a scale matrix and represents the inputted scale.
	// @param InScale - The scale this matrix will be multiplied by.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 ToScale(const SVector& NewScale);

	// Makes this matrix a scale matrix and represents the inputted scale.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @param Z - The scaling amount along the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	INLINE SMatrix4 ToScale(const float& X, const float& Y, const float& Z);

	// Transforms this matrix using a location, rotation and scale.
	// @param InTransform - The location, rotation and scale to multiply against this matrix.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix4 Transform(const STransform& InTransform) const;

	// Transforms this matrix using a location, rotation and scale.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The 
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix4 Transform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale) const;

	// Transforms this matrix using a location, rotationa and scale and set this matrix to the resulting matrix.
	// @param InTransform - The location, rotation and scale to multiply against this matrix.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix4 SetTransform(const STransform& InTransform);

	// Transforms this matrix using a location, rotationa and scale and set this matrix to the resulting matrix.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The 
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix4 SetTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale);

	// Makes this matrix a transform matrix and represents the inputted transformation.
	// @param InTransform - The location, rotation and scale to multiply against this matrix.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix4 ToTransform(const STransform& NewTransform);

	// Makes this matrix a transform matrix and represents the inputted transformation.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The scale to be scaled in.
	// @return - The resulting matrix after the transformation.
	INLINE SMatrix4 ToTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale);

	// Transforms a vector by this matrix.
	// @param Vec - The vector to transform by.
	// @return - The resulting vector after the transformation.
	INLINE SVector4 VectorTransform(SVector4 Vec) const;

	// 
	// @param Other - 
	// @return - 
	INLINE SMatrix4 Merge(SMatrix4 Other);



	/// Getters


	// Returns this matrix's location.
	INLINE SVector GetLocation() const;

	// Returns this matrix's rotation.
	INLINE SQuaternion GetRotation() const;

	// Returns thsi matrix's scale.
	INLINE SVector GetScale() const;



	/// Statics

	// 
	// @param FovAngleY - 
	// @param AspectRatio - 
	// @param NearZ - The near clip plane.
	// @param FarZ - The far clip plane.
	// @return - 
	static INLINE SMatrix4 PerspectiveFOV(float FovAngleY, float AspectRatio, float NearZ, float FarZ)
	{
		assert(NearZ > 0.0f && FarZ > 0.0f);
		assert(!TMath::ScalarNearEqual(FovAngleY, 0.0f, 0.00001f * 2.0f));
		assert(!TMath::ScalarNearEqual(AspectRatio, 0.0f, 0.00001f));
		assert(!TMath::ScalarNearEqual(FarZ, NearZ, 0.00001f));

		float SinFOV;
		float CosFOV;
		TMath::SinCos(&SinFOV, &CosFOV, 0.5f * FovAngleY);

		float Height = CosFOV / SinFOV;
		float Width = Height / AspectRatio;
		float FRange = FarZ / (FarZ - NearZ);

		SMatrix4 Result{ 0.0f };
		Result[0][0] = Width;
		Result[1][1] = Height;
		Result[2][2] = FRange;
		Result[2][3] = 1.0f;
		Result[3][2] = -FRange * NearZ;

		return Result;
	}


	// 
	// @param EyePosition - 
	// @param FocusPosition - 
	// @param UpDirection - 
	// @return - 
	static INLINE SMatrix4 LookAt(SVector4 EyePosition, SVector4 FocusPosition, SVector4 UpDirection)
	{
		SVector4 EyeDirection{ FocusPosition - EyePosition };
		return LookTo(EyePosition, EyeDirection, UpDirection);
	}


	// 
	// @param EyePosition - 
	// @param EyeDirection - 
	// @param UpDirection - 
	// @return - 
	static INLINE SMatrix4 LookTo(SVector4 EyePosition, SVector4 EyeDirection, SVector4 UpDirection)
	{
		assert(!(EyeDirection == 0.0f));
		assert(!(UpDirection == 0.0f));

		SVector4 R2{ EyeDirection.GetNormal() };
		SVector4 R0{ UpDirection | R2 };

		R0.Normalize();

		SVector4 R1{ SVector4::CrossProduct(R2, R0) };
		SVector4 NegEyePos{ -EyePosition };

		SVector4 D0 = SVector4::DotProduct(R0, NegEyePos);
		SVector4 D1 = SVector4::DotProduct(R1, NegEyePos);
		SVector4 D2 = SVector4::DotProduct(R2, NegEyePos);

		STVector<4, bool> Control{ false, false, false, true };

		SMatrix4 Mat;
		Mat.SetColumn(0, SVector4::Select(D0, R0, Control));
		Mat.SetColumn(1, SVector4::Select(D1, R1, Control));
		Mat.SetColumn(2, SVector4::Select(D2, R2, Control));
		Mat.SetColumn(3, SVector4{ 0.0f, 0.0f, 0.0f, 1.0f });

		Mat.Merge(Mat);
		return Mat;
	}


	// 
	// @param Other - 
	// @return -
	static INLINE SMatrix4 Transpose(SMatrix4 Other)
	{
		SMatrix4 Mat;
		Mat[0] = SVector4::Merge(Other[0], Other[2], EAxis::X, EAxis::Y);
		Mat[1] = SVector4::Merge(Other[1], Other[3], EAxis::X, EAxis::Y);
		Mat[2] = SVector4::Merge(Other[0], Other[2], EAxis::Z, EAxis::W);
		Mat[3] = SVector4::Merge(Other[1], Other[3], EAxis::Z, EAxis::W);

		SMatrix4 Result;
		Result[0] = SVector4::Merge(Mat[0], Mat[1], EAxis::X, EAxis::Y);
		Result[1] = SVector4::Merge(Mat[0], Mat[1], EAxis::Z, EAxis::W);
		Result[2] = SVector4::Merge(Mat[2], Mat[3], EAxis::X, EAxis::Y);
		Result[3] = SVector4::Merge(Mat[2], Mat[3], EAxis::Z, EAxis::W);
		return Result;
	}


	// Returns a matrix translated to a location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The translation to be moved by.
	// @return - The translation matrix.
	static INLINE SMatrix4 MatrixTranslate(SVector4 Translation)
	{
		SMatrix4 Result;
		return Result.ToTranslation(Translation);
	}


	// Returns a matrix translated to a location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param Translation - The translation to be moved by.
	// @param W - The W value.
	// @return - The translation matrix.
	static INLINE SMatrix4 MatrixTranslate(const SVector& Translation, const float& W = 1.0f)
	{
		SMatrix4 Result;
		return Result.ToTranslation(Translation, W);
	}


	// Returns a matrix translated to a location.
	// @note - If the W component is set to 0, this matrix will be set to the location instead of a movement.
	// @param X - The X location to move in.
	// @param Y - The Y location to move in.
	// @param Z - The Z location to move in.
	// @param W - The W value. 
	// @return - The translation matrix.
	static INLINE SMatrix4 MatrixTranslate(const float& X, const float& Y, const float& Z, const float& W = 1.0f)
	{
		SMatrix4 Result;
		return Result.ToTranslation(X, Y, Z, W);
	}


	// Returns a matrix rotated in the X axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotatePitch(const float& Angle)
	{
		SMatrix4 Result;
		return Result.ToRotationPitch(Angle);
	}


	// Returns a matrix rotated in the Y axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotateYaw(const float& Angle)
	{
		SMatrix4 Result;
		return Result.ToRotationYaw(Angle);
	}


	// Returns a matrix rotated in the Z axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotateRoll(const float& Angle)
	{
		SMatrix4 Result;
		return Result.ToRotationRoll(Angle);
	}


	// Returns a matrix rotated in the X axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotateX(const float& Angle)
	{
		SMatrix4 Result;
		return Result.ToRotationPitch(Angle);
	}


	// Returns a matrix rotated in the Y axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotateY(const float& Angle)
	{
		SMatrix4 Result;
		return Result.ToRotationYaw(Angle);
	}


	// Returns a matrix rotated in the Z axis by an inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Angle - The angle to rotate by.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotateZ(const float& Angle)
	{
		SMatrix4 Result;
		return Result.ToRotationRoll(Angle);
	}


	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotate(const SQuaternion& Rotation)
	{
		SMatrix4 Result;
		return Result.ToRotation(Rotation);
	}


	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotate(const SVector4& Rotation)
	{
		SMatrix4 Result;
		return Result.ToRotation(Rotation);
	}


	// Makes this matrix a rotation matrix and represents the inputted angle.
	// @note - This assumes the input is in Radians.
	// @param Rotation - The rotations in the X, Y and Z axis' to rotate.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotate(const SVector& Rotation)
	{
		SMatrix4 Result;
		return Result.ToRotation(Rotation);
	}


	// Returns a matrix rotated by an angle in all dimensions.
	// @note - This assumes the input is in Radians.
	// @param X - The angle to rotate in the X axis.
	// @param Y - The angle to rotate in the Y axis.
	// @param Z - The angle to rotate in the Z axis.
	// @return - The resulting matrix from the rotation.
	static INLINE SMatrix4 MatrixRotate(const float& X, const float& Y, const float& Z)
	{
		SMatrix4 Result;
		return Result.ToRotation(X, Y, Z);
	}


	// Returns a matrix scaled by a value in all dimensions.
	// @param InScale - The amount to be scaled.
	// @return - The resulting matrix from the scaling.
	static INLINE SMatrix4 MatrixScale(const SVector4& Scale)
	{
		SMatrix4 Result;
		return Result.ToScale(Scale);
	}


	// Returns a matrix scaled by a value in all dimensions.
	// @param InScale - The amount to be scaled.
	// @return - The resulting matrix from the scaling.
	static INLINE SMatrix4 MatrixScale(const SVector& Scale)
	{
		SMatrix4 Result;
		return Result.ToScale(Scale);
	}


	// Returns a matrix scaled by a value in all dimensions.
	// @param X - The scaling amount along the X axis.
	// @param Y - The scaling amount along the Y axis.
	// @param Z - The scaling amount along the Z axis.
	// @param W - The W component.
	// @return - The resulting matrix from the scaling.
	static INLINE SMatrix4 MatrixScale(const float& X, const float& Y, const float& Z)
	{
		SMatrix4 Result;
		return Result.ToScale(X, Y, Z);
	}


	// Returns a matrix scaled, rotated and translated by an inputted transformation.
	// @param Transform - The location, rotation and scale this matrix will be transformed.
	// @return - The resulting matrix after the transformation.
	static INLINE SMatrix4 MatrixTransform(const STransform& Transform)
	{
		SMatrix4 Result;
		return Result.ToTransform(Transform);
	}


	// Returns a matrix scaled, rotated and translated by an inputted transformation.
	// @param Location - The location to be moved in.
	// @param Rotation - The rotation to be rotated in.
	// @param InScale - The scale to be scaled in.
	// @return - The resulting matrix after the transformation.
	static INLINE SMatrix4 MatrixTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& Scale)
	{
		SMatrix4 Result;
		return Result.ToTransform(Location, Rotation, Scale);
	}


	// Gets the world transform as a matrix based on an inputted transform.
	// @param Transform - The object to get the world location of.
	// @return - Returns the world transform as a matrix4.
	static INLINE SMatrix4 GetWorldTransform(const STransform Transform)
	{
		SMatrix4 ParentMat;
		ParentMat.Identity();
		if (Transform.GetParent()) ParentMat = GetWorldTransform(*Transform.GetParent());

		SMatrix4 ScaleMat;
		SMatrix4 RotationMat;
		SMatrix4 LocationMat;

		ScaleMat.ToScale(Transform.Scale);
		RotationMat.ToRotation(Transform.Rotation);
		LocationMat.ToTranslation(Transform.Location);

		SMatrix4 Local{ ScaleMat * RotationMat * LocationMat };
		//SMatrix4 Local{ LocationMat * RotationMat * ScaleMat };

		//return ParentMat * Local;
		return Local * ParentMat;
	}
};


INLINE SMatrix4::SMatrix4(SMatrix<4, 4, float> Other)
{
	for (uint y = 0; y < 4; ++y)
	{
		for (uint x = 0; x < 4; ++x)
		{
			Data[y][x] = Other[y][x];
		}
	}
}


INLINE SMatrix4::SMatrix4(STransform Transform)
{
	ToTransform(Transform);
}


INLINE SMatrix4::SMatrix4(SVector Location, SQuaternion Rotation, SVector Scale)
{
	ToTransform(Location, Rotation, Scale);
}


INLINE STransform SMatrix4::GetTransform() const
{
	return STransform{ GetLocation(), GetRotation(), GetScale() };
}


#ifdef INCLUDE_DIRECTX_MATH
INLINE DirectX::XMMATRIX SMatrix4::ToXMMatrix() const
{
	DirectX::XMMATRIX Result;
#ifdef _XM_NO_INTRINSICS_
	for (uint y = 0; y < 4; ++y)
	{
		for (uint x = 0; x < 4; ++x)
		{
			Result(y, x) = Data[y][x];
		}
	}

#else
	for (uint i = 0; i < 4; ++i)
	{
		Result.r[i] = Data[i].ToXMVector();
	}
#endif // _XM_NO_INTRINSICS_
	return Result;

}
#endif


INLINE SMatrix4 SMatrix4::Translate(const SVector4& Translation) const
{
	SMatrix4 Result;
	Result.Identity();
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Result[i][3] = Translation[i];
	}
	return *this * Result;
}


INLINE SMatrix4 SMatrix4::Translate(const SVector& Translation, const float& W) const
{
	return Translate(SVector4{ Translation, W });
}


INLINE SMatrix4 SMatrix4::Translate(const float& X, const float& Y, const float& Z, const float& W) const
{
	return Translate(SVector4{ X, Y, Z, W });
}


INLINE SMatrix4 SMatrix4::SetTranslate(const SVector4& Translation)
{
	*this = Translate(Translation);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetTranslate(const SVector& Translation, const float& W)
{
	*this = Translate(Translation, W);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetTranslate(const float& X, const float& Y, const float& Z, const float& W)
{
	*this = Translate(X, Y, Z, W);
	return *this;
}


INLINE SMatrix4 SMatrix4::ToTranslation(const SVector4& Translation)
{
	Identity();
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Data[i][3] = Translation[i];
	}
	return *this;
}


INLINE SMatrix4 SMatrix4::ToTranslation(const SVector& Translation, const float& W)
{
	return ToTranslation(SVector4{ Translation, W });
}


INLINE SMatrix4 SMatrix4::ToTranslation(const float& X, const float& Y, const float& Z, const float& W)
{
	return ToTranslation(SVector4{ X, Y, Z, W });
}


INLINE SMatrix4 SMatrix4::RotatePitch(const float& Angle) const
{
	SMatrix4 Result;
	Result.Identity();

	Result[1][1] = TMath::Cos(Angle);
	Result[1][2] = TMath::Sin(Angle);
	Result[2][1] = -TMath::Sin(Angle);
	Result[2][2] = TMath::Cos(Angle);

	return Result;
}


INLINE SMatrix4 SMatrix4::RotateYaw(const float& Angle) const
{
	SMatrix4 Result;
	Result.Identity();

	Result[0][0] = TMath::Cos(Angle);
	Result[0][2] = -TMath::Sin(Angle);
	Result[2][0] = TMath::Sin(Angle);
	Result[2][2] = TMath::Cos(Angle);

	return Result;
}


INLINE SMatrix4 SMatrix4::RotateRoll(const float& Angle) const
{
	SMatrix4 Result;
	Result.Identity();

	Result[0][0] = TMath::Cos(Angle);
	Result[0][1] = TMath::Sin(Angle);
	Result[1][0] = -TMath::Sin(Angle);
	Result[1][1] = TMath::Cos(Angle);

	return Result;
}


INLINE SMatrix4 SMatrix4::SetRotatePitch(const float& Angle)
{
	*this = RotatePitch(Angle);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetRotateYaw(const float& Angle)
{
	*this = RotateYaw(Angle);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetRotateRoll(const float& Angle)
{
	*this = RotateRoll(Angle);
	return *this;
}


INLINE SMatrix4 SMatrix4::ToRotationPitch(const float& Angle)
{
	Identity();
	Data[1][1] = TMath::Cos(Angle);
	Data[1][2] = TMath::Sin(Angle);
	Data[2][1] = -TMath::Sin(Angle);
	Data[2][2] = TMath::Cos(Angle);

	return *this;
}


INLINE SMatrix4 SMatrix4::ToRotationYaw(const float& Angle)
{
	Identity();
	Data[0][0] = TMath::Cos(Angle);
	Data[0][2] = -TMath::Sin(Angle);
	Data[2][0] = TMath::Sin(Angle);
	Data[2][2] = TMath::Cos(Angle);

	return *this;
}


INLINE SMatrix4 SMatrix4::ToRotationRoll(const float& Angle)
{
	Identity();
	Data[0][0] = TMath::Cos(Angle);
	Data[0][1] = TMath::Sin(Angle);
	Data[1][0] = -TMath::Sin(Angle);
	Data[1][1] = TMath::Cos(Angle);

	return *this;
}


INLINE SMatrix4 SMatrix4::RotateX(const float& Angle) const
{
	return RotatePitch(Angle);
}


INLINE SMatrix4 SMatrix4::RotateY(const float& Angle) const
{
	return RotateYaw(Angle);
}


INLINE SMatrix4 SMatrix4::RotateZ(const float& Angle) const
{
	return RotateRoll(Angle);
}


INLINE SMatrix4 SMatrix4::SetRotateX(const float& Angle)
{
	return SetRotatePitch(Angle);
}


INLINE SMatrix4 SMatrix4::SetRotateY(const float& Angle)
{
	return SetRotateYaw(Angle);
}


INLINE SMatrix4 SMatrix4::SetRotateZ(const float& Angle)
{
	return SetRotateRoll(Angle);
}


INLINE SMatrix4 SMatrix4::ToRotationX(const float& Angle)
{
	return ToRotationPitch(Angle);
}


INLINE SMatrix4 SMatrix4::ToRotationY(const float& Angle)
{
	return ToRotationYaw(Angle);
}


INLINE SMatrix4 SMatrix4::ToRotationZ(const float& Angle)
{
	return ToRotationRoll(Angle);
}


INLINE SMatrix4 SMatrix4::Rotate(const SQuaternion& Rotation) const
{
	SMatrix4 Result;
	Result = RotatePitch(Rotation.X);
	Result *= RotateYaw(Rotation.Y);
	Result *= RotateRoll(Rotation.Z);

	return *this * Result;
}


INLINE SMatrix4 SMatrix4::Rotate(const SVector4& Rotation) const
{
	return Rotate(SQuaternion{ Rotation });
}


INLINE SMatrix4 SMatrix4::Rotate(const SVector& Rotation) const
{
	return Rotate(SQuaternion{ Rotation, 1.0f });
}


INLINE SMatrix4 SMatrix4::Rotate(const float& Pitch, const float& Yaw, const float& Roll) const
{
	return Rotate(SQuaternion{ Pitch, Yaw, Roll, 1.0f });
}


INLINE SMatrix4 SMatrix4::SetRotate(const SQuaternion& Rotation)
{
	*this = Rotate(Rotation);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetRotate(const SVector4& Rotation)
{
	*this = Rotate(Rotation);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetRotate(const SVector& Rotation)
{
	*this = Rotate(Rotation);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetRotate(const float& Pitch, const float& Yaw, const float& Roll)
{
	*this = Rotate(Pitch, Yaw, Roll);
	return *this;
}


INLINE SMatrix4 SMatrix4::ToRotation(const SQuaternion& Rotation)
{
	SMatrix4 Matrix;
	Matrix = RotatePitch(Rotation.X);
	Matrix *= RotateYaw(Rotation.Y);
	Matrix *= RotateRoll(Rotation.Z);
	*this = Matrix;
	return *this;
}


INLINE SMatrix4 SMatrix4::ToRotation(const SVector4& Rotation)
{
	return ToRotation(SQuaternion{ Rotation });
}


INLINE SMatrix4 SMatrix4::ToRotation(const SVector& Rotation)
{
	return ToRotation(SQuaternion{ Rotation });
}


INLINE SMatrix4 SMatrix4::ToRotation(const float& Pitch, const float& Yaw, const float& Roll)
{
	return ToRotation(SQuaternion{ Pitch, Yaw, Roll });
}


INLINE SMatrix4 SMatrix4::Scale(const SVector4& InScale) const
{
	SMatrix4 Result{ 0.0f };
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Result[i][i] = InScale[i];
	}
	return *this * Result;
}


INLINE SMatrix4 SMatrix4::Scale(const SVector& InScale) const
{
	return Scale(SVector4{ InScale });
}


INLINE SMatrix4 SMatrix4::Scale(const float& X, const float& Y, const float& Z) const
{
	return Scale(SVector{ X, Y, Z });
}


INLINE SMatrix4 SMatrix4::SetScale(const SVector4& InScale)
{
	*this = Scale(InScale);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetScale(const SVector& InScale)
{
	*this = Scale(InScale);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetScale(const float& X, const float& Y, const float& Z)
{
	*this = Scale(X, Y, Z);
	return *this;
}


INLINE SMatrix4 SMatrix4::ToScale(const SVector4& NewScale)
{
	Identity();
	for (uint i = 0; i < GetRowCount(); ++i)
	{
		Data[i][i] = NewScale[i];
	}
	Data[3][3] = 1.0f;
	return *this;
}


INLINE SMatrix4 SMatrix4::ToScale(const SVector& NewScale)
{
	return ToScale(SVector4{ NewScale });
}


INLINE SMatrix4 SMatrix4::ToScale(const float& X, const float& Y, const float& Z)
{
	return ToScale(SVector{ X, Y, Z });
}


INLINE SMatrix4 SMatrix4::Transform(const STransform& InTransform) const
{
	SMatrix4 Result{ Scale(InTransform.Scale) * Rotate(InTransform.Rotation) * Translate(InTransform.Location) };
	return *this * Result;
}


INLINE SMatrix4 SMatrix4::Transform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale) const
{
	return Transform(STransform{ Location, Rotation, InScale });
}


INLINE SMatrix4 SMatrix4::SetTransform(const STransform& InTransform)
{
	*this = Transform(InTransform);
	return *this;
}


INLINE SMatrix4 SMatrix4::SetTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale)
{
	*this = Transform(Location, Rotation, InScale);
	return *this;
}


INLINE SMatrix4 SMatrix4::ToTransform(const STransform& NewTransform)
{
	SMatrix4 TranslateMat;
	SMatrix4 RotationMat;
	SMatrix4 ScaleMat;

	TranslateMat.ToTranslation(NewTransform.Location);
	RotationMat.ToRotation(NewTransform.Rotation);
	ScaleMat.ToScale(NewTransform.Scale);

	*this = ScaleMat * RotationMat * TranslateMat;
	return *this;
}


INLINE SMatrix4 SMatrix4::ToTransform(const SVector& Location, const SQuaternion& Rotation, const SVector& InScale)
{
	return ToTransform(STransform{ Location, Rotation, InScale });
}


INLINE SVector4 SMatrix4::VectorTransform(SVector4 Vec) const
{
	SVector4 VecX{ Vec[0] };
	SVector4 VecY{ Vec[1] };
	SVector4 VecZ{ Vec[2] };

	SVector4 Result{ (VecZ * Data[2]) + Data[3] };
	Result = (VecY * Data[1]) + Result;
	Result = (VecX * Data[0]) + Result;
	return Result;
}


INLINE SMatrix4 SMatrix4::Merge(SMatrix4 Mat)
{
	SMatrix4 P;
	P[0] = SVector4::Merge(Mat[0], Mat[2], EAxis::X, EAxis::Y);
	P[1] = SVector4::Merge(Mat[1], Mat[3], EAxis::X, EAxis::Y);
	P[2] = SVector4::Merge(Mat[0], Mat[2], EAxis::Z, EAxis::W);
	P[3] = SVector4::Merge(Mat[1], Mat[3], EAxis::Z, EAxis::W);


	SMatrix4 MT;
	MT[0] = SVector4::Merge(P[0], P[1], EAxis::X, EAxis::Y);
	MT[1] = SVector4::Merge(P[0], P[1], EAxis::Z, EAxis::W);
	MT[2] = SVector4::Merge(P[2], P[3], EAxis::X, EAxis::Y);
	MT[3] = SVector4::Merge(P[2], P[3], EAxis::Z, EAxis::W);
	return MT;
}


INLINE SVector SMatrix4::GetLocation() const
{
	return GetColumn(3);
}


INLINE SQuaternion SMatrix4::GetRotation() const
{


	/*SQuaternion Result;
	float Trace{ Data[0][0] + Data[1][1] + Data[2][2] };
	if (Trace > 0.0f)
	{
		float S{ 0.5f / TMath::Sqrt(Trace + 1.0f) };
		Result.W = 0.25f / S;
		Result.X = (Data[2][1] - Data[1][2]) * S;
		Result.Y = (Data[0][2] - Data[2][0]) * S;
		Result.Z = (Data[1][0] - Data[0][1]) * S;
	}
	else
	{
		if (Data[0][0] > Data[1][1] && Data[0][0] > Data[2][2])
		{
			float S{ 2.0f * TMath::Sqrt(1.0f + Data[0][0] - Data[1][1] - Data[2][2]) };
			Result.W = (Data[2][1] - Data[1][0]) / S;
			Result.X = 0.25f * S;
			Result.Y = (Data[0][1] + Data[1][0]) / S;
			Result.Z = (Data[1][0] + Data[1][0]) / S;
		}
		else if (Data[1][1] > Data[2][2])
		{
			float S{ 2.0f * TMath::Sqrt(1.0f + Data[1][1] - Data[0][0] - Data[2][2]) };
			Result.W = (Data[0][2] - Data[2][0]) / S;
			Result.X = (Data[0][1] + Data[1][0]) / S;
			Result.Y = 0.25f * S;
			Result.Z = (Data[1][2] + Data[2][1]) / S;
		}
		else
		{
			float S{ 2.0f * TMath::Sqrt(1.0f + Data[2][2] - Data[0][0] - Data[1][1]) };
			Result.W = (Data[1][0] - Data[0][1]) / S;
			Result.X = (Data[0][2] + Data[2][0]) / S;
			Result.Y = (Data[1][2] + Data[2][1]) / S;
			Result.Z = 0.25f * S;
		}
	}
	return Result;*/

	// If the matrix is NULL, return Identity quaternion. Additionally if any of the Axis equals 0, then a rotation can't be made.
	/*if (GetScaledAxis(EAxis::X).IsZero() || GetScaledAxis(EAxis::Y).IsZero() || GetScaledAxis(EAxis::Z).IsZero())
	{
		return SQuaternion::Identity;
	}

	if ((TMath::Abs(1.0f - GetScaledAxis(EAxis::X).SizeSquared()) <= SMALL_NUMBER) && (TMath::Abs(1.0f - GetScaledAxis(EAxis::Y).SizeSquared()) <= SMALL_NUMBER) && TMath::Abs(1.0f - GetScaledAxis(EAxis::Z).SizeSquared()) <= SMALL_NUMBER)
	{
		return SQuaternion::Identity;
	}

	SQuaternion Result;
	float S;

	const float Trace{ Data[0][0] + Data[1][1] + Data[2][2] };
	if (Trace > 0.0f)
	{
		float InvS{ TMath::InvSqrt(Trace + 1.0f) };
		Result.W = 0.5f * (1.0f / InvS);
		S = 0.5f * InvS;

		Result.X = (Data[1][2] - Data[2][1]) * S;
		Result.Y = (Data[2][0] - Data[0][2]) * S;
		Result.Z = (Data[0][1] - Data[1][0]) * S;
	}
	else
	{
		int32 i{ 0 };

		if (Data[1][1] > Data[0][0]) i = 1;
		if (Data[2][2] > Data[i][i]) i = 2;

		static const int32 Next[3] = { 1, 2, 3 };
		const int32 j = Next[i];
		const int32 k = Next[j];

		S = Data[i][i] - Data[j][j] - Data[k][k] + 1.0f;

		float InvS{ TMath::InvSqrt(S) };
		float QT[4];
		QT[i] = 0.5f * (1.0f / InvS);
		S = 0.5f * InvS;

		QT[3] = (Data[j][k] - Data[k][j]) * S;
		QT[j] = (Data[i][j] + Data[j][i]) * S;
		QT[k] = (Data[i][k] + Data[k][i]) * S;

		Result.X = QT[0];
		Result.Y = QT[1];
		Result.Z = QT[2];
		Result.W = QT[3];
		Result.CheckNaN();
	}
	return Result;*/

	// http://www.iri.upc.edu/files/scidoc/2068-Accurate-Computation-of-Quaternions-from-Rotation-Matrices.pdf

	SQuaternion Quat;
	float R22{ Data[2][2] };
	if (R22 <= 0.0f)
	{
		float Dif10{ Data[1][1] - Data[0][0] };
		float Omr22{ 1.0f - R22 };
		if (Dif10 <= 0.0f)
		{
			float FourXSqr{ Omr22 - Dif10 };
			float Inv4x{ 0.5f / TMath::Sqrt(FourXSqr) };
			Quat.X = FourXSqr * Inv4x;
			Quat.Y = (Data[0][1] + Data[1][0]) * Inv4x;
			Quat.Z = (Data[0][2] + Data[2][0]) * Inv4x;
			Quat.W = (Data[1][2] - Data[2][1]) * Inv4x;
		}
		else
		{
			float FourYSqr{ Omr22 + Dif10 };
			float Inv4y{ 0.5f / TMath::Sqrt(FourYSqr) };
			Quat.X = (Data[0][1] + Data[1][0]) * Inv4y;
			Quat.Y = FourYSqr * Inv4y;
			Quat.Z = (Data[1][2] + Data[2][1]) * Inv4y;
			Quat.W = (Data[2][0] - Data[0][2]) * Inv4y;
		}
	}
	else
	{
		float Sum10{ Data[1][1] + Data[0][0] };
		float Opr22{ 1.0f + R22 };
		if (Sum10 <= 0.0f)
		{
			float FourZSqr = Opr22 - Sum10;
			float Inv4z = 0.5f / TMath::Sqrt(FourZSqr);
			Quat.X = (Data[0][2] + Data[2][0]) * Inv4z;
			Quat.Y = (Data[1][2] + Data[2][1]) * Inv4z;
			Quat.Z = FourZSqr * Inv4z;
			Quat.W = (Data[0][1] - Data[1][0]) * Inv4z;
		}
		else
		{
			float FourWSqr{ Opr22 + Sum10 };
			float Inv4w{ 0.5f / sqrtf(FourWSqr) };
			Quat.X = (Data[1][2] - Data[2][1]) * Inv4w;
			Quat.Y = (Data[2][0] - Data[0][2]) * Inv4w;
			Quat.Z = (Data[0][1] - Data[1][0]) * Inv4w;
			Quat.W = FourWSqr * Inv4w;
		}
	}
	return Quat;
}


INLINE SVector SMatrix4::GetScale() const
{
	SVector Result;
	for (uint i = 0; i < Result.GetSize(); ++i)
	{
		SVector Axis{ Data[i][0], Data[i][1], Data[i][2] };
		Result[i] = SVector::Distance(0.0f, Axis);
	}
	return Result;
}