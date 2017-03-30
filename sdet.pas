{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit sdet;
interface
uses Math, Sysutils, Ap, ldlt;

function SMatrixLDLTDet(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function SMatrixDet(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function DeterminantLDLT(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function DeterminantSymmetric(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;

implementation

(*************************************************************************
Determinant calculation of the matrix given by LDLT decomposition.

Input parameters:
    A       -   LDLT-decomposition of the matrix,
                output of subroutine SMatrixLDLT.
    Pivots  -   table of permutations which were made during
                LDLT decomposition, output of subroutine SMatrixLDLT.
    N       -   size of matrix A.
    IsUpper -   matrix storage format. The value is equal to the input
                parameter of subroutine SMatrixLDLT.

Result:
    matrix determinant.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
function SMatrixLDLTDet(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    K : AlglibInteger;
begin
    Result := 1;
    if IsUpper then
    begin
        K := 0;
        while K<N do
        begin
            if Pivots[K]>=0 then
            begin
                Result := Result*A[K,K];
                K := K+1;
            end
            else
            begin
                Result := Result*(A[K,K]*A[K+1,K+1]-A[K,K+1]*A[K,K+1]);
                K := K+2;
            end;
        end;
    end
    else
    begin
        K := N-1;
        while K>=0 do
        begin
            if Pivots[K]>=0 then
            begin
                Result := Result*A[K,K];
                K := K-1;
            end
            else
            begin
                Result := Result*(A[K-1,K-1]*A[K,K]-A[K,K-1]*A[K,K-1]);
                K := K-2;
            end;
        end;
    end;
end;


(*************************************************************************
Determinant calculation of the symmetric matrix

Input parameters:
    A       -   matrix. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper = True, then symmetric matrix A is given by its
                upper triangle, and the lower triangle isn’t used by
                subroutine. Similarly, if IsUpper = False, then A is given
                by its lower triangle.

Result:
    determinant of matrix A.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
function SMatrixDet(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    SMatrixLDLT(A, N, IsUpper, Pivots);
    Result := SMatrixLDLTDet(A, Pivots, N, IsUpper);
end;


function DeterminantLDLT(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    K : AlglibInteger;
begin
    Result := 1;
    if IsUpper then
    begin
        K := 1;
        while K<=N do
        begin
            if Pivots[K]>0 then
            begin
                Result := Result*A[K,K];
                K := K+1;
            end
            else
            begin
                Result := Result*(A[K,K]*A[K+1,K+1]-A[K,K+1]*A[K,K+1]);
                K := K+2;
            end;
        end;
    end
    else
    begin
        K := N;
        while K>=1 do
        begin
            if Pivots[K]>0 then
            begin
                Result := Result*A[K,K];
                K := K-1;
            end
            else
            begin
                Result := Result*(A[K-1,K-1]*A[K,K]-A[K,K-1]*A[K,K-1]);
                K := K-2;
            end;
        end;
    end;
end;


function DeterminantSymmetric(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    LDLTDecomposition(A, N, IsUpper, Pivots);
    Result := DeterminantLDLT(A, Pivots, N, IsUpper);
end;


end.