{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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
unit rcond;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve;

function RMatrixRCond1(A : TReal2DArray; N : AlglibInteger):Double;
function RMatrixRCondInf(A : TReal2DArray; N : AlglibInteger):Double;
function SPDMatrixRCond(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function RMatrixTRRCond1(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
function RMatrixTRRCondInf(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
function HPDMatrixRCond(A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function CMatrixRCond1(A : TComplex2DArray; N : AlglibInteger):Double;
function CMatrixRCondInf(A : TComplex2DArray; N : AlglibInteger):Double;
function RMatrixLURCond1(const LUA : TReal2DArray; N : AlglibInteger):Double;
function RMatrixLURCondInf(const LUA : TReal2DArray; N : AlglibInteger):Double;
function SPDMatrixCholeskyRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function HPDMatrixCholeskyRCond(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function CMatrixLURCond1(const LUA : TComplex2DArray;
     N : AlglibInteger):Double;
function CMatrixLURCondInf(const LUA : TComplex2DArray;
     N : AlglibInteger):Double;
function CMatrixTRRCond1(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
function CMatrixTRRCondInf(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
function RCondThreshold():Double;

implementation

procedure RMatrixRCondTRInternal(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OneNorm : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure CMatrixRCondTRInternal(const A : TComplex2DArray;
     const N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OneNorm : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure SPDMatrixRCondCholeskyInternal(const CHA : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure HPDMatrixRCondCholeskyInternal(const CHA : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure RMatrixRCondLUInternal(const LUA : TReal2DArray;
     N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure CMatrixRCondLUInternal(const LUA : TComplex2DArray;
     const N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure RMatrixEstimateNorm(N : AlglibInteger;
     var V : TReal1DArray;
     var X : TReal1DArray;
     var ISGN : TInteger1DArray;
     var EST : Double;
     var KASE : AlglibInteger);forward;
procedure CMatrixEstimateNorm(const N : AlglibInteger;
     var V : TComplex1DArray;
     var X : TComplex1DArray;
     var EST : Double;
     var KASE : AlglibInteger;
     var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray);forward;
function InternalComplexRCondSCSUM1(const X : TComplex1DArray;
     N : AlglibInteger):Double;forward;
function InternalComplexRCondICMAX1(const X : TComplex1DArray;
     N : AlglibInteger):AlglibInteger;forward;
procedure InternalComplexRCondSaveAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);forward;
procedure InternalComplexRCondLoadAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);forward;


(*************************************************************************
Estimate of a matrix condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function RMatrixRCond1(A : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
    T : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert(N>=1, 'RMatrixRCond1: N<1!');
    SetLength(T, N);
    I:=0;
    while I<=N-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[J] := T[J]+AbsReal(A[I,J]);
            Inc(J);
        end;
        Inc(I);
    end;
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        Nrm := Max(Nrm, T[I]);
        Inc(I);
    end;
    RMatrixLU(A, N, N, Pivots);
    RMatrixRCondLUInternal(A, N, True, True, Nrm, V);
    Result := V;
end;


(*************************************************************************
Estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function RMatrixRCondInf(A : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert(N>=1, 'RMatrixRCondInf: N<1!');
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := 0;
        J:=0;
        while J<=N-1 do
        begin
            V := V+AbsReal(A[I,J]);
            Inc(J);
        end;
        Nrm := Max(Nrm, V);
        Inc(I);
    end;
    RMatrixLU(A, N, N, Pivots);
    RMatrixRCondLUInternal(A, N, False, True, Nrm, V);
    Result := V;
end;


(*************************************************************************
Condition number estimate of a symmetric positive definite matrix.

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm of condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   symmetric positive definite matrix which is given by its
                upper or lower triangle depending on the value of
                IsUpper. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A)), if matrix A is positive definite,
   -1, if matrix A is not positive definite, and its condition number
    could not be found by this algorithm.

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function SPDMatrixRCond(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    V : Double;
    Nrm : Double;
    T : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    SetLength(T, N);
    I:=0;
    while I<=N-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        J:=J1;
        while J<=J2 do
        begin
            if I=J then
            begin
                T[I] := T[I]+AbsReal(A[I,I]);
            end
            else
            begin
                T[I] := T[I]+AbsReal(A[I,J]);
                T[J] := T[J]+AbsReal(A[I,J]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        Nrm := Max(Nrm, T[I]);
        Inc(I);
    end;
    if SPDMatrixCholesky(A, N, IsUpper) then
    begin
        SPDMatrixRCondCholeskyInternal(A, N, IsUpper, True, Nrm, V);
        Result := V;
    end
    else
    begin
        Result := -1;
    end;
end;


(*************************************************************************
Triangular matrix: estimate of a condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A       -   matrix. Array[0..N-1, 0..N-1].
    N       -   size of A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function RMatrixTRRCond1(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
    T : TReal1DArray;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    Assert(N>=1, 'RMatrixTRRCond1: N<1!');
    SetLength(T, N);
    I:=0;
    while I<=N-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I+1;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I-1;
        end;
        J:=J1;
        while J<=J2 do
        begin
            T[J] := T[J]+AbsReal(A[I,J]);
            Inc(J);
        end;
        if IsUnit then
        begin
            T[I] := T[I]+1;
        end
        else
        begin
            T[I] := T[I]+AbsReal(A[I,I]);
        end;
        Inc(I);
    end;
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        Nrm := Max(Nrm, T[I]);
        Inc(I);
    end;
    RMatrixRCondTRInternal(A, N, IsUpper, IsUnit, True, Nrm, V);
    Result := V;
end;


(*************************************************************************
Triangular matrix: estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function RMatrixTRRCondInf(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    Assert(N>=1, 'RMatrixTRRCondInf: N<1!');
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I+1;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I-1;
        end;
        V := 0;
        J:=J1;
        while J<=J2 do
        begin
            V := V+AbsReal(A[I,J]);
            Inc(J);
        end;
        if IsUnit then
        begin
            V := V+1;
        end
        else
        begin
            V := V+AbsReal(A[I,I]);
        end;
        Nrm := Max(Nrm, V);
        Inc(I);
    end;
    RMatrixRCondTRInternal(A, N, IsUpper, IsUnit, False, Nrm, V);
    Result := V;
end;


(*************************************************************************
Condition number estimate of a Hermitian positive definite matrix.

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm of condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   Hermitian positive definite matrix which is given by its
                upper or lower triangle depending on the value of
                IsUpper. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A)), if matrix A is positive definite,
   -1, if matrix A is not positive definite, and its condition number
    could not be found by this algorithm.

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function HPDMatrixRCond(A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    V : Double;
    Nrm : Double;
    T : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    SetLength(T, N);
    I:=0;
    while I<=N-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        J:=J1;
        while J<=J2 do
        begin
            if I=J then
            begin
                T[I] := T[I]+AbsComplex(A[I,I]);
            end
            else
            begin
                T[I] := T[I]+AbsComplex(A[I,J]);
                T[J] := T[J]+AbsComplex(A[I,J]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        Nrm := Max(Nrm, T[I]);
        Inc(I);
    end;
    if HPDMatrixCholesky(A, N, IsUpper) then
    begin
        HPDMatrixRCondCholeskyInternal(A, N, IsUpper, True, Nrm, V);
        Result := V;
    end
    else
    begin
        Result := -1;
    end;
end;


(*************************************************************************
Estimate of a matrix condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function CMatrixRCond1(A : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
    T : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert(N>=1, 'CMatrixRCond1: N<1!');
    SetLength(T, N);
    I:=0;
    while I<=N-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[J] := T[J]+AbsComplex(A[I,J]);
            Inc(J);
        end;
        Inc(I);
    end;
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        Nrm := Max(Nrm, T[I]);
        Inc(I);
    end;
    CMatrixLU(A, N, N, Pivots);
    CMatrixRCondLUInternal(A, N, True, True, Nrm, V);
    Result := V;
end;


(*************************************************************************
Estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function CMatrixRCondInf(A : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert(N>=1, 'CMatrixRCondInf: N<1!');
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := 0;
        J:=0;
        while J<=N-1 do
        begin
            V := V+AbsComplex(A[I,J]);
            Inc(J);
        end;
        Nrm := Max(Nrm, V);
        Inc(I);
    end;
    CMatrixLU(A, N, N, Pivots);
    CMatrixRCondLUInternal(A, N, False, True, Nrm, V);
    Result := V;
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA         -   LU decomposition of a matrix in compact form. Output of
                    the RMatrixLU subroutine.
    N           -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function RMatrixLURCond1(const LUA : TReal2DArray; N : AlglibInteger):Double;
var
    V : Double;
begin
    RMatrixRCondLUInternal(LUA, N, True, False, 0, V);
    Result := V;
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition
(infinity norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA     -   LU decomposition of a matrix in compact form. Output of
                the RMatrixLU subroutine.
    N       -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function RMatrixLURCondInf(const LUA : TReal2DArray; N : AlglibInteger):Double;
var
    V : Double;
begin
    RMatrixRCondLUInternal(LUA, N, False, False, 0, V);
    Result := V;
end;


(*************************************************************************
Condition number estimate of a symmetric positive definite matrix given by
Cholesky decomposition.

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    CD  - Cholesky decomposition of matrix A,
          output of SMatrixCholesky subroutine.
    N   - size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function SPDMatrixCholeskyRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    V : Double;
begin
    SPDMatrixRCondCholeskyInternal(A, N, IsUpper, False, 0, V);
    Result := V;
end;


(*************************************************************************
Condition number estimate of a Hermitian positive definite matrix given by
Cholesky decomposition.

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    CD  - Cholesky decomposition of matrix A,
          output of SMatrixCholesky subroutine.
    N   - size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function HPDMatrixCholeskyRCond(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    V : Double;
begin
    HPDMatrixRCondCholeskyInternal(A, N, IsUpper, False, 0, V);
    Result := V;
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA         -   LU decomposition of a matrix in compact form. Output of
                    the CMatrixLU subroutine.
    N           -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function CMatrixLURCond1(const LUA : TComplex2DArray;
     N : AlglibInteger):Double;
var
    V : Double;
begin
    Assert(N>=1, 'CMatrixLURCond1: N<1!');
    CMatrixRCondLUInternal(LUA, N, True, False, Double(0.0), V);
    Result := V;
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition
(infinity norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA     -   LU decomposition of a matrix in compact form. Output of
                the CMatrixLU subroutine.
    N       -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function CMatrixLURCondInf(const LUA : TComplex2DArray;
     N : AlglibInteger):Double;
var
    V : Double;
begin
    Assert(N>=1, 'CMatrixLURCondInf: N<1!');
    CMatrixRCondLUInternal(LUA, N, False, False, Double(0.0), V);
    Result := V;
end;


(*************************************************************************
Triangular matrix: estimate of a condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A       -   matrix. Array[0..N-1, 0..N-1].
    N       -   size of A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function CMatrixTRRCond1(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
    T : TReal1DArray;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    Assert(N>=1, 'RMatrixTRRCond1: N<1!');
    SetLength(T, N);
    I:=0;
    while I<=N-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I+1;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I-1;
        end;
        J:=J1;
        while J<=J2 do
        begin
            T[J] := T[J]+AbsComplex(A[I,J]);
            Inc(J);
        end;
        if IsUnit then
        begin
            T[I] := T[I]+1;
        end
        else
        begin
            T[I] := T[I]+AbsComplex(A[I,I]);
        end;
        Inc(I);
    end;
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        Nrm := Max(Nrm, T[I]);
        Inc(I);
    end;
    CMatrixRCondTRInternal(A, N, IsUpper, IsUnit, True, Nrm, V);
    Result := V;
end;


(*************************************************************************
Triangular matrix: estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************)
function CMatrixTRRCondInf(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    Assert(N>=1, 'RMatrixTRRCondInf: N<1!');
    Nrm := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I+1;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I-1;
        end;
        V := 0;
        J:=J1;
        while J<=J2 do
        begin
            V := V+AbsComplex(A[I,J]);
            Inc(J);
        end;
        if IsUnit then
        begin
            V := V+1;
        end
        else
        begin
            V := V+AbsComplex(A[I,I]);
        end;
        Nrm := Max(Nrm, V);
        Inc(I);
    end;
    CMatrixRCondTRInternal(A, N, IsUpper, IsUnit, False, Nrm, V);
    Result := V;
end;


(*************************************************************************
Threshold for rcond: matrices with condition number beyond this  threshold
are considered singular.

Threshold must be far enough from underflow, at least Sqr(Threshold)  must
be greater than underflow.
*************************************************************************)
function RCondThreshold():Double;
begin
    Result := Sqrt(Sqrt(MinRealNumber));
end;


(*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
procedure RMatrixRCondTRInternal(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OneNorm : Boolean;
     ANORM : Double;
     var RC : Double);
var
    EX : TReal1DArray;
    EV : TReal1DArray;
    IWORK : TInteger1DArray;
    Tmp : TReal1DArray;
    V : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    KASE : AlglibInteger;
    KASE1 : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    AINVNM : Double;
    MaxGrowth : Double;
    S : Double;
    MUpper : Boolean;
    MTrans : Boolean;
    Munit : Boolean;
begin
    
    //
    // RC=0 if something happens
    //
    RC := 0;
    
    //
    // init
    //
    if OneNorm then
    begin
        KASE1 := 1;
    end
    else
    begin
        KASE1 := 2;
    end;
    MUpper := True;
    MTrans := True;
    Munit := True;
    SetLength(IWORK, N+1);
    SetLength(Tmp, N);
    
    //
    // prepare parameters for triangular solver
    //
    MaxGrowth := 1/RCondThreshold;
    S := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I+1;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I-1;
        end;
        J:=J1;
        while J<=J2 do
        begin
            S := Max(S, AbsReal(A[I,J]));
            Inc(J);
        end;
        if IsUnit then
        begin
            S := Max(S, 1);
        end
        else
        begin
            S := Max(S, AbsReal(A[I,I]));
        end;
        Inc(I);
    end;
    if AP_FP_Eq(S,0) then
    begin
        S := 1;
    end;
    S := 1/S;
    
    //
    // Scale according to S
    //
    ANORM := ANORM*S;
    
    //
    // Quick return if possible
    // We assume that ANORM<>0 after this block
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    if N=1 then
    begin
        RC := 1;
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    AINVNM := 0;
    KASE := 0;
    while True do
    begin
        RMatrixEstimateNorm(N, EV, EX, IWORK, AINVNM, KASE);
        if KASE=0 then
        begin
            Break;
        end;
        
        //
        // from 1-based array to 0-based
        //
        I:=0;
        while I<=N-1 do
        begin
            EX[I] := EX[I+1];
            Inc(I);
        end;
        
        //
        // multiply by inv(A) or inv(A')
        //
        if KASE=KASE1 then
        begin
            
            //
            // multiply by inv(A)
            //
            if  not RMatrixScaledTRSafeSolve(A, S, N, EX, IsUpper, 0, IsUnit, MaxGrowth) then
            begin
                Exit;
            end;
        end
        else
        begin
            
            //
            // multiply by inv(A')
            //
            if  not RMatrixScaledTRSafeSolve(A, S, N, EX, IsUpper, 1, IsUnit, MaxGrowth) then
            begin
                Exit;
            end;
        end;
        
        //
        // from 0-based array to 1-based
        //
        I:=N-1;
        while I>=0 do
        begin
            EX[I+1] := EX[I];
            Dec(I);
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RC := 1/AINVNM;
        RC := RC/ANORM;
        if AP_FP_Less(RC,RCondThreshold) then
        begin
            RC := 0;
        end;
    end;
end;


(*************************************************************************
Condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************)
procedure CMatrixRCondTRInternal(const A : TComplex2DArray;
     const N : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OneNorm : Boolean;
     ANORM : Double;
     var RC : Double);
var
    EX : TComplex1DArray;
    CWORK2 : TComplex1DArray;
    CWORK3 : TComplex1DArray;
    CWORK4 : TComplex1DArray;
    ISAVE : TInteger1DArray;
    RSAVE : TReal1DArray;
    KASE : AlglibInteger;
    KASE1 : AlglibInteger;
    AINVNM : Double;
    V : Complex;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    S : Double;
    MaxGrowth : Double;
begin
    
    //
    // RC=0 if something happens
    //
    RC := 0;
    
    //
    // init
    //
    if N<=0 then
    begin
        Exit;
    end;
    if N=0 then
    begin
        RC := 1;
        Exit;
    end;
    SetLength(CWORK2, N+1);
    
    //
    // prepare parameters for triangular solver
    //
    MaxGrowth := 1/RCondThreshold;
    S := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I+1;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I-1;
        end;
        J:=J1;
        while J<=J2 do
        begin
            S := Max(S, AbsComplex(A[I,J]));
            Inc(J);
        end;
        if IsUnit then
        begin
            S := Max(S, 1);
        end
        else
        begin
            S := Max(S, AbsComplex(A[I,I]));
        end;
        Inc(I);
    end;
    if AP_FP_Eq(S,0) then
    begin
        S := 1;
    end;
    S := 1/S;
    
    //
    // Scale according to S
    //
    ANORM := ANORM*S;
    
    //
    // Quick return if possible
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    AINVNM := 0;
    if OneNorm then
    begin
        KASE1 := 1;
    end
    else
    begin
        KASE1 := 2;
    end;
    KASE := 0;
    while True do
    begin
        CMatrixEstimateNorm(N, CWORK4, EX, AINVNM, KASE, ISAVE, RSAVE);
        if KASE=0 then
        begin
            Break;
        end;
        
        //
        // From 1-based to 0-based
        //
        I:=0;
        while I<=N-1 do
        begin
            EX[I] := EX[I+1];
            Inc(I);
        end;
        
        //
        // multiply by inv(A) or inv(A')
        //
        if KASE=KASE1 then
        begin
            
            //
            // multiply by inv(A)
            //
            if  not CMatrixScaledTRSafeSolve(A, S, N, EX, IsUpper, 0, IsUnit, MaxGrowth) then
            begin
                Exit;
            end;
        end
        else
        begin
            
            //
            // multiply by inv(A')
            //
            if  not CMatrixScaledTRSafeSolve(A, S, N, EX, IsUpper, 2, IsUnit, MaxGrowth) then
            begin
                Exit;
            end;
        end;
        
        //
        // from 0-based to 1-based
        //
        I:=N-1;
        while I>=0 do
        begin
            EX[I+1] := EX[I];
            Dec(I);
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RC := 1/AINVNM;
        RC := RC/ANORM;
        if AP_FP_Less(RC,RCondThreshold) then
        begin
            RC := 0;
        end;
    end;
end;


(*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
procedure SPDMatrixRCondCholeskyInternal(const CHA : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RC : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    KASE : AlglibInteger;
    AINVNM : Double;
    EX : TReal1DArray;
    EV : TReal1DArray;
    Tmp : TReal1DArray;
    IWORK : TInteger1DArray;
    SA : Double;
    V : Double;
    MaxGrowth : Double;
begin
    Assert(N>=1);
    SetLength(Tmp, N);
    
    //
    // RC=0 if something happens
    //
    RC := 0;
    
    //
    // prepare parameters for triangular solver
    //
    MaxGrowth := 1/RCondThreshold;
    SA := 0;
    if IsUpper then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=I;
            while J<=N-1 do
            begin
                SA := Max(SA, AbsComplex(C_Complex(CHA[I,J])));
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=I do
            begin
                SA := Max(SA, AbsComplex(C_Complex(CHA[I,J])));
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if AP_FP_Eq(SA,0) then
    begin
        SA := 1;
    end;
    SA := 1/SA;
    
    //
    // Estimate the norm of A.
    //
    if  not IsNormProvided then
    begin
        KASE := 0;
        ANORM := 0;
        while True do
        begin
            RMatrixEstimateNorm(N, EV, EX, IWORK, ANORM, KASE);
            if KASE=0 then
            begin
                Break;
            end;
            if IsUpper then
            begin
                
                //
                // Multiply by U
                //
                I:=1;
                while I<=N do
                begin
                    V := APVDotProduct(@CHA[I-1][0], I-1, N-1, @EX[0], I, N);
                    EX[I] := V;
                    Inc(I);
                end;
                APVMul(@EX[0], 1, N, SA);
                
                //
                // Multiply by U'
                //
                I:=0;
                while I<=N-1 do
                begin
                    Tmp[I] := 0;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := EX[I+1];
                    APVAdd(@Tmp[0], I, N-1, @CHA[I][0], I, N-1, V);
                    Inc(I);
                end;
                APVMove(@EX[0], 1, N, @Tmp[0], 0, N-1);
                APVMul(@EX[0], 1, N, SA);
            end
            else
            begin
                
                //
                // Multiply by L'
                //
                I:=0;
                while I<=N-1 do
                begin
                    Tmp[I] := 0;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := EX[I+1];
                    APVAdd(@Tmp[0], 0, I, @CHA[I][0], 0, I, V);
                    Inc(I);
                end;
                APVMove(@EX[0], 1, N, @Tmp[0], 0, N-1);
                APVMul(@EX[0], 1, N, SA);
                
                //
                // Multiply by L
                //
                I:=N;
                while I>=1 do
                begin
                    V := APVDotProduct(@CHA[I-1][0], 0, I-1, @EX[0], 1, I);
                    EX[I] := V;
                    Dec(I);
                end;
                APVMul(@EX[0], 1, N, SA);
            end;
        end;
    end;
    
    //
    // Quick return if possible
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    if N=1 then
    begin
        RC := 1;
        Exit;
    end;
    
    //
    // Estimate the 1-norm of inv(A).
    //
    KASE := 0;
    while True do
    begin
        RMatrixEstimateNorm(N, EV, EX, IWORK, AINVNM, KASE);
        if KASE=0 then
        begin
            Break;
        end;
        I:=0;
        while I<=N-1 do
        begin
            EX[I] := EX[I+1];
            Inc(I);
        end;
        if IsUpper then
        begin
            
            //
            // Multiply by inv(U').
            //
            if  not RMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 1, False, MaxGrowth) then
            begin
                Exit;
            end;
            
            //
            // Multiply by inv(U).
            //
            if  not RMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 0, False, MaxGrowth) then
            begin
                Exit;
            end;
        end
        else
        begin
            
            //
            // Multiply by inv(L).
            //
            if  not RMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 0, False, MaxGrowth) then
            begin
                Exit;
            end;
            
            //
            // Multiply by inv(L').
            //
            if  not RMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 1, False, MaxGrowth) then
            begin
                Exit;
            end;
        end;
        I:=N-1;
        while I>=0 do
        begin
            EX[I+1] := EX[I];
            Dec(I);
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        V := 1/AINVNM;
        RC := V/ANORM;
        if AP_FP_Less(RC,RCondThreshold) then
        begin
            RC := 0;
        end;
    end;
end;


(*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
procedure HPDMatrixRCondCholeskyInternal(const CHA : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RC : Double);
var
    ISAVE : TInteger1DArray;
    RSAVE : TReal1DArray;
    EX : TComplex1DArray;
    EV : TComplex1DArray;
    Tmp : TComplex1DArray;
    KASE : AlglibInteger;
    AINVNM : Double;
    V : Complex;
    I : AlglibInteger;
    J : AlglibInteger;
    SA : Double;
    MaxGrowth : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=1);
    SetLength(Tmp, N);
    
    //
    // RC=0 if something happens
    //
    RC := 0;
    
    //
    // prepare parameters for triangular solver
    //
    MaxGrowth := 1/RCondThreshold;
    SA := 0;
    if IsUpper then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=I;
            while J<=N-1 do
            begin
                SA := Max(SA, AbsComplex(CHA[I,J]));
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=I do
            begin
                SA := Max(SA, AbsComplex(CHA[I,J]));
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if AP_FP_Eq(SA,0) then
    begin
        SA := 1;
    end;
    SA := 1/SA;
    
    //
    // Estimate the norm of A
    //
    if  not IsNormProvided then
    begin
        ANORM := 0;
        KASE := 0;
        while True do
        begin
            CMatrixEstimateNorm(N, EV, EX, ANORM, KASE, ISAVE, RSAVE);
            if KASE=0 then
            begin
                Break;
            end;
            if IsUpper then
            begin
                
                //
                // Multiply by U
                //
                I:=1;
                while I<=N do
                begin
                    i1_ := (I)-(I-1);
                    V := C_Complex(0.0);
                    for i_ := I-1 to N-1 do
                    begin
                        V := C_Add(V,C_Mul(CHA[I-1,i_],EX[i_+i1_]));
                    end;
                    EX[I] := V;
                    Inc(I);
                end;
                for i_ := 1 to N do
                begin
                    EX[i_] := C_MulR(EX[i_],SA);
                end;
                
                //
                // Multiply by U'
                //
                I:=0;
                while I<=N-1 do
                begin
                    Tmp[I] := C_Complex(0);
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := EX[I+1];
                    for i_ := I to N-1 do
                    begin
                        Tmp[i_] := C_Add(Tmp[i_], C_Mul(V, Conj(CHA[I,i_])));
                    end;
                    Inc(I);
                end;
                i1_ := (0) - (1);
                for i_ := 1 to N do
                begin
                    EX[i_] := Tmp[i_+i1_];
                end;
                for i_ := 1 to N do
                begin
                    EX[i_] := C_MulR(EX[i_],SA);
                end;
            end
            else
            begin
                
                //
                // Multiply by L'
                //
                I:=0;
                while I<=N-1 do
                begin
                    Tmp[I] := C_Complex(0);
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := EX[I+1];
                    for i_ := 0 to I do
                    begin
                        Tmp[i_] := C_Add(Tmp[i_], C_Mul(V, Conj(CHA[I,i_])));
                    end;
                    Inc(I);
                end;
                i1_ := (0) - (1);
                for i_ := 1 to N do
                begin
                    EX[i_] := Tmp[i_+i1_];
                end;
                for i_ := 1 to N do
                begin
                    EX[i_] := C_MulR(EX[i_],SA);
                end;
                
                //
                // Multiply by L
                //
                I:=N;
                while I>=1 do
                begin
                    i1_ := (1)-(0);
                    V := C_Complex(0.0);
                    for i_ := 0 to I-1 do
                    begin
                        V := C_Add(V,C_Mul(CHA[I-1,i_],EX[i_+i1_]));
                    end;
                    EX[I] := V;
                    Dec(I);
                end;
                for i_ := 1 to N do
                begin
                    EX[i_] := C_MulR(EX[i_],SA);
                end;
            end;
        end;
    end;
    
    //
    // Quick return if possible
    // After this block we assume that ANORM<>0
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    if N=1 then
    begin
        RC := 1;
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    AINVNM := 0;
    KASE := 0;
    while True do
    begin
        CMatrixEstimateNorm(N, EV, EX, AINVNM, KASE, ISAVE, RSAVE);
        if KASE=0 then
        begin
            Break;
        end;
        I:=0;
        while I<=N-1 do
        begin
            EX[I] := EX[I+1];
            Inc(I);
        end;
        if IsUpper then
        begin
            
            //
            // Multiply by inv(U').
            //
            if  not CMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 2, False, MaxGrowth) then
            begin
                Exit;
            end;
            
            //
            // Multiply by inv(U).
            //
            if  not CMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 0, False, MaxGrowth) then
            begin
                Exit;
            end;
        end
        else
        begin
            
            //
            // Multiply by inv(L).
            //
            if  not CMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 0, False, MaxGrowth) then
            begin
                Exit;
            end;
            
            //
            // Multiply by inv(L').
            //
            if  not CMatrixScaledTRSafeSolve(CHA, SA, N, EX, IsUpper, 2, False, MaxGrowth) then
            begin
                Exit;
            end;
        end;
        I:=N-1;
        while I>=0 do
        begin
            EX[I+1] := EX[I];
            Dec(I);
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RC := 1/AINVNM;
        RC := RC/ANORM;
        if AP_FP_Less(RC,RCondThreshold) then
        begin
            RC := 0;
        end;
    end;
end;


(*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
procedure RMatrixRCondLUInternal(const LUA : TReal2DArray;
     N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RC : Double);
var
    EX : TReal1DArray;
    EV : TReal1DArray;
    IWORK : TInteger1DArray;
    Tmp : TReal1DArray;
    V : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    KASE : AlglibInteger;
    KASE1 : AlglibInteger;
    AINVNM : Double;
    MaxGrowth : Double;
    SU : Double;
    SL : Double;
    MUpper : Boolean;
    MTrans : Boolean;
    Munit : Boolean;
begin
    
    //
    // RC=0 if something happens
    //
    RC := 0;
    
    //
    // init
    //
    if OneNorm then
    begin
        KASE1 := 1;
    end
    else
    begin
        KASE1 := 2;
    end;
    MUpper := True;
    MTrans := True;
    Munit := True;
    SetLength(IWORK, N+1);
    SetLength(Tmp, N);
    
    //
    // prepare parameters for triangular solver
    //
    MaxGrowth := 1/RCondThreshold;
    SU := 0;
    SL := 1;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I-1 do
        begin
            SL := Max(SL, AbsReal(LUA[I,J]));
            Inc(J);
        end;
        J:=I;
        while J<=N-1 do
        begin
            SU := Max(SU, AbsReal(LUA[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(SU,0) then
    begin
        SU := 1;
    end;
    SU := 1/SU;
    SL := 1/SL;
    
    //
    // Estimate the norm of A.
    //
    if  not IsANormProvided then
    begin
        KASE := 0;
        ANORM := 0;
        while True do
        begin
            RMatrixEstimateNorm(N, EV, EX, IWORK, ANORM, KASE);
            if KASE=0 then
            begin
                Break;
            end;
            if KASE=KASE1 then
            begin
                
                //
                // Multiply by U
                //
                I:=1;
                while I<=N do
                begin
                    V := APVDotProduct(@LUA[I-1][0], I-1, N-1, @EX[0], I, N);
                    EX[I] := V;
                    Inc(I);
                end;
                
                //
                // Multiply by L
                //
                I:=N;
                while I>=1 do
                begin
                    if I>1 then
                    begin
                        V := APVDotProduct(@LUA[I-1][0], 0, I-2, @EX[0], 1, I-1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    EX[I] := EX[I]+V;
                    Dec(I);
                end;
            end
            else
            begin
                
                //
                // Multiply by L'
                //
                I:=0;
                while I<=N-1 do
                begin
                    Tmp[I] := 0;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := EX[I+1];
                    if I>=1 then
                    begin
                        APVAdd(@Tmp[0], 0, I-1, @LUA[I][0], 0, I-1, V);
                    end;
                    Tmp[I] := Tmp[I]+V;
                    Inc(I);
                end;
                APVMove(@EX[0], 1, N, @Tmp[0], 0, N-1);
                
                //
                // Multiply by U'
                //
                I:=0;
                while I<=N-1 do
                begin
                    Tmp[I] := 0;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := EX[I+1];
                    APVAdd(@Tmp[0], I, N-1, @LUA[I][0], I, N-1, V);
                    Inc(I);
                end;
                APVMove(@EX[0], 1, N, @Tmp[0], 0, N-1);
            end;
        end;
    end;
    
    //
    // Scale according to SU/SL
    //
    ANORM := ANORM*SU*SL;
    
    //
    // Quick return if possible
    // We assume that ANORM<>0 after this block
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    if N=1 then
    begin
        RC := 1;
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    AINVNM := 0;
    KASE := 0;
    while True do
    begin
        RMatrixEstimateNorm(N, EV, EX, IWORK, AINVNM, KASE);
        if KASE=0 then
        begin
            Break;
        end;
        
        //
        // from 1-based array to 0-based
        //
        I:=0;
        while I<=N-1 do
        begin
            EX[I] := EX[I+1];
            Inc(I);
        end;
        
        //
        // multiply by inv(A) or inv(A')
        //
        if KASE=KASE1 then
        begin
            
            //
            // Multiply by inv(L).
            //
            if  not RMatrixScaledTRSafeSolve(LUA, SL, N, EX,  not MUpper, 0, Munit, MaxGrowth) then
            begin
                Exit;
            end;
            
            //
            // Multiply by inv(U).
            //
            if  not RMatrixScaledTRSafeSolve(LUA, SU, N, EX, MUpper, 0,  not Munit, MaxGrowth) then
            begin
                Exit;
            end;
        end
        else
        begin
            
            //
            // Multiply by inv(U').
            //
            if  not RMatrixScaledTRSafeSolve(LUA, SU, N, EX, MUpper, 1,  not Munit, MaxGrowth) then
            begin
                Exit;
            end;
            
            //
            // Multiply by inv(L').
            //
            if  not RMatrixScaledTRSafeSolve(LUA, SL, N, EX,  not MUpper, 1, Munit, MaxGrowth) then
            begin
                Exit;
            end;
        end;
        
        //
        // from 0-based array to 1-based
        //
        I:=N-1;
        while I>=0 do
        begin
            EX[I+1] := EX[I];
            Dec(I);
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RC := 1/AINVNM;
        RC := RC/ANORM;
        if AP_FP_Less(RC,RCondThreshold) then
        begin
            RC := 0;
        end;
    end;
end;


(*************************************************************************
Condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************)
procedure CMatrixRCondLUInternal(const LUA : TComplex2DArray;
     const N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RC : Double);
var
    EX : TComplex1DArray;
    CWORK2 : TComplex1DArray;
    CWORK3 : TComplex1DArray;
    CWORK4 : TComplex1DArray;
    ISAVE : TInteger1DArray;
    RSAVE : TReal1DArray;
    KASE : AlglibInteger;
    KASE1 : AlglibInteger;
    AINVNM : Double;
    V : Complex;
    I : AlglibInteger;
    J : AlglibInteger;
    SU : Double;
    SL : Double;
    MaxGrowth : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(CWORK2, N+1);
    RC := 0;
    if N=0 then
    begin
        RC := 1;
        Exit;
    end;
    
    //
    // prepare parameters for triangular solver
    //
    MaxGrowth := 1/RCondThreshold;
    SU := 0;
    SL := 1;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I-1 do
        begin
            SL := Max(SL, AbsComplex(LUA[I,J]));
            Inc(J);
        end;
        J:=I;
        while J<=N-1 do
        begin
            SU := Max(SU, AbsComplex(LUA[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(SU,0) then
    begin
        SU := 1;
    end;
    SU := 1/SU;
    SL := 1/SL;
    
    //
    // Estimate the norm of SU*SL*A.
    //
    if  not IsANormProvided then
    begin
        ANORM := 0;
        if OneNorm then
        begin
            KASE1 := 1;
        end
        else
        begin
            KASE1 := 2;
        end;
        KASE := 0;
        repeat
            CMatrixEstimateNorm(N, CWORK4, EX, ANORM, KASE, ISAVE, RSAVE);
            if KASE<>0 then
            begin
                if KASE=KASE1 then
                begin
                    
                    //
                    // Multiply by U
                    //
                    I:=1;
                    while I<=N do
                    begin
                        i1_ := (I)-(I-1);
                        V := C_Complex(0.0);
                        for i_ := I-1 to N-1 do
                        begin
                            V := C_Add(V,C_Mul(LUA[I-1,i_],EX[i_+i1_]));
                        end;
                        EX[I] := V;
                        Inc(I);
                    end;
                    
                    //
                    // Multiply by L
                    //
                    I:=N;
                    while I>=1 do
                    begin
                        V := C_Complex(0);
                        if I>1 then
                        begin
                            i1_ := (1)-(0);
                            V := C_Complex(0.0);
                            for i_ := 0 to I-2 do
                            begin
                                V := C_Add(V,C_Mul(LUA[I-1,i_],EX[i_+i1_]));
                            end;
                        end;
                        EX[I] := C_Add(V,EX[I]);
                        Dec(I);
                    end;
                end
                else
                begin
                    
                    //
                    // Multiply by L'
                    //
                    I:=1;
                    while I<=N do
                    begin
                        CWORK2[I] := C_Complex(0);
                        Inc(I);
                    end;
                    I:=1;
                    while I<=N do
                    begin
                        V := EX[I];
                        if I>1 then
                        begin
                            i1_ := (0) - (1);
                            for i_ := 1 to I-1 do
                            begin
                                CWORK2[i_] := C_Add(CWORK2[i_], C_Mul(V, Conj(LUA[I-1,i_+i1_])));
                            end;
                        end;
                        CWORK2[I] := C_Add(CWORK2[I],V);
                        Inc(I);
                    end;
                    
                    //
                    // Multiply by U'
                    //
                    I:=1;
                    while I<=N do
                    begin
                        EX[I] := C_Complex(0);
                        Inc(I);
                    end;
                    I:=1;
                    while I<=N do
                    begin
                        V := CWORK2[I];
                        i1_ := (I-1) - (I);
                        for i_ := I to N do
                        begin
                            EX[i_] := C_Add(EX[i_], C_Mul(V, Conj(LUA[I-1,i_+i1_])));
                        end;
                        Inc(I);
                    end;
                end;
            end;
        until KASE=0;
    end;
    
    //
    // Scale according to SU/SL
    //
    ANORM := ANORM*SU*SL;
    
    //
    // Quick return if possible
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    AINVNM := 0;
    if OneNorm then
    begin
        KASE1 := 1;
    end
    else
    begin
        KASE1 := 2;
    end;
    KASE := 0;
    while True do
    begin
        CMatrixEstimateNorm(N, CWORK4, EX, AINVNM, KASE, ISAVE, RSAVE);
        if KASE=0 then
        begin
            Break;
        end;
        
        //
        // From 1-based to 0-based
        //
        I:=0;
        while I<=N-1 do
        begin
            EX[I] := EX[I+1];
            Inc(I);
        end;
        
        //
        // multiply by inv(A) or inv(A')
        //
        if KASE=KASE1 then
        begin
            
            //
            // Multiply by inv(L).
            //
            if  not CMatrixScaledTRSafeSolve(LUA, SL, N, EX, False, 0, True, MaxGrowth) then
            begin
                RC := 0;
                Exit;
            end;
            
            //
            // Multiply by inv(U).
            //
            if  not CMatrixScaledTRSafeSolve(LUA, SU, N, EX, True, 0, False, MaxGrowth) then
            begin
                RC := 0;
                Exit;
            end;
        end
        else
        begin
            
            //
            // Multiply by inv(U').
            //
            if  not CMatrixScaledTRSafeSolve(LUA, SU, N, EX, True, 2, False, MaxGrowth) then
            begin
                RC := 0;
                Exit;
            end;
            
            //
            // Multiply by inv(L').
            //
            if  not CMatrixScaledTRSafeSolve(LUA, SL, N, EX, False, 2, True, MaxGrowth) then
            begin
                RC := 0;
                Exit;
            end;
        end;
        
        //
        // from 0-based to 1-based
        //
        I:=N-1;
        while I>=0 do
        begin
            EX[I+1] := EX[I];
            Dec(I);
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RC := 1/AINVNM;
        RC := RC/ANORM;
        if AP_FP_Less(RC,RCondThreshold) then
        begin
            RC := 0;
        end;
    end;
end;


(*************************************************************************
Internal subroutine for matrix norm estimation

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
procedure RMatrixEstimateNorm(N : AlglibInteger;
     var V : TReal1DArray;
     var X : TReal1DArray;
     var ISGN : TInteger1DArray;
     var EST : Double;
     var KASE : AlglibInteger);
var
    ITMAX : AlglibInteger;
    I : AlglibInteger;
    T : Double;
    Flg : Boolean;
    PosITER : AlglibInteger;
    PosJ : AlglibInteger;
    PosJLAST : AlglibInteger;
    PosJUMP : AlglibInteger;
    PosALTSGN : AlglibInteger;
    PosESTOLD : AlglibInteger;
    PosTEMP : AlglibInteger;
begin
    ITMAX := 5;
    PosALTSGN := N+1;
    PosESTOLD := N+2;
    PosTEMP := N+3;
    PosITER := N+1;
    PosJ := N+2;
    PosJLAST := N+3;
    PosJUMP := N+4;
    if KASE=0 then
    begin
        SetLength(V, N+4);
        SetLength(X, N+1);
        SetLength(ISGN, N+5);
        T := AP_Double(1)/N;
        I:=1;
        while I<=N do
        begin
            X[I] := T;
            Inc(I);
        end;
        KASE := 1;
        ISGN[PosJUMP] := 1;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 1)
    //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    //
    if ISGN[PosJUMP]=1 then
    begin
        if N=1 then
        begin
            V[1] := X[1];
            EST := ABSReal(V[1]);
            KASE := 0;
            Exit;
        end;
        EST := 0;
        I:=1;
        while I<=N do
        begin
            EST := EST+AbsReal(X[I]);
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            if AP_FP_Greater_Eq(X[I],0) then
            begin
                X[I] := 1;
            end
            else
            begin
                X[I] := -1;
            end;
            ISGN[I] := Sign(X[I]);
            Inc(I);
        end;
        KASE := 2;
        ISGN[PosJUMP] := 2;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 2)
    //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
    //
    if ISGN[PosJUMP]=2 then
    begin
        ISGN[PosJ] := 1;
        I:=2;
        while I<=N do
        begin
            if AP_FP_Greater(AbsReal(X[I]),AbsReal(X[ISGN[PosJ]])) then
            begin
                ISGN[PosJ] := I;
            end;
            Inc(I);
        end;
        ISGN[PosITER] := 2;
        
        //
        // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
        //
        I:=1;
        while I<=N do
        begin
            X[I] := 0;
            Inc(I);
        end;
        X[ISGN[PosJ]] := 1;
        KASE := 1;
        ISGN[PosJUMP] := 3;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 3)
    //     X HAS BEEN OVERWRITTEN BY A*X.
    //
    if ISGN[PosJUMP]=3 then
    begin
        APVMove(@V[0], 1, N, @X[0], 1, N);
        V[PosESTOLD] := EST;
        EST := 0;
        I:=1;
        while I<=N do
        begin
            EST := EST+AbsReal(V[I]);
            Inc(I);
        end;
        Flg := False;
        I:=1;
        while I<=N do
        begin
            if AP_FP_Greater_Eq(X[I],0) and (ISGN[I]<0) or AP_FP_Less(X[I],0) and (ISGN[I]>=0) then
            begin
                Flg := True;
            end;
            Inc(I);
        end;
        
        //
        // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
        // OR MAY BE CYCLING.
        //
        if  not Flg or AP_FP_Less_Eq(EST,V[PosESTOLD]) then
        begin
            V[PosALTSGN] := 1;
            I:=1;
            while I<=N do
            begin
                X[I] := V[PosALTSGN]*(1+AP_Double((I-1))/(N-1));
                V[PosALTSGN] := -V[PosALTSGN];
                Inc(I);
            end;
            KASE := 1;
            ISGN[PosJUMP] := 5;
            Exit;
        end;
        I:=1;
        while I<=N do
        begin
            if AP_FP_Greater_Eq(X[I],0) then
            begin
                X[I] := 1;
                ISGN[I] := 1;
            end
            else
            begin
                X[I] := -1;
                ISGN[I] := -1;
            end;
            Inc(I);
        end;
        KASE := 2;
        ISGN[PosJUMP] := 4;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 4)
    //     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
    //
    if ISGN[PosJUMP]=4 then
    begin
        ISGN[PosJLAST] := ISGN[PosJ];
        ISGN[PosJ] := 1;
        I:=2;
        while I<=N do
        begin
            if AP_FP_Greater(AbsReal(X[I]),AbsReal(X[ISGN[PosJ]])) then
            begin
                ISGN[PosJ] := I;
            end;
            Inc(I);
        end;
        if AP_FP_Neq(X[ISGN[PosJLAST]],ABSReal(X[ISGN[PosJ]])) and (ISGN[PosITER]<ITMAX) then
        begin
            ISGN[PosITER] := ISGN[PosITER]+1;
            I:=1;
            while I<=N do
            begin
                X[I] := 0;
                Inc(I);
            end;
            X[ISGN[PosJ]] := 1;
            KASE := 1;
            ISGN[PosJUMP] := 3;
            Exit;
        end;
        
        //
        // ITERATION COMPLETE.  FINAL STAGE.
        //
        V[PosALTSGN] := 1;
        I:=1;
        while I<=N do
        begin
            X[I] := V[PosALTSGN]*(1+AP_Double((I-1))/(N-1));
            V[PosALTSGN] := -V[PosALTSGN];
            Inc(I);
        end;
        KASE := 1;
        ISGN[PosJUMP] := 5;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 5)
    //     X HAS BEEN OVERWRITTEN BY A*X.
    //
    if ISGN[PosJUMP]=5 then
    begin
        V[PosTEMP] := 0;
        I:=1;
        while I<=N do
        begin
            V[PosTEMP] := V[PosTEMP]+AbsReal(X[I]);
            Inc(I);
        end;
        V[PosTEMP] := 2*V[PosTEMP]/(3*N);
        if AP_FP_Greater(V[PosTEMP],EST) then
        begin
            APVMove(@V[0], 1, N, @X[0], 1, N);
            EST := V[PosTEMP];
        end;
        KASE := 0;
        Exit;
    end;
end;


procedure CMatrixEstimateNorm(const N : AlglibInteger;
     var V : TComplex1DArray;
     var X : TComplex1DArray;
     var EST : Double;
     var KASE : AlglibInteger;
     var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray);
var
    ITMAX : AlglibInteger;
    I : AlglibInteger;
    ITER : AlglibInteger;
    J : AlglibInteger;
    JLAST : AlglibInteger;
    JUMP : AlglibInteger;
    ABSXI : Double;
    ALTSGN : Double;
    ESTOLD : Double;
    SAFMIN : Double;
    TEMP : Double;
    i_ : AlglibInteger;
begin
    
    //
    //Executable Statements ..
    //
    ITMAX := 5;
    SAFMIN := MinRealNumber;
    if KASE=0 then
    begin
        SetLength(V, N+1);
        SetLength(X, N+1);
        SetLength(ISAVE, 5);
        SetLength(RSAVE, 4);
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(AP_Double(1)/N);
            Inc(I);
        end;
        KASE := 1;
        JUMP := 1;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    InternalComplexRCondLoadAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
    
    //
    // ENTRY   (JUMP = 1)
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    //
    if JUMP=1 then
    begin
        if N=1 then
        begin
            V[1] := X[1];
            EST := AbsComplex(V[1]);
            KASE := 0;
            InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
            Exit;
        end;
        EST := InternalComplexRCondSCSUM1(X, N);
        I:=1;
        while I<=N do
        begin
            ABSXI := AbsComplex(X[I]);
            if AP_FP_Greater(ABSXI,SAFMIN) then
            begin
                X[I] := C_DivR(X[I],ABSXI);
            end
            else
            begin
                X[I] := C_Complex(1);
            end;
            Inc(I);
        end;
        KASE := 2;
        JUMP := 2;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 2)
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
    //
    if JUMP=2 then
    begin
        J := InternalComplexRCondICMAX1(X, N);
        ITER := 2;
        
        //
        // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
        //
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(0);
            Inc(I);
        end;
        X[J] := C_Complex(1);
        KASE := 1;
        JUMP := 3;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 3)
    // X HAS BEEN OVERWRITTEN BY A*X.
    //
    if JUMP=3 then
    begin
        for i_ := 1 to N do
        begin
            V[i_] := X[i_];
        end;
        ESTOLD := EST;
        EST := InternalComplexRCondSCSUM1(V, N);
        
        //
        // TEST FOR CYCLING.
        //
        if AP_FP_Less_Eq(EST,ESTOLD) then
        begin
            
            //
            // ITERATION COMPLETE.  FINAL STAGE.
            //
            ALTSGN := 1;
            I:=1;
            while I<=N do
            begin
                X[I] := C_Complex(ALTSGN*(1+AP_Double((I-1))/(N-1)));
                ALTSGN := -ALTSGN;
                Inc(I);
            end;
            KASE := 1;
            JUMP := 5;
            InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
            Exit;
        end;
        I:=1;
        while I<=N do
        begin
            ABSXI := AbsComplex(X[I]);
            if AP_FP_Greater(ABSXI,SAFMIN) then
            begin
                X[I] := C_DivR(X[I],ABSXI);
            end
            else
            begin
                X[I] := C_Complex(1);
            end;
            Inc(I);
        end;
        KASE := 2;
        JUMP := 4;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 4)
    // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
    //
    if JUMP=4 then
    begin
        JLAST := J;
        J := InternalComplexRCondICMAX1(X, N);
        if AP_FP_Neq(AbsComplex(X[JLAST]),AbsComplex(X[J])) and (ITER<ITMAX) then
        begin
            ITER := ITER+1;
            
            //
            // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
            //
            I:=1;
            while I<=N do
            begin
                X[I] := C_Complex(0);
                Inc(I);
            end;
            X[J] := C_Complex(1);
            KASE := 1;
            JUMP := 3;
            InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
            Exit;
        end;
        
        //
        // ITERATION COMPLETE.  FINAL STAGE.
        //
        ALTSGN := 1;
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(ALTSGN*(1+AP_Double((I-1))/(N-1)));
            ALTSGN := -ALTSGN;
            Inc(I);
        end;
        KASE := 1;
        JUMP := 5;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 5)
    // X HAS BEEN OVERWRITTEN BY A*X.
    //
    if JUMP=5 then
    begin
        TEMP := 2*(InternalComplexRCondSCSUM1(X, N)/(3*N));
        if AP_FP_Greater(TEMP,EST) then
        begin
            for i_ := 1 to N do
            begin
                V[i_] := X[i_];
            end;
            EST := TEMP;
        end;
        KASE := 0;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
end;


function InternalComplexRCondSCSUM1(const X : TComplex1DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
begin
    Result := 0;
    I:=1;
    while I<=N do
    begin
        Result := Result+AbsComplex(X[I]);
        Inc(I);
    end;
end;


function InternalComplexRCondICMAX1(const X : TComplex1DArray;
     N : AlglibInteger):AlglibInteger;
var
    I : AlglibInteger;
    M : Double;
begin
    Result := 1;
    M := AbsComplex(X[1]);
    I:=2;
    while I<=N do
    begin
        if AP_FP_Greater(AbsComplex(X[I]),M) then
        begin
            Result := I;
            M := AbsComplex(X[I]);
        end;
        Inc(I);
    end;
end;


procedure InternalComplexRCondSaveAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);
begin
    ISAVE[0] := I;
    ISAVE[1] := ITER;
    ISAVE[2] := J;
    ISAVE[3] := JLAST;
    ISAVE[4] := JUMP;
    RSAVE[0] := ABSXI;
    RSAVE[1] := ALTSGN;
    RSAVE[2] := ESTOLD;
    RSAVE[3] := TEMP;
end;


procedure InternalComplexRCondLoadAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);
begin
    I := ISAVE[0];
    ITER := ISAVE[1];
    J := ISAVE[2];
    JLAST := ISAVE[3];
    JUMP := ISAVE[4];
    ABSXI := RSAVE[0];
    ALTSGN := RSAVE[1];
    ESTOLD := RSAVE[2];
    TEMP := RSAVE[3];
end;


end.