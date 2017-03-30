{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

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
unit trfac;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas;

procedure RMatrixLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure CMatrixLU(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
function HPDMatrixCholesky(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function SPDMatrixCholesky(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
procedure RMatrixLUP(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure CMatrixLUP(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure RMatrixPLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure CMatrixPLU(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);

implementation

procedure CMatrixLUPRec(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);forward;
procedure RMatrixLUPRec(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);forward;
procedure CMatrixPLURec(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);forward;
procedure RMatrixPLURec(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);forward;
procedure CMatrixLUP2(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);forward;
procedure RMatrixLUP2(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);forward;
procedure CMatrixPLU2(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);forward;
procedure RMatrixPLU2(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);forward;
function HPDMatrixCholeskyRec(var A : TComplex2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TComplex1DArray):Boolean;forward;
function SPDMatrixCholeskyRec(var A : TReal2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TReal1DArray):Boolean;forward;
function HPDMatrixCholesky2(var AAA : TComplex2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TComplex1DArray):Boolean;forward;
function SPDMatrixCholesky2(var AAA : TReal2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TReal1DArray):Boolean;forward;


(*************************************************************************
LU decomposition of a general real matrix with row pivoting

A is represented as A = P*L*U, where:
* L is lower unitriangular matrix
* U is upper triangular matrix
* P = P0*P1*...*PK, K=min(M,N)-1,
  Pi - permutation matrix for I and Pivots[I]

This is cache-oblivous implementation of LU decomposition.
It is optimized for square matrices. As for rectangular matrices:
* best case - M>>N
* worst case - N>>M, small M, large N, matrix does not fit in CPU cache

INPUT PARAMETERS:
    A       -   array[0..M-1, 0..N-1].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.


OUTPUT PARAMETERS:
    A       -   matrices L and U in compact form:
                * L is stored under main diagonal
                * U is stored on and above main diagonal
    Pivots  -   permutation matrix in compact form.
                array[0..Min(M-1,N-1)].

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
begin
    Assert(M>0, 'RMatrixLU: incorrect M!');
    Assert(N>0, 'RMatrixLU: incorrect N!');
    RMatrixPLU(A, M, N, Pivots);
end;


(*************************************************************************
LU decomposition of a general complex matrix with row pivoting

A is represented as A = P*L*U, where:
* L is lower unitriangular matrix
* U is upper triangular matrix
* P = P0*P1*...*PK, K=min(M,N)-1,
  Pi - permutation matrix for I and Pivots[I]

This is cache-oblivous implementation of LU decomposition. It is optimized
for square matrices. As for rectangular matrices:
* best case - M>>N
* worst case - N>>M, small M, large N, matrix does not fit in CPU cache

INPUT PARAMETERS:
    A       -   array[0..M-1, 0..N-1].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.


OUTPUT PARAMETERS:
    A       -   matrices L and U in compact form:
                * L is stored under main diagonal
                * U is stored on and above main diagonal
    Pivots  -   permutation matrix in compact form.
                array[0..Min(M-1,N-1)].

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixLU(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
begin
    Assert(M>0, 'CMatrixLU: incorrect M!');
    Assert(N>0, 'CMatrixLU: incorrect N!');
    CMatrixPLU(A, M, N, Pivots);
end;


(*************************************************************************
Cache-oblivious Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  Hermitian  positive-
definite matrix. The result of an algorithm is a representation  of  A  as
A=U'*U  or A=L*L' (here X' detones conj(X^T)).

INPUT PARAMETERS:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

OUTPUT PARAMETERS:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U'*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

RESULT:
    If  the  matrix  is  positive-definite,  the  function  returns  True.
    Otherwise, the function returns False. Contents of A is not determined
    in such case.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
function HPDMatrixCholesky(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    Tmp : TComplex1DArray;
begin
    if N<1 then
    begin
        Result := False;
        Exit;
    end;
    SetLength(Tmp, 2*N);
    Result := HPDMatrixCholeskyRec(A, 0, N, IsUpper, Tmp);
end;


(*************************************************************************
Cache-oblivious Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
definite matrix. The result of an algorithm is a representation  of  A  as
A=U^T*U  or A=L*L^T

INPUT PARAMETERS:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

OUTPUT PARAMETERS:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U^T*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

RESULT:
    If  the  matrix  is  positive-definite,  the  function  returns  True.
    Otherwise, the function returns False. Contents of A is not determined
    in such case.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
function SPDMatrixCholesky(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    Tmp : TReal1DArray;
begin
    if N<1 then
    begin
        Result := False;
        Exit;
    end;
    SetLength(Tmp, 2*N);
    Result := SPDMatrixCholeskyRec(A, 0, N, IsUpper, Tmp);
end;


procedure RMatrixLUP(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    Tmp : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MX : Double;
    V : Double;
begin
    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    Assert(M>0, 'RMatrixLUP: incorrect M!');
    Assert(N>0, 'RMatrixLUP: incorrect N!');
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    MX := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            MX := Max(MX, AbsReal(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Neq(MX,0) then
    begin
        V := 1/MX;
        I:=0;
        while I<=M-1 do
        begin
            APVMul(@A[I][0], 0, N-1, V);
            Inc(I);
        end;
    end;
    SetLength(Pivots, Min(M, N));
    SetLength(Tmp, 2*Max(M, N));
    RMatrixLUPRec(A, 0, M, N, Pivots, Tmp);
    if AP_FP_Neq(MX,0) then
    begin
        V := MX;
        I:=0;
        while I<=M-1 do
        begin
            APVMul(@A[I][0], 0, Min(I, N-1), V);
            Inc(I);
        end;
    end;
end;


procedure CMatrixLUP(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    Tmp : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MX : Double;
    V : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    Assert(M>0, 'CMatrixLUP: incorrect M!');
    Assert(N>0, 'CMatrixLUP: incorrect N!');
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    MX := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            MX := Max(MX, AbsComplex(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Neq(MX,0) then
    begin
        V := 1/MX;
        I:=0;
        while I<=M-1 do
        begin
            for i_ := 0 to N-1 do
            begin
                A[I,i_] := C_MulR(A[I,i_],V);
            end;
            Inc(I);
        end;
    end;
    SetLength(Pivots, Min(M, N));
    SetLength(Tmp, 2*Max(M, N));
    CMatrixLUPRec(A, 0, M, N, Pivots, Tmp);
    if AP_FP_Neq(MX,0) then
    begin
        V := MX;
        I:=0;
        while I<=M-1 do
        begin
            for i_ := 0 to Min(I, N-1) do
            begin
                A[I,i_] := C_MulR(A[I,i_],V);
            end;
            Inc(I);
        end;
    end;
end;


procedure RMatrixPLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    Tmp : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MX : Double;
    V : Double;
begin
    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    Assert(M>0, 'RMatrixPLU: incorrect M!');
    Assert(N>0, 'RMatrixPLU: incorrect N!');
    SetLength(Tmp, 2*Max(M, N));
    SetLength(Pivots, Min(M, N));
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    MX := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            MX := Max(MX, AbsReal(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Neq(MX,0) then
    begin
        V := 1/MX;
        I:=0;
        while I<=M-1 do
        begin
            APVMul(@A[I][0], 0, N-1, V);
            Inc(I);
        end;
    end;
    RMatrixPLURec(A, 0, M, N, Pivots, Tmp);
    if AP_FP_Neq(MX,0) then
    begin
        V := MX;
        I:=0;
        while I<=Min(M, N)-1 do
        begin
            APVMul(@A[I][0], I, N-1, V);
            Inc(I);
        end;
    end;
end;


procedure CMatrixPLU(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    Tmp : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MX : Double;
    V : Complex;
    i_ : AlglibInteger;
begin
    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    Assert(M>0, 'CMatrixPLU: incorrect M!');
    Assert(N>0, 'CMatrixPLU: incorrect N!');
    SetLength(Tmp, 2*Max(M, N));
    SetLength(Pivots, Min(M, N));
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    MX := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            MX := Max(MX, AbsComplex(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Neq(MX,0) then
    begin
        V := C_Complex(1/MX);
        I:=0;
        while I<=M-1 do
        begin
            for i_ := 0 to N-1 do
            begin
                A[I,i_] := C_Mul(V, A[I,i_]);
            end;
            Inc(I);
        end;
    end;
    CMatrixPLURec(A, 0, M, N, Pivots, Tmp);
    if AP_FP_Neq(MX,0) then
    begin
        V := C_Complex(MX);
        I:=0;
        while I<=Min(M, N)-1 do
        begin
            for i_ := I to N-1 do
            begin
                A[I,i_] := C_Mul(V, A[I,i_]);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Recurrent complex LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixLUPRec(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);
var
    I : AlglibInteger;
    M1 : AlglibInteger;
    M2 : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Kernel case
    //
    if Min(M, N)<=ABLASComplexBlockSize(A) then
    begin
        CMatrixLUP2(A, Offs, M, N, Pivots, Tmp);
        Exit;
    end;
    
    //
    // Preliminary step, make N>=M
    //
    //     ( A1 )
    // A = (    ), where A1 is square
    //     ( A2 )
    //
    // Factorize A1, update A2
    //
    if M>N then
    begin
        CMatrixLUPRec(A, Offs, N, N, Pivots, Tmp);
        I:=0;
        while I<=N-1 do
        begin
            i1_ := (Offs+N) - (0);
            for i_ := 0 to M-N-1 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+I];
            end;
            for i_ := Offs+N to Offs+M-1 do
            begin
                A[i_,Offs+I] := A[i_,Pivots[Offs+I]];
            end;
            i1_ := (0) - (Offs+N);
            for i_ := Offs+N to Offs+M-1 do
            begin
                A[i_,Pivots[Offs+I]] := Tmp[i_+i1_];
            end;
            Inc(I);
        end;
        CMatrixRightTRSM(M-N, N, A, Offs, Offs, True, True, 0, A, Offs+N, Offs);
        Exit;
    end;
    
    //
    // Non-kernel case
    //
    ABLASComplexSplitLength(A, M, M1, M2);
    CMatrixLUPRec(A, Offs, M1, N, Pivots, Tmp);
    if M2>0 then
    begin
        I:=0;
        while I<=M1-1 do
        begin
            if Offs+I<>Pivots[Offs+I] then
            begin
                i1_ := (Offs+M1) - (0);
                for i_ := 0 to M2-1 do
                begin
                    Tmp[i_] := A[i_+i1_,Offs+I];
                end;
                for i_ := Offs+M1 to Offs+M-1 do
                begin
                    A[i_,Offs+I] := A[i_,Pivots[Offs+I]];
                end;
                i1_ := (0) - (Offs+M1);
                for i_ := Offs+M1 to Offs+M-1 do
                begin
                    A[i_,Pivots[Offs+I]] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
        CMatrixRightTRSM(M2, M1, A, Offs, Offs, True, True, 0, A, Offs+M1, Offs);
        CMatrixGEMM(M-M1, N-M1, M1, C_Complex(-Double(1.0)), A, Offs+M1, Offs, 0, A, Offs, Offs+M1, 0, C_Complex(+Double(1.0)), A, Offs+M1, Offs+M1);
        CMatrixLUPRec(A, Offs+M1, M-M1, N-M1, Pivots, Tmp);
        I:=0;
        while I<=M2-1 do
        begin
            if Offs+M1+I<>Pivots[Offs+M1+I] then
            begin
                i1_ := (Offs) - (0);
                for i_ := 0 to M1-1 do
                begin
                    Tmp[i_] := A[i_+i1_,Offs+M1+I];
                end;
                for i_ := Offs to Offs+M1-1 do
                begin
                    A[i_,Offs+M1+I] := A[i_,Pivots[Offs+M1+I]];
                end;
                i1_ := (0) - (Offs);
                for i_ := Offs to Offs+M1-1 do
                begin
                    A[i_,Pivots[Offs+M1+I]] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Recurrent real LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixLUPRec(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);
var
    I : AlglibInteger;
    M1 : AlglibInteger;
    M2 : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Kernel case
    //
    if Min(M, N)<=ABLASBlockSize(A) then
    begin
        RMatrixLUP2(A, Offs, M, N, Pivots, Tmp);
        Exit;
    end;
    
    //
    // Preliminary step, make N>=M
    //
    //     ( A1 )
    // A = (    ), where A1 is square
    //     ( A2 )
    //
    // Factorize A1, update A2
    //
    if M>N then
    begin
        RMatrixLUPRec(A, Offs, N, N, Pivots, Tmp);
        I:=0;
        while I<=N-1 do
        begin
            if Offs+I<>Pivots[Offs+I] then
            begin
                i1_ := (Offs+N) - (0);
                for i_ := 0 to M-N-1 do
                begin
                    Tmp[i_] := A[i_+i1_,Offs+I];
                end;
                for i_ := Offs+N to Offs+M-1 do
                begin
                    A[i_,Offs+I] := A[i_,Pivots[Offs+I]];
                end;
                i1_ := (0) - (Offs+N);
                for i_ := Offs+N to Offs+M-1 do
                begin
                    A[i_,Pivots[Offs+I]] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
        RMatrixRightTRSM(M-N, N, A, Offs, Offs, True, True, 0, A, Offs+N, Offs);
        Exit;
    end;
    
    //
    // Non-kernel case
    //
    ABLASSplitLength(A, M, M1, M2);
    RMatrixLUPRec(A, Offs, M1, N, Pivots, Tmp);
    if M2>0 then
    begin
        I:=0;
        while I<=M1-1 do
        begin
            if Offs+I<>Pivots[Offs+I] then
            begin
                i1_ := (Offs+M1) - (0);
                for i_ := 0 to M2-1 do
                begin
                    Tmp[i_] := A[i_+i1_,Offs+I];
                end;
                for i_ := Offs+M1 to Offs+M-1 do
                begin
                    A[i_,Offs+I] := A[i_,Pivots[Offs+I]];
                end;
                i1_ := (0) - (Offs+M1);
                for i_ := Offs+M1 to Offs+M-1 do
                begin
                    A[i_,Pivots[Offs+I]] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
        RMatrixRightTRSM(M2, M1, A, Offs, Offs, True, True, 0, A, Offs+M1, Offs);
        RMatrixGEMM(M-M1, N-M1, M1, -Double(1.0), A, Offs+M1, Offs, 0, A, Offs, Offs+M1, 0, +Double(1.0), A, Offs+M1, Offs+M1);
        RMatrixLUPRec(A, Offs+M1, M-M1, N-M1, Pivots, Tmp);
        I:=0;
        while I<=M2-1 do
        begin
            if Offs+M1+I<>Pivots[Offs+M1+I] then
            begin
                i1_ := (Offs) - (0);
                for i_ := 0 to M1-1 do
                begin
                    Tmp[i_] := A[i_+i1_,Offs+M1+I];
                end;
                for i_ := Offs to Offs+M1-1 do
                begin
                    A[i_,Offs+M1+I] := A[i_,Pivots[Offs+M1+I]];
                end;
                i1_ := (0) - (Offs);
                for i_ := Offs to Offs+M1-1 do
                begin
                    A[i_,Pivots[Offs+M1+I]] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Recurrent complex LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixPLURec(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);
var
    I : AlglibInteger;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Kernel case
    //
    if Min(M, N)<=ABLASComplexBlockSize(A) then
    begin
        CMatrixPLU2(A, Offs, M, N, Pivots, Tmp);
        Exit;
    end;
    
    //
    // Preliminary step, make M>=N.
    //
    // A = (A1 A2), where A1 is square
    // Factorize A1, update A2
    //
    if N>M then
    begin
        CMatrixPLURec(A, Offs, M, M, Pivots, Tmp);
        I:=0;
        while I<=M-1 do
        begin
            i1_ := (Offs+M) - (0);
            for i_ := 0 to N-M-1 do
            begin
                Tmp[i_] := A[Offs+I,i_+i1_];
            end;
            for i_ := Offs+M to Offs+N-1 do
            begin
                A[Offs+I,i_] := A[Pivots[Offs+I],i_];
            end;
            i1_ := (0) - (Offs+M);
            for i_ := Offs+M to Offs+N-1 do
            begin
                A[Pivots[Offs+I],i_] := Tmp[i_+i1_];
            end;
            Inc(I);
        end;
        CMatrixLeftTRSM(M, N-M, A, Offs, Offs, False, True, 0, A, Offs, Offs+M);
        Exit;
    end;
    
    //
    // Non-kernel case
    //
    ABLASComplexSplitLength(A, N, N1, N2);
    CMatrixPLURec(A, Offs, M, N1, Pivots, Tmp);
    if N2>0 then
    begin
        I:=0;
        while I<=N1-1 do
        begin
            if Offs+I<>Pivots[Offs+I] then
            begin
                i1_ := (Offs+N1) - (0);
                for i_ := 0 to N2-1 do
                begin
                    Tmp[i_] := A[Offs+I,i_+i1_];
                end;
                for i_ := Offs+N1 to Offs+N-1 do
                begin
                    A[Offs+I,i_] := A[Pivots[Offs+I],i_];
                end;
                i1_ := (0) - (Offs+N1);
                for i_ := Offs+N1 to Offs+N-1 do
                begin
                    A[Pivots[Offs+I],i_] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
        CMatrixLeftTRSM(N1, N2, A, Offs, Offs, False, True, 0, A, Offs, Offs+N1);
        CMatrixGEMM(M-N1, N-N1, N1, C_Complex(-Double(1.0)), A, Offs+N1, Offs, 0, A, Offs, Offs+N1, 0, C_Complex(+Double(1.0)), A, Offs+N1, Offs+N1);
        CMatrixPLURec(A, Offs+N1, M-N1, N-N1, Pivots, Tmp);
        I:=0;
        while I<=N2-1 do
        begin
            if Offs+N1+I<>Pivots[Offs+N1+I] then
            begin
                i1_ := (Offs) - (0);
                for i_ := 0 to N1-1 do
                begin
                    Tmp[i_] := A[Offs+N1+I,i_+i1_];
                end;
                for i_ := Offs to Offs+N1-1 do
                begin
                    A[Offs+N1+I,i_] := A[Pivots[Offs+N1+I],i_];
                end;
                i1_ := (0) - (Offs);
                for i_ := Offs to Offs+N1-1 do
                begin
                    A[Pivots[Offs+N1+I],i_] := Tmp[i_+i1_];
                end;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Recurrent real LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixPLURec(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);
var
    I : AlglibInteger;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
begin
    
    //
    // Kernel case
    //
    if Min(M, N)<=ABLASBlockSize(A) then
    begin
        RMatrixPLU2(A, Offs, M, N, Pivots, Tmp);
        Exit;
    end;
    
    //
    // Preliminary step, make M>=N.
    //
    // A = (A1 A2), where A1 is square
    // Factorize A1, update A2
    //
    if N>M then
    begin
        RMatrixPLURec(A, Offs, M, M, Pivots, Tmp);
        I:=0;
        while I<=M-1 do
        begin
            APVMove(@Tmp[0], 0, N-M-1, @A[Offs+I][0], Offs+M, Offs+N-1);
            APVMove(@A[Offs+I][0], Offs+M, Offs+N-1, @A[Pivots[Offs+I]][0], Offs+M, Offs+N-1);
            APVMove(@A[Pivots[Offs+I]][0], Offs+M, Offs+N-1, @Tmp[0], 0, N-M-1);
            Inc(I);
        end;
        RMatrixLeftTRSM(M, N-M, A, Offs, Offs, False, True, 0, A, Offs, Offs+M);
        Exit;
    end;
    
    //
    // Non-kernel case
    //
    ABLASSplitLength(A, N, N1, N2);
    RMatrixPLURec(A, Offs, M, N1, Pivots, Tmp);
    if N2>0 then
    begin
        I:=0;
        while I<=N1-1 do
        begin
            if Offs+I<>Pivots[Offs+I] then
            begin
                APVMove(@Tmp[0], 0, N2-1, @A[Offs+I][0], Offs+N1, Offs+N-1);
                APVMove(@A[Offs+I][0], Offs+N1, Offs+N-1, @A[Pivots[Offs+I]][0], Offs+N1, Offs+N-1);
                APVMove(@A[Pivots[Offs+I]][0], Offs+N1, Offs+N-1, @Tmp[0], 0, N2-1);
            end;
            Inc(I);
        end;
        RMatrixLeftTRSM(N1, N2, A, Offs, Offs, False, True, 0, A, Offs, Offs+N1);
        RMatrixGEMM(M-N1, N-N1, N1, -Double(1.0), A, Offs+N1, Offs, 0, A, Offs, Offs+N1, 0, +Double(1.0), A, Offs+N1, Offs+N1);
        RMatrixPLURec(A, Offs+N1, M-N1, N-N1, Pivots, Tmp);
        I:=0;
        while I<=N2-1 do
        begin
            if Offs+N1+I<>Pivots[Offs+N1+I] then
            begin
                APVMove(@Tmp[0], 0, N1-1, @A[Offs+N1+I][0], Offs, Offs+N1-1);
                APVMove(@A[Offs+N1+I][0], Offs, Offs+N1-1, @A[Pivots[Offs+N1+I]][0], Offs, Offs+N1-1);
                APVMove(@A[Pivots[Offs+N1+I]][0], Offs, Offs+N1-1, @Tmp[0], 0, N1-1);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Complex LUP kernel

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixLUP2(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    S : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    
    //
    // main cycle
    //
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot, swap columns
        //
        JP := J;
        I:=J+1;
        while I<=N-1 do
        begin
            if AP_FP_Greater(AbsComplex(A[Offs+J,Offs+I]),AbsComplex(A[Offs+J,Offs+JP])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[Offs+J] := Offs+JP;
        if JP<>J then
        begin
            i1_ := (Offs) - (0);
            for i_ := 0 to M-1 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+J];
            end;
            for i_ := Offs to Offs+M-1 do
            begin
                A[i_,Offs+J] := A[i_,Offs+JP];
            end;
            i1_ := (0) - (Offs);
            for i_ := Offs to Offs+M-1 do
            begin
                A[i_,Offs+JP] := Tmp[i_+i1_];
            end;
        end;
        
        //
        // LU decomposition of 1x(N-J) matrix
        //
        if C_NotEqualR(A[Offs+J,Offs+J],0) and (J+1<=N-1) then
        begin
            S := C_RDiv(1,A[Offs+J,Offs+J]);
            for i_ := Offs+J+1 to Offs+N-1 do
            begin
                A[Offs+J,i_] := C_Mul(S, A[Offs+J,i_]);
            end;
        end;
        
        //
        // Update trailing (M-J-1)x(N-J-1) matrix
        //
        if J<Min(M-1, N-1) then
        begin
            i1_ := (Offs+J+1) - (0);
            for i_ := 0 to M-J-2 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+J];
            end;
            i1_ := (Offs+J+1) - (M);
            for i_ := M to M+N-J-2 do
            begin
                Tmp[i_] := C_Opposite(A[Offs+J,i_+i1_]);
            end;
            CMatrixRank1(M-J-1, N-J-1, A, Offs+J+1, Offs+J+1, Tmp, 0, Tmp, M);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Real LUP kernel

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixLUP2(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    S : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    
    //
    // main cycle
    //
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot, swap columns
        //
        JP := J;
        I:=J+1;
        while I<=N-1 do
        begin
            if AP_FP_Greater(AbsReal(A[Offs+J,Offs+I]),AbsReal(A[Offs+J,Offs+JP])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[Offs+J] := Offs+JP;
        if JP<>J then
        begin
            i1_ := (Offs) - (0);
            for i_ := 0 to M-1 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+J];
            end;
            for i_ := Offs to Offs+M-1 do
            begin
                A[i_,Offs+J] := A[i_,Offs+JP];
            end;
            i1_ := (0) - (Offs);
            for i_ := Offs to Offs+M-1 do
            begin
                A[i_,Offs+JP] := Tmp[i_+i1_];
            end;
        end;
        
        //
        // LU decomposition of 1x(N-J) matrix
        //
        if AP_FP_Neq(A[Offs+J,Offs+J],0) and (J+1<=N-1) then
        begin
            S := 1/A[Offs+J,Offs+J];
            APVMul(@A[Offs+J][0], Offs+J+1, Offs+N-1, S);
        end;
        
        //
        // Update trailing (M-J-1)x(N-J-1) matrix
        //
        if J<Min(M-1, N-1) then
        begin
            i1_ := (Offs+J+1) - (0);
            for i_ := 0 to M-J-2 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+J];
            end;
            APVMoveNeg(@Tmp[0], M, M+N-J-2, @A[Offs+J][0], Offs+J+1, Offs+N-1);
            RMatrixRank1(M-J-1, N-J-1, A, Offs+J+1, Offs+J+1, Tmp, 0, Tmp, M);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Complex PLU kernel

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************)
procedure CMatrixPLU2(var A : TComplex2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TComplex1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    S : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M-1 do
        begin
            if AP_FP_Greater(AbsComplex(A[Offs+I,Offs+J]),AbsComplex(A[Offs+JP,Offs+J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[Offs+J] := Offs+JP;
        if C_NotEqualR(A[Offs+JP,Offs+J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    S := A[Offs+J,Offs+I];
                    A[Offs+J,Offs+I] := A[Offs+JP,Offs+I];
                    A[Offs+JP,Offs+I] := S;
                    Inc(I);
                end;
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J+1<=M-1 then
            begin
                S := C_RDiv(1,A[Offs+J,Offs+J]);
                for i_ := Offs+J+1 to Offs+M-1 do
                begin
                    A[i_,Offs+J] := C_Mul(S, A[i_,Offs+J]);
                end;
            end;
        end;
        if J<Min(M, N)-1 then
        begin
            
            //
            //Update trailing submatrix.
            //
            i1_ := (Offs+J+1) - (0);
            for i_ := 0 to M-J-2 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+J];
            end;
            i1_ := (Offs+J+1) - (M);
            for i_ := M to M+N-J-2 do
            begin
                Tmp[i_] := C_Opposite(A[Offs+J,i_+i1_]);
            end;
            CMatrixRank1(M-J-1, N-J-1, A, Offs+J+1, Offs+J+1, Tmp, 0, Tmp, M);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Real PLU kernel

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************)
procedure RMatrixPLU2(var A : TReal2DArray;
     Offs : AlglibInteger;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray;
     var Tmp : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    S : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M-1 do
        begin
            if AP_FP_Greater(AbsReal(A[Offs+I,Offs+J]),AbsReal(A[Offs+JP,Offs+J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[Offs+J] := Offs+JP;
        if AP_FP_Neq(A[Offs+JP,Offs+J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    S := A[Offs+J,Offs+I];
                    A[Offs+J,Offs+I] := A[Offs+JP,Offs+I];
                    A[Offs+JP,Offs+I] := S;
                    Inc(I);
                end;
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J+1<=M-1 then
            begin
                S := 1/A[Offs+J,Offs+J];
                for i_ := Offs+J+1 to Offs+M-1 do
                begin
                    A[i_,Offs+J] := S*A[i_,Offs+J];
                end;
            end;
        end;
        if J<Min(M, N)-1 then
        begin
            
            //
            //Update trailing submatrix.
            //
            i1_ := (Offs+J+1) - (0);
            for i_ := 0 to M-J-2 do
            begin
                Tmp[i_] := A[i_+i1_,Offs+J];
            end;
            APVMoveNeg(@Tmp[0], M, M+N-J-2, @A[Offs+J][0], Offs+J+1, Offs+N-1);
            RMatrixRank1(M-J-1, N-J-1, A, Offs+J+1, Offs+J+1, Tmp, 0, Tmp, M);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Recursive computational subroutine for HPDMatrixCholesky

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
function HPDMatrixCholeskyRec(var A : TComplex2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TComplex1DArray):Boolean;
var
    N1 : AlglibInteger;
    N2 : AlglibInteger;
begin
    
    //
    // check N
    //
    if N<1 then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // special cases
    //
    if N=1 then
    begin
        if AP_FP_Greater(A[Offs,Offs].X,0) then
        begin
            A[Offs,Offs] := C_Complex(Sqrt(A[Offs,Offs].X));
            Result := True;
        end
        else
        begin
            Result := False;
        end;
        Exit;
    end;
    if N<=ABLASComplexBlockSize(A) then
    begin
        Result := HPDMatrixCholesky2(A, Offs, N, IsUpper, Tmp);
        Exit;
    end;
    
    //
    // general case: split task in cache-oblivious manner
    //
    Result := True;
    ABLASComplexSplitLength(A, N, N1, N2);
    Result := HPDMatrixCholeskyRec(A, Offs, N1, IsUpper, Tmp);
    if  not Result then
    begin
        Exit;
    end;
    if N2>0 then
    begin
        if IsUpper then
        begin
            CMatrixLeftTRSM(N1, N2, A, Offs, Offs, IsUpper, False, 2, A, Offs, Offs+N1);
            CMatrixSYRK(N2, N1, -Double(1.0), A, Offs, Offs+N1, 2, +Double(1.0), A, Offs+N1, Offs+N1, IsUpper);
        end
        else
        begin
            CMatrixRightTRSM(N2, N1, A, Offs, Offs, IsUpper, False, 2, A, Offs+N1, Offs);
            CMatrixSYRK(N2, N1, -Double(1.0), A, Offs+N1, Offs, 0, +Double(1.0), A, Offs+N1, Offs+N1, IsUpper);
        end;
        Result := HPDMatrixCholeskyRec(A, Offs+N1, N2, IsUpper, Tmp);
        if  not Result then
        begin
            Exit;
        end;
    end;
end;


(*************************************************************************
Recursive computational subroutine for SPDMatrixCholesky

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
function SPDMatrixCholeskyRec(var A : TReal2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TReal1DArray):Boolean;
var
    N1 : AlglibInteger;
    N2 : AlglibInteger;
begin
    
    //
    // check N
    //
    if N<1 then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // special cases
    //
    if N=1 then
    begin
        if AP_FP_Greater(A[Offs,Offs],0) then
        begin
            A[Offs,Offs] := Sqrt(A[Offs,Offs]);
            Result := True;
        end
        else
        begin
            Result := False;
        end;
        Exit;
    end;
    if N<=ABLASBlockSize(A) then
    begin
        Result := SPDMatrixCholesky2(A, Offs, N, IsUpper, Tmp);
        Exit;
    end;
    
    //
    // general case: split task in cache-oblivious manner
    //
    Result := True;
    ABLASSplitLength(A, N, N1, N2);
    Result := SPDMatrixCholeskyRec(A, Offs, N1, IsUpper, Tmp);
    if  not Result then
    begin
        Exit;
    end;
    if N2>0 then
    begin
        if IsUpper then
        begin
            RMatrixLeftTRSM(N1, N2, A, Offs, Offs, IsUpper, False, 1, A, Offs, Offs+N1);
            RMatrixSYRK(N2, N1, -Double(1.0), A, Offs, Offs+N1, 1, +Double(1.0), A, Offs+N1, Offs+N1, IsUpper);
        end
        else
        begin
            RMatrixRightTRSM(N2, N1, A, Offs, Offs, IsUpper, False, 1, A, Offs+N1, Offs);
            RMatrixSYRK(N2, N1, -Double(1.0), A, Offs+N1, Offs, 0, +Double(1.0), A, Offs+N1, Offs+N1, IsUpper);
        end;
        Result := SPDMatrixCholeskyRec(A, Offs+N1, N2, IsUpper, Tmp);
        if  not Result then
        begin
            Exit;
        end;
    end;
end;


(*************************************************************************
Level-2 Hermitian Cholesky subroutine.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function HPDMatrixCholesky2(var AAA : TComplex2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TComplex1DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    AJJ : Double;
    V : Complex;
    R : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Result := True;
    if N<0 then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            V := C_Complex(0.0);
            for i_ := Offs to Offs+J-1 do
            begin
                V := C_Add(V,C_Mul(Conj(AAA[i_,Offs+J]),AAA[i_,Offs+J]));
            end;
            AJJ := C_Sub(AAA[Offs+J,Offs+J],V).X;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                AAA[Offs+J,Offs+J] := C_Complex(AJJ);
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            AAA[Offs+J,Offs+J] := C_Complex(AJJ);
            
            //
            // Compute elements J+1:N-1 of row J.
            //
            if J<N-1 then
            begin
                if J>0 then
                begin
                    i1_ := (Offs) - (0);
                    for i_ := 0 to J-1 do
                    begin
                        Tmp[i_] := C_Opposite(Conj(AAA[i_+i1_,Offs+J]));
                    end;
                    CMatrixMV(N-J-1, J, AAA, Offs, Offs+J+1, 1, Tmp, 0, Tmp, N);
                    i1_ := (N) - (Offs+J+1);
                    for i_ := Offs+J+1 to Offs+N-1 do
                    begin
                        AAA[Offs+J,i_] := C_Add(AAA[Offs+J,i_], Tmp[i_+i1_]);
                    end;
                end;
                R := 1/AJJ;
                for i_ := Offs+J+1 to Offs+N-1 do
                begin
                    AAA[Offs+J,i_] := C_MulR(AAA[Offs+J,i_],R);
                end;
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute L(J+1,J+1) and test for non-positive-definiteness.
            //
            V := C_Complex(0.0);
            for i_ := Offs to Offs+J-1 do
            begin
                V := C_Add(V,C_Mul(Conj(AAA[Offs+J,i_]),AAA[Offs+J,i_]));
            end;
            AJJ := C_Sub(AAA[Offs+J,Offs+J],V).X;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                AAA[Offs+J,Offs+J] := C_Complex(AJJ);
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            AAA[Offs+J,Offs+J] := C_Complex(AJJ);
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N-1 then
            begin
                if J>0 then
                begin
                    i1_ := (Offs) - (0);
                    for i_ := 0 to J-1 do
                    begin
                        Tmp[i_] := Conj(AAA[Offs+J,i_+i1_]);
                    end;
                    CMatrixMV(N-J-1, J, AAA, Offs+J+1, Offs, 0, Tmp, 0, Tmp, N);
                    I:=0;
                    while I<=N-J-2 do
                    begin
                        AAA[Offs+J+1+I,Offs+J] := C_DivR(C_Sub(AAA[Offs+J+1+I,Offs+J],Tmp[N+I]),AJJ);
                        Inc(I);
                    end;
                end
                else
                begin
                    I:=0;
                    while I<=N-J-2 do
                    begin
                        AAA[Offs+J+1+I,Offs+J] := C_DivR(AAA[Offs+J+1+I,Offs+J],AJJ);
                        Inc(I);
                    end;
                end;
            end;
            Inc(J);
        end;
    end;
end;


(*************************************************************************
Level-2 Cholesky subroutine

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function SPDMatrixCholesky2(var AAA : TReal2DArray;
     Offs : AlglibInteger;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tmp : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    AJJ : Double;
    V : Double;
    R : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Result := True;
    if N<0 then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            V := 0.0;
            for i_ := Offs to Offs+J-1 do
            begin
                V := V + AAA[i_,Offs+J]*AAA[i_,Offs+J];
            end;
            AJJ := AAA[Offs+J,Offs+J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                AAA[Offs+J,Offs+J] := AJJ;
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            AAA[Offs+J,Offs+J] := AJJ;
            
            //
            // Compute elements J+1:N-1 of row J.
            //
            if J<N-1 then
            begin
                if J>0 then
                begin
                    i1_ := (Offs) - (0);
                    for i_ := 0 to J-1 do
                    begin
                        Tmp[i_] := -AAA[i_+i1_,Offs+J];
                    end;
                    RMatrixMV(N-J-1, J, AAA, Offs, Offs+J+1, 1, Tmp, 0, Tmp, N);
                    APVAdd(@AAA[Offs+J][0], Offs+J+1, Offs+N-1, @Tmp[0], N, 2*N-J-2);
                end;
                R := 1/AJJ;
                APVMul(@AAA[Offs+J][0], Offs+J+1, Offs+N-1, R);
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute L(J+1,J+1) and test for non-positive-definiteness.
            //
            V := APVDotProduct(@AAA[Offs+J][0], Offs, Offs+J-1, @AAA[Offs+J][0], Offs, Offs+J-1);
            AJJ := AAA[Offs+J,Offs+J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                AAA[Offs+J,Offs+J] := AJJ;
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            AAA[Offs+J,Offs+J] := AJJ;
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N-1 then
            begin
                if J>0 then
                begin
                    APVMove(@Tmp[0], 0, J-1, @AAA[Offs+J][0], Offs, Offs+J-1);
                    RMatrixMV(N-J-1, J, AAA, Offs+J+1, Offs, 0, Tmp, 0, Tmp, N);
                    I:=0;
                    while I<=N-J-2 do
                    begin
                        AAA[Offs+J+1+I,Offs+J] := (AAA[Offs+J+1+I,Offs+J]-Tmp[N+I])/AJJ;
                        Inc(I);
                    end;
                end
                else
                begin
                    I:=0;
                    while I<=N-J-2 do
                    begin
                        AAA[Offs+J+1+I,Offs+J] := AAA[Offs+J+1+I,Offs+J]/AJJ;
                        Inc(I);
                    end;
                end;
            end;
            Inc(J);
        end;
    end;
end;


end.