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
unit schur;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur;

function RMatrixSchur(var A : TReal2DArray;
     N : AlglibInteger;
     var S : TReal2DArray):Boolean;

implementation

(*************************************************************************
Subroutine performing the Schur decomposition of a general matrix by using
the QR algorithm with multiple shifts.

The source matrix A is represented as S'*A*S = T, where S is an orthogonal
matrix (Schur vectors), T - upper quasi-triangular matrix (with blocks of
sizes 1x1 and 2x2 on the main diagonal).

Input parameters:
    A   -   matrix to be decomposed.
            Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of A, N>=0.


Output parameters:
    A   -   contains matrix T.
            Array whose indexes range within [0..N-1, 0..N-1].
    S   -   contains Schur vectors.
            Array whose indexes range within [0..N-1, 0..N-1].

Note 1:
    The block structure of matrix T can be easily recognized: since all
    the elements below the blocks are zeros, the elements a[i+1,i] which
    are equal to 0 show the block border.

Note 2:
    The algorithm performance depends on the value of the internal parameter
    NS of the InternalSchurDecomposition subroutine which defines the number
    of shifts in the QR algorithm (similarly to the block width in block-matrix
    algorithms in linear algebra). If you require maximum performance on
    your machine, it is recommended to adjust this parameter manually.

Result:
    True,
        if the algorithm has converged and parameters A and S contain the result.
    False,
        if the algorithm has not converged.

Algorithm implemented on the basis of the DHSEQR subroutine (LAPACK 3.0 library).
*************************************************************************)
function RMatrixSchur(var A : TReal2DArray;
     N : AlglibInteger;
     var S : TReal2DArray):Boolean;
var
    Tau : TReal1DArray;
    WI : TReal1DArray;
    WR : TReal1DArray;
    A1 : TReal2DArray;
    S1 : TReal2DArray;
    INFO : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
begin
    
    //
    // Upper Hessenberg form of the 0-based matrix
    //
    RMatrixHessenberg(A, N, Tau);
    RMatrixHessenbergUnpackQ(A, N, Tau, S);
    
    //
    // Convert from 0-based arrays to 1-based,
    // then call InternalSchurDecomposition
    // Awkward, of course, but Schur decompisiton subroutine
    // is too complex to fix it.
    //
    //
    SetLength(A1, N+1, N+1);
    SetLength(S1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A1[I,J] := A[I-1,J-1];
            S1[I,J] := S[I-1,J-1];
            Inc(J);
        end;
        Inc(I);
    end;
    InternalSchurDecomposition(A1, N, 1, 1, WR, WI, S1, INFO);
    Result := INFO=0;
    
    //
    // convert from 1-based arrays to -based
    //
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A[I-1,J-1] := A1[I,J];
            S[I-1,J-1] := S1[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


end.