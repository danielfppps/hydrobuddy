{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2005-2010 Sergey Bochkanov.

Additional copyrights:
    1992-2007 The University of Tennessee (as indicated in subroutines
    comments).

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
unit ortfac;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas;

procedure RMatrixQR(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure RMatrixLQ(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure CMatrixQR(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
procedure CMatrixLQ(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
procedure RMatrixQRUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
procedure RMatrixQRUnpackR(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TReal2DArray);
procedure RMatrixLQUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QRows : AlglibInteger;
     var Q : TReal2DArray);
procedure RMatrixLQUnpackL(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray);
procedure CMatrixQRUnpackQ(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QColumns : AlglibInteger;
     var Q : TComplex2DArray);
procedure CMatrixQRUnpackR(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TComplex2DArray);
procedure CMatrixLQUnpackQ(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QRows : AlglibInteger;
     var Q : TComplex2DArray);
procedure CMatrixLQUnpackL(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TComplex2DArray);
procedure RMatrixBD(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var TauQ : TReal1DArray;
     var TauP : TReal1DArray);
procedure RMatrixBDUnpackQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
procedure RMatrixBDMultiplyByQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
procedure RMatrixBDUnpackPT(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     PTRows : AlglibInteger;
     var PT : TReal2DArray);
procedure RMatrixBDMultiplyByP(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
procedure RMatrixBDUnpackDiagonals(const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var IsUpper : Boolean;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure RMatrixHessenberg(var A : TReal2DArray;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure RMatrixHessenbergUnpackQ(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
procedure RMatrixHessenbergUnpackH(const A : TReal2DArray;
     N : AlglibInteger;
     var H : TReal2DArray);
procedure SMatrixTD(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TReal1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure SMatrixTDUnpackQ(const A : TReal2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
procedure HMatrixTD(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TComplex1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure HMatrixTDUnpackQ(const A : TComplex2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TComplex1DArray;
     var Q : TComplex2DArray);

implementation

procedure RMatrixQRBaseCase(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TReal1DArray;
     var T : TReal1DArray;
     var Tau : TReal1DArray);forward;
procedure RMatrixLQBaseCase(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TReal1DArray;
     var T : TReal1DArray;
     var Tau : TReal1DArray);forward;
procedure CMatrixQRBaseCase(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TComplex1DArray;
     var T : TComplex1DArray;
     var Tau : TComplex1DArray);forward;
procedure CMatrixLQBaseCase(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TComplex1DArray;
     var T : TComplex1DArray;
     var Tau : TComplex1DArray);forward;
procedure RMatrixBlockReflector(var A : TReal2DArray;
     var Tau : TReal1DArray;
     ColumnwiseA : Boolean;
     LengthA : AlglibInteger;
     BlockSize : AlglibInteger;
     var T : TReal2DArray;
     var WORK : TReal1DArray);forward;
procedure CMatrixBlockReflector(var A : TComplex2DArray;
     var Tau : TComplex1DArray;
     ColumnwiseA : Boolean;
     LengthA : AlglibInteger;
     BlockSize : AlglibInteger;
     var T : TComplex2DArray;
     var WORK : TComplex1DArray);forward;


(*************************************************************************
QR decomposition of a rectangular matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and R in compact form (see below).
    Tau -   array of scalar factors which are used to form
            matrix Q. Array whose index ranges within [0.. Min(M-1,N-1)].

Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
MxM, R - upper triangular (or upper trapezoid) matrix of size M x N.

The elements of matrix R are located on and above the main diagonal of
matrix A. The elements which are located in Tau array and below the main
diagonal of matrix A are used to form matrix Q as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(0)*H(2)*...*H(k-1),

where k = min(m,n), and each H(i) is in the form

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - real vector,
so that v(0:i-1) = 0, v(i) = 1, v(i+1:m-1) stored in A(i+1:m-1,i).

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixQR(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    TauBuf : TReal1DArray;
    MinMN : AlglibInteger;
    TmpA : TReal2DArray;
    TmpT : TReal2DArray;
    TmpR : TReal2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    RowsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(Tau, MinMN);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, M, ABLASBlockSize(A));
    SetLength(TmpT, ABLASBlockSize(A), 2*ABLASBlockSize(A));
    SetLength(TmpR, 2*ABLASBlockSize(A), N);
    
    //
    // Blocked code
    //
    BlockStart := 0;
    while BlockStart<>MinMN do
    begin
        
        //
        // Determine block size
        //
        BlockSize := MinMN-BlockStart;
        if BlockSize>ABLASBlockSize(A) then
        begin
            BlockSize := ABLASBlockSize(A);
        end;
        RowsCount := M-BlockStart;
        
        //
        // QR decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        RMatrixCopy(RowsCount, BlockSize, A, BlockStart, BlockStart, TmpA, 0, 0);
        RMatrixQRBaseCase(TmpA, RowsCount, BlockSize, WORK, T, TauBuf);
        RMatrixCopy(RowsCount, BlockSize, TmpA, 0, 0, A, BlockStart, BlockStart);
        APVMove(@Tau[0], BlockStart, BlockStart+BlockSize-1, @TauBuf[0], 0, BlockSize-1);
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if BlockStart+BlockSize<=N-1 then
        begin
            if (N-BlockStart-BlockSize>=2*ABLASBlockSize(A)) or (RowsCount>=4*ABLASBlockSize(A)) then
            begin
                
                //
                // Prepare block reflector
                //
                RMatrixBlockReflector(TmpA, TauBuf, True, RowsCount, BlockSize, TmpT, WORK);
                
                //
                // Multiply the rest of A by Q'.
                //
                // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
                // Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
                //
                RMatrixGEMM(BlockSize, N-BlockStart-BlockSize, RowsCount, Double(1.0), TmpA, 0, 0, 1, A, BlockStart, BlockStart+BlockSize, 0, Double(0.0), TmpR, 0, 0);
                RMatrixGEMM(BlockSize, N-BlockStart-BlockSize, BlockSize, Double(1.0), TmpT, 0, 0, 1, TmpR, 0, 0, 0, Double(0.0), TmpR, BlockSize, 0);
                RMatrixGEMM(RowsCount, N-BlockStart-BlockSize, BlockSize, Double(1.0), TmpA, 0, 0, 0, TmpR, BlockSize, 0, 0, Double(1.0), A, BlockStart, BlockStart+BlockSize);
            end
            else
            begin
                
                //
                // Level 2 algorithm
                //
                I:=0;
                while I<=BlockSize-1 do
                begin
                    i1_ := (I) - (1);
                    for i_ := 1 to RowsCount-I do
                    begin
                        T[i_] := TmpA[i_+i1_,I];
                    end;
                    T[1] := 1;
                    ApplyReflectionFromTheLeft(A, TauBuf[I], T, BlockStart+I, M-1, BlockStart+BlockSize, N-1, WORK);
                    Inc(I);
                end;
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart+BlockSize;
    end;
end;


(*************************************************************************
LQ decomposition of a rectangular matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices L and Q in compact form (see below)
    Tau -   array of scalar factors which are used to form
            matrix Q. Array whose index ranges within [0..Min(M,N)-1].

Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
MxM, L - lower triangular (or lower trapezoid) matrix of size M x N.

The elements of matrix L are located on and below  the  main  diagonal  of
matrix A. The elements which are located in Tau array and above  the  main
diagonal of matrix A are used to form matrix Q as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(k-1)*H(k-2)*...*H(1)*H(0),

where k = min(m,n), and each H(i) is of the form

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - real vector, so that v(0:i-1)=0,
v(i) = 1, v(i+1:n-1) stored in A(i,i+1:n-1).

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixLQ(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    TauBuf : TReal1DArray;
    MinMN : AlglibInteger;
    TmpA : TReal2DArray;
    TmpT : TReal2DArray;
    TmpR : TReal2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    ColumnsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(Tau, MinMN);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, ABLASBlockSize(A), N);
    SetLength(TmpT, ABLASBlockSize(A), 2*ABLASBlockSize(A));
    SetLength(TmpR, M, 2*ABLASBlockSize(A));
    
    //
    // Blocked code
    //
    BlockStart := 0;
    while BlockStart<>MinMN do
    begin
        
        //
        // Determine block size
        //
        BlockSize := MinMN-BlockStart;
        if BlockSize>ABLASBlockSize(A) then
        begin
            BlockSize := ABLASBlockSize(A);
        end;
        ColumnsCount := N-BlockStart;
        
        //
        // LQ decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        RMatrixCopy(BlockSize, ColumnsCount, A, BlockStart, BlockStart, TmpA, 0, 0);
        RMatrixLQBaseCase(TmpA, BlockSize, ColumnsCount, WORK, T, TauBuf);
        RMatrixCopy(BlockSize, ColumnsCount, TmpA, 0, 0, A, BlockStart, BlockStart);
        APVMove(@Tau[0], BlockStart, BlockStart+BlockSize-1, @TauBuf[0], 0, BlockSize-1);
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if BlockStart+BlockSize<=M-1 then
        begin
            if M-BlockStart-BlockSize>=2*ABLASBlockSize(A) then
            begin
                
                //
                // Prepare block reflector
                //
                RMatrixBlockReflector(TmpA, TauBuf, False, ColumnsCount, BlockSize, TmpT, WORK);
                
                //
                // Multiply the rest of A by Q.
                //
                // Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
                //
                RMatrixGEMM(M-BlockStart-BlockSize, BlockSize, ColumnsCount, Double(1.0), A, BlockStart+BlockSize, BlockStart, 0, TmpA, 0, 0, 1, Double(0.0), TmpR, 0, 0);
                RMatrixGEMM(M-BlockStart-BlockSize, BlockSize, BlockSize, Double(1.0), TmpR, 0, 0, 0, TmpT, 0, 0, 0, Double(0.0), TmpR, 0, BlockSize);
                RMatrixGEMM(M-BlockStart-BlockSize, ColumnsCount, BlockSize, Double(1.0), TmpR, 0, BlockSize, 0, TmpA, 0, 0, 0, Double(1.0), A, BlockStart+BlockSize, BlockStart);
            end
            else
            begin
                
                //
                // Level 2 algorithm
                //
                I:=0;
                while I<=BlockSize-1 do
                begin
                    APVMove(@T[0], 1, ColumnsCount-I, @TmpA[I][0], I, ColumnsCount-1);
                    T[1] := 1;
                    ApplyReflectionFromTheRight(A, TauBuf[I], T, BlockStart+BlockSize, M-1, BlockStart+I, N-1, WORK);
                    Inc(I);
                end;
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart+BlockSize;
    end;
end;


(*************************************************************************
QR decomposition of a rectangular complex matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and R in compact form
    Tau -   array of scalar factors which are used to form matrix Q. Array
            whose indexes range within [0.. Min(M,N)-1]

Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
MxM, R - upper triangular (or upper trapezoid) matrix of size MxN.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure CMatrixQR(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
var
    WORK : TComplex1DArray;
    T : TComplex1DArray;
    TauBuf : TComplex1DArray;
    MinMN : AlglibInteger;
    TmpA : TComplex2DArray;
    TmpT : TComplex2DArray;
    TmpR : TComplex2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    RowsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(Tau, MinMN);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, M, ABLASComplexBlockSize(A));
    SetLength(TmpT, ABLASComplexBlockSize(A), ABLASComplexBlockSize(A));
    SetLength(TmpR, 2*ABLASComplexBlockSize(A), N);
    
    //
    // Blocked code
    //
    BlockStart := 0;
    while BlockStart<>MinMN do
    begin
        
        //
        // Determine block size
        //
        BlockSize := MinMN-BlockStart;
        if BlockSize>ABLASComplexBlockSize(A) then
        begin
            BlockSize := ABLASComplexBlockSize(A);
        end;
        RowsCount := M-BlockStart;
        
        //
        // QR decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        CMatrixCopy(RowsCount, BlockSize, A, BlockStart, BlockStart, TmpA, 0, 0);
        CMatrixQRBaseCase(TmpA, RowsCount, BlockSize, WORK, T, TauBuf);
        CMatrixCopy(RowsCount, BlockSize, TmpA, 0, 0, A, BlockStart, BlockStart);
        i1_ := (0) - (BlockStart);
        for i_ := BlockStart to BlockStart+BlockSize-1 do
        begin
            Tau[i_] := TauBuf[i_+i1_];
        end;
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if BlockStart+BlockSize<=N-1 then
        begin
            if N-BlockStart-BlockSize>=2*ABLASComplexBlockSize(A) then
            begin
                
                //
                // Prepare block reflector
                //
                CMatrixBlockReflector(TmpA, TauBuf, True, RowsCount, BlockSize, TmpT, WORK);
                
                //
                // Multiply the rest of A by Q'.
                //
                // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
                // Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
                //
                CMatrixGEMM(BlockSize, N-BlockStart-BlockSize, RowsCount, C_Complex(Double(1.0)), TmpA, 0, 0, 2, A, BlockStart, BlockStart+BlockSize, 0, C_Complex(Double(0.0)), TmpR, 0, 0);
                CMatrixGEMM(BlockSize, N-BlockStart-BlockSize, BlockSize, C_Complex(Double(1.0)), TmpT, 0, 0, 2, TmpR, 0, 0, 0, C_Complex(Double(0.0)), TmpR, BlockSize, 0);
                CMatrixGEMM(RowsCount, N-BlockStart-BlockSize, BlockSize, C_Complex(Double(1.0)), TmpA, 0, 0, 0, TmpR, BlockSize, 0, 0, C_Complex(Double(1.0)), A, BlockStart, BlockStart+BlockSize);
            end
            else
            begin
                
                //
                // Level 2 algorithm
                //
                I:=0;
                while I<=BlockSize-1 do
                begin
                    i1_ := (I) - (1);
                    for i_ := 1 to RowsCount-I do
                    begin
                        T[i_] := TmpA[i_+i1_,I];
                    end;
                    T[1] := C_Complex(1);
                    ComplexApplyReflectionFromTheLeft(A, Conj(TauBuf[I]), T, BlockStart+I, M-1, BlockStart+BlockSize, N-1, WORK);
                    Inc(I);
                end;
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart+BlockSize;
    end;
end;


(*************************************************************************
LQ decomposition of a rectangular complex matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and L in compact form
    Tau -   array of scalar factors which are used to form matrix Q. Array
            whose indexes range within [0.. Min(M,N)-1]

Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
MxM, L - lower triangular (or lower trapezoid) matrix of size MxN.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure CMatrixLQ(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
var
    WORK : TComplex1DArray;
    T : TComplex1DArray;
    TauBuf : TComplex1DArray;
    MinMN : AlglibInteger;
    TmpA : TComplex2DArray;
    TmpT : TComplex2DArray;
    TmpR : TComplex2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    ColumnsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(Tau, MinMN);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, ABLASComplexBlockSize(A), N);
    SetLength(TmpT, ABLASComplexBlockSize(A), ABLASComplexBlockSize(A));
    SetLength(TmpR, M, 2*ABLASComplexBlockSize(A));
    
    //
    // Blocked code
    //
    BlockStart := 0;
    while BlockStart<>MinMN do
    begin
        
        //
        // Determine block size
        //
        BlockSize := MinMN-BlockStart;
        if BlockSize>ABLASComplexBlockSize(A) then
        begin
            BlockSize := ABLASComplexBlockSize(A);
        end;
        ColumnsCount := N-BlockStart;
        
        //
        // LQ decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        CMatrixCopy(BlockSize, ColumnsCount, A, BlockStart, BlockStart, TmpA, 0, 0);
        CMatrixLQBaseCase(TmpA, BlockSize, ColumnsCount, WORK, T, TauBuf);
        CMatrixCopy(BlockSize, ColumnsCount, TmpA, 0, 0, A, BlockStart, BlockStart);
        i1_ := (0) - (BlockStart);
        for i_ := BlockStart to BlockStart+BlockSize-1 do
        begin
            Tau[i_] := TauBuf[i_+i1_];
        end;
        
        //
        // Update the rest, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if BlockStart+BlockSize<=M-1 then
        begin
            if M-BlockStart-BlockSize>=2*ABLASComplexBlockSize(A) then
            begin
                
                //
                // Prepare block reflector
                //
                CMatrixBlockReflector(TmpA, TauBuf, False, ColumnsCount, BlockSize, TmpT, WORK);
                
                //
                // Multiply the rest of A by Q.
                //
                // Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
                //
                CMatrixGEMM(M-BlockStart-BlockSize, BlockSize, ColumnsCount, C_Complex(Double(1.0)), A, BlockStart+BlockSize, BlockStart, 0, TmpA, 0, 0, 2, C_Complex(Double(0.0)), TmpR, 0, 0);
                CMatrixGEMM(M-BlockStart-BlockSize, BlockSize, BlockSize, C_Complex(Double(1.0)), TmpR, 0, 0, 0, TmpT, 0, 0, 0, C_Complex(Double(0.0)), TmpR, 0, BlockSize);
                CMatrixGEMM(M-BlockStart-BlockSize, ColumnsCount, BlockSize, C_Complex(Double(1.0)), TmpR, 0, BlockSize, 0, TmpA, 0, 0, 0, C_Complex(Double(1.0)), A, BlockStart+BlockSize, BlockStart);
            end
            else
            begin
                
                //
                // Level 2 algorithm
                //
                I:=0;
                while I<=BlockSize-1 do
                begin
                    i1_ := (I) - (1);
                    for i_ := 1 to ColumnsCount-I do
                    begin
                        T[i_] := Conj(TmpA[I,i_+i1_]);
                    end;
                    T[1] := C_Complex(1);
                    ComplexApplyReflectionFromTheRight(A, TauBuf[I], T, BlockStart+BlockSize, M-1, BlockStart+I, N-1, WORK);
                    Inc(I);
                end;
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart+BlockSize;
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of RMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.
    Tau     -   scalar factors which are used to form Q.
                Output of the RMatrixQR subroutine.
    QColumns -  required number of columns of matrix Q. M>=QColumns>=0.

Output parameters:
    Q       -   first QColumns columns of matrix Q.
                Array whose indexes range within [0..M-1, 0..QColumns-1].
                If QColumns=0, the array remains unchanged.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixQRUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    TauBuf : TReal1DArray;
    MinMN : AlglibInteger;
    RefCnt : AlglibInteger;
    TmpA : TReal2DArray;
    TmpT : TReal2DArray;
    TmpR : TReal2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    RowsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(QColumns<=M, 'UnpackQFromQR: QColumns>M!');
    if (M<=0) or (N<=0) or (QColumns<=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    RefCnt := Min(MinMN, QColumns);
    SetLength(Q, M, QColumns);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=QColumns-1 do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(WORK, Max(M, QColumns)+1);
    SetLength(T, Max(M, QColumns)+1);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, M, ABLASBlockSize(A));
    SetLength(TmpT, ABLASBlockSize(A), 2*ABLASBlockSize(A));
    SetLength(TmpR, 2*ABLASBlockSize(A), QColumns);
    
    //
    // Blocked code
    //
    BlockStart := ABLASBlockSize(A)*(RefCnt div ABLASBlockSize(A));
    BlockSize := RefCnt-BlockStart;
    while BlockStart>=0 do
    begin
        RowsCount := M-BlockStart;
        
        //
        // Copy current block
        //
        RMatrixCopy(RowsCount, BlockSize, A, BlockStart, BlockStart, TmpA, 0, 0);
        APVMove(@TauBuf[0], 0, BlockSize-1, @Tau[0], BlockStart, BlockStart+BlockSize-1);
        
        //
        // Update, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if QColumns>=2*ABLASBlockSize(A) then
        begin
            
            //
            // Prepare block reflector
            //
            RMatrixBlockReflector(TmpA, TauBuf, True, RowsCount, BlockSize, TmpT, WORK);
            
            //
            // Multiply matrix by Q.
            //
            // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            //
            RMatrixGEMM(BlockSize, QColumns, RowsCount, Double(1.0), TmpA, 0, 0, 1, Q, BlockStart, 0, 0, Double(0.0), TmpR, 0, 0);
            RMatrixGEMM(BlockSize, QColumns, BlockSize, Double(1.0), TmpT, 0, 0, 0, TmpR, 0, 0, 0, Double(0.0), TmpR, BlockSize, 0);
            RMatrixGEMM(RowsCount, QColumns, BlockSize, Double(1.0), TmpA, 0, 0, 0, TmpR, BlockSize, 0, 0, Double(1.0), Q, BlockStart, 0);
        end
        else
        begin
            
            //
            // Level 2 algorithm
            //
            I:=BlockSize-1;
            while I>=0 do
            begin
                i1_ := (I) - (1);
                for i_ := 1 to RowsCount-I do
                begin
                    T[i_] := TmpA[i_+i1_,I];
                end;
                T[1] := 1;
                ApplyReflectionFromTheLeft(Q, TauBuf[I], T, BlockStart+I, M-1, 0, QColumns-1, WORK);
                Dec(I);
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart-ABLASBlockSize(A);
        BlockSize := ABLASBlockSize(A);
    end;
end;


(*************************************************************************
Unpacking of matrix R from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of RMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    R       -   matrix R, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixQRUnpackR(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TReal2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    K := Min(M, N);
    SetLength(R, M, N);
    I:=0;
    while I<=N-1 do
    begin
        R[0,I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        APVMove(@R[I][0], 0, N-1, @R[0][0], 0, N-1);
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        APVMove(@R[I][0], I, N-1, @A[I][0], I, N-1);
        Inc(I);
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices L and Q in compact form.
                Output of RMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.
    Tau     -   scalar factors which are used to form Q.
                Output of the RMatrixLQ subroutine.
    QRows   -   required number of rows in matrix Q. N>=QRows>=0.

Output parameters:
    Q       -   first QRows rows of matrix Q. Array whose indexes range
                within [0..QRows-1, 0..N-1]. If QRows=0, the array remains
                unchanged.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixLQUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QRows : AlglibInteger;
     var Q : TReal2DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    TauBuf : TReal1DArray;
    MinMN : AlglibInteger;
    RefCnt : AlglibInteger;
    TmpA : TReal2DArray;
    TmpT : TReal2DArray;
    TmpR : TReal2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    ColumnsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    Assert(QRows<=N, 'RMatrixLQUnpackQ: QRows>N!');
    if (M<=0) or (N<=0) or (QRows<=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    RefCnt := Min(MinMN, QRows);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, ABLASBlockSize(A), N);
    SetLength(TmpT, ABLASBlockSize(A), 2*ABLASBlockSize(A));
    SetLength(TmpR, QRows, 2*ABLASBlockSize(A));
    SetLength(Q, QRows, N);
    I:=0;
    while I<=QRows-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Blocked code
    //
    BlockStart := ABLASBlockSize(A)*(RefCnt div ABLASBlockSize(A));
    BlockSize := RefCnt-BlockStart;
    while BlockStart>=0 do
    begin
        ColumnsCount := N-BlockStart;
        
        //
        // Copy submatrix
        //
        RMatrixCopy(BlockSize, ColumnsCount, A, BlockStart, BlockStart, TmpA, 0, 0);
        APVMove(@TauBuf[0], 0, BlockSize-1, @Tau[0], BlockStart, BlockStart+BlockSize-1);
        
        //
        // Update matrix, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if QRows>=2*ABLASBlockSize(A) then
        begin
            
            //
            // Prepare block reflector
            //
            RMatrixBlockReflector(TmpA, TauBuf, False, ColumnsCount, BlockSize, TmpT, WORK);
            
            //
            // Multiply the rest of A by Q'.
            //
            // Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
            //
            RMatrixGEMM(QRows, BlockSize, ColumnsCount, Double(1.0), Q, 0, BlockStart, 0, TmpA, 0, 0, 1, Double(0.0), TmpR, 0, 0);
            RMatrixGEMM(QRows, BlockSize, BlockSize, Double(1.0), TmpR, 0, 0, 0, TmpT, 0, 0, 1, Double(0.0), TmpR, 0, BlockSize);
            RMatrixGEMM(QRows, ColumnsCount, BlockSize, Double(1.0), TmpR, 0, BlockSize, 0, TmpA, 0, 0, 0, Double(1.0), Q, 0, BlockStart);
        end
        else
        begin
            
            //
            // Level 2 algorithm
            //
            I:=BlockSize-1;
            while I>=0 do
            begin
                APVMove(@T[0], 1, ColumnsCount-I, @TmpA[I][0], I, ColumnsCount-1);
                T[1] := 1;
                ApplyReflectionFromTheRight(Q, TauBuf[I], T, 0, QRows-1, BlockStart+I, N-1, WORK);
                Dec(I);
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart-ABLASBlockSize(A);
        BlockSize := ABLASBlockSize(A);
    end;
end;


(*************************************************************************
Unpacking of matrix L from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices Q and L in compact form.
                Output of RMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    L       -   matrix L, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixLQUnpackL(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    SetLength(L, M, N);
    I:=0;
    while I<=N-1 do
    begin
        L[0,I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        APVMove(@L[I][0], 0, N-1, @L[0][0], 0, N-1);
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        K := Min(I, N-1);
        APVMove(@L[I][0], 0, K, @A[I][0], 0, K);
        Inc(I);
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from QR decomposition of a complex matrix A.

Input parameters:
    A           -   matrices Q and R in compact form.
                    Output of CMatrixQR subroutine .
    M           -   number of rows in matrix A. M>=0.
    N           -   number of columns in matrix A. N>=0.
    Tau         -   scalar factors which are used to form Q.
                    Output of CMatrixQR subroutine .
    QColumns    -   required number of columns in matrix Q. M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array whose index ranges within [0..M-1, 0..QColumns-1].
                    If QColumns=0, array isn't changed.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixQRUnpackQ(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QColumns : AlglibInteger;
     var Q : TComplex2DArray);
var
    WORK : TComplex1DArray;
    T : TComplex1DArray;
    TauBuf : TComplex1DArray;
    MinMN : AlglibInteger;
    RefCnt : AlglibInteger;
    TmpA : TComplex2DArray;
    TmpT : TComplex2DArray;
    TmpR : TComplex2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    RowsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(QColumns<=M, 'UnpackQFromQR: QColumns>M!');
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    RefCnt := Min(MinMN, QColumns);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, M, ABLASComplexBlockSize(A));
    SetLength(TmpT, ABLASComplexBlockSize(A), ABLASComplexBlockSize(A));
    SetLength(TmpR, 2*ABLASComplexBlockSize(A), QColumns);
    SetLength(Q, M, QColumns);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=QColumns-1 do
        begin
            if I=J then
            begin
                Q[I,J] := C_Complex(1);
            end
            else
            begin
                Q[I,J] := C_Complex(0);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Blocked code
    //
    BlockStart := ABLASComplexBlockSize(A)*(RefCnt div ABLASComplexBlockSize(A));
    BlockSize := RefCnt-BlockStart;
    while BlockStart>=0 do
    begin
        RowsCount := M-BlockStart;
        
        //
        // QR decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        CMatrixCopy(RowsCount, BlockSize, A, BlockStart, BlockStart, TmpA, 0, 0);
        i1_ := (BlockStart) - (0);
        for i_ := 0 to BlockSize-1 do
        begin
            TauBuf[i_] := Tau[i_+i1_];
        end;
        
        //
        // Update matrix, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if QColumns>=2*ABLASComplexBlockSize(A) then
        begin
            
            //
            // Prepare block reflector
            //
            CMatrixBlockReflector(TmpA, TauBuf, True, RowsCount, BlockSize, TmpT, WORK);
            
            //
            // Multiply the rest of A by Q.
            //
            // Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            //
            CMatrixGEMM(BlockSize, QColumns, RowsCount, C_Complex(Double(1.0)), TmpA, 0, 0, 2, Q, BlockStart, 0, 0, C_Complex(Double(0.0)), TmpR, 0, 0);
            CMatrixGEMM(BlockSize, QColumns, BlockSize, C_Complex(Double(1.0)), TmpT, 0, 0, 0, TmpR, 0, 0, 0, C_Complex(Double(0.0)), TmpR, BlockSize, 0);
            CMatrixGEMM(RowsCount, QColumns, BlockSize, C_Complex(Double(1.0)), TmpA, 0, 0, 0, TmpR, BlockSize, 0, 0, C_Complex(Double(1.0)), Q, BlockStart, 0);
        end
        else
        begin
            
            //
            // Level 2 algorithm
            //
            I:=BlockSize-1;
            while I>=0 do
            begin
                i1_ := (I) - (1);
                for i_ := 1 to RowsCount-I do
                begin
                    T[i_] := TmpA[i_+i1_,I];
                end;
                T[1] := C_Complex(1);
                ComplexApplyReflectionFromTheLeft(Q, TauBuf[I], T, BlockStart+I, M-1, 0, QColumns-1, WORK);
                Dec(I);
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart-ABLASComplexBlockSize(A);
        BlockSize := ABLASComplexBlockSize(A);
    end;
end;


(*************************************************************************
Unpacking of matrix R from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of CMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    R       -   matrix R, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixQRUnpackR(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TComplex2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    i_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    K := Min(M, N);
    SetLength(R, M, N);
    I:=0;
    while I<=N-1 do
    begin
        R[0,I] := C_Complex(0);
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        for i_ := 0 to N-1 do
        begin
            R[I,i_] := R[0,i_];
        end;
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        for i_ := I to N-1 do
        begin
            R[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from LQ decomposition of a complex matrix A.

Input parameters:
    A           -   matrices Q and R in compact form.
                    Output of CMatrixLQ subroutine .
    M           -   number of rows in matrix A. M>=0.
    N           -   number of columns in matrix A. N>=0.
    Tau         -   scalar factors which are used to form Q.
                    Output of CMatrixLQ subroutine .
    QRows       -   required number of rows in matrix Q. N>=QColumns>=0.

Output parameters:
    Q           -   first QRows rows of matrix Q.
                    Array whose index ranges within [0..QRows-1, 0..N-1].
                    If QRows=0, array isn't changed.

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixLQUnpackQ(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QRows : AlglibInteger;
     var Q : TComplex2DArray);
var
    WORK : TComplex1DArray;
    T : TComplex1DArray;
    TauBuf : TComplex1DArray;
    MinMN : AlglibInteger;
    RefCnt : AlglibInteger;
    TmpA : TComplex2DArray;
    TmpT : TComplex2DArray;
    TmpR : TComplex2DArray;
    BlockStart : AlglibInteger;
    BlockSize : AlglibInteger;
    ColumnsCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    
    //
    // Init
    //
    MinMN := Min(M, N);
    RefCnt := Min(MinMN, QRows);
    SetLength(WORK, Max(M, N)+1);
    SetLength(T, Max(M, N)+1);
    SetLength(TauBuf, MinMN);
    SetLength(TmpA, ABLASComplexBlockSize(A), N);
    SetLength(TmpT, ABLASComplexBlockSize(A), ABLASComplexBlockSize(A));
    SetLength(TmpR, QRows, 2*ABLASComplexBlockSize(A));
    SetLength(Q, QRows, N);
    I:=0;
    while I<=QRows-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                Q[I,J] := C_Complex(1);
            end
            else
            begin
                Q[I,J] := C_Complex(0);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Blocked code
    //
    BlockStart := ABLASComplexBlockSize(A)*(RefCnt div ABLASComplexBlockSize(A));
    BlockSize := RefCnt-BlockStart;
    while BlockStart>=0 do
    begin
        ColumnsCount := N-BlockStart;
        
        //
        // LQ decomposition of submatrix.
        // Matrix is copied to temporary storage to solve
        // some TLB issues arising from non-contiguous memory
        // access pattern.
        //
        CMatrixCopy(BlockSize, ColumnsCount, A, BlockStart, BlockStart, TmpA, 0, 0);
        i1_ := (BlockStart) - (0);
        for i_ := 0 to BlockSize-1 do
        begin
            TauBuf[i_] := Tau[i_+i1_];
        end;
        
        //
        // Update matrix, choose between:
        // a) Level 2 algorithm (when the rest of the matrix is small enough)
        // b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
        //    representation for products of Householder transformations',
        //    by R. Schreiber and C. Van Loan.
        //
        if QRows>=2*ABLASComplexBlockSize(A) then
        begin
            
            //
            // Prepare block reflector
            //
            CMatrixBlockReflector(TmpA, TauBuf, False, ColumnsCount, BlockSize, TmpT, WORK);
            
            //
            // Multiply the rest of A by Q'.
            //
            // Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
            //
            CMatrixGEMM(QRows, BlockSize, ColumnsCount, C_Complex(Double(1.0)), Q, 0, BlockStart, 0, TmpA, 0, 0, 2, C_Complex(Double(0.0)), TmpR, 0, 0);
            CMatrixGEMM(QRows, BlockSize, BlockSize, C_Complex(Double(1.0)), TmpR, 0, 0, 0, TmpT, 0, 0, 2, C_Complex(Double(0.0)), TmpR, 0, BlockSize);
            CMatrixGEMM(QRows, ColumnsCount, BlockSize, C_Complex(Double(1.0)), TmpR, 0, BlockSize, 0, TmpA, 0, 0, 0, C_Complex(Double(1.0)), Q, 0, BlockStart);
        end
        else
        begin
            
            //
            // Level 2 algorithm
            //
            I:=BlockSize-1;
            while I>=0 do
            begin
                i1_ := (I) - (1);
                for i_ := 1 to ColumnsCount-I do
                begin
                    T[i_] := Conj(TmpA[I,i_+i1_]);
                end;
                T[1] := C_Complex(1);
                ComplexApplyReflectionFromTheRight(Q, Conj(TauBuf[I]), T, 0, QRows-1, BlockStart+I, N-1, WORK);
                Dec(I);
            end;
        end;
        
        //
        // Advance
        //
        BlockStart := BlockStart-ABLASComplexBlockSize(A);
        BlockSize := ABLASComplexBlockSize(A);
    end;
end;


(*************************************************************************
Unpacking of matrix L from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices Q and L in compact form.
                Output of CMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    L       -   matrix L, array[0..M-1, 0..N-1].

  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixLQUnpackL(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TComplex2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    i_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    SetLength(L, M, N);
    I:=0;
    while I<=N-1 do
    begin
        L[0,I] := C_Complex(0);
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        for i_ := 0 to N-1 do
        begin
            L[I,i_] := L[0,i_];
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        K := Min(I, N-1);
        for i_ := 0 to K do
        begin
            L[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Reduction of a rectangular matrix to  bidiagonal form

The algorithm reduces the rectangular matrix A to  bidiagonal form by
orthogonal transformations P and Q: A = Q*B*P.

Input parameters:
    A       -   source matrix. array[0..M-1, 0..N-1]
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.

Output parameters:
    A       -   matrices Q, B, P in compact form (see below).
    TauQ    -   scalar factors which are used to form matrix Q.
    TauP    -   scalar factors which are used to form matrix P.

The main diagonal and one of the  secondary  diagonals  of  matrix  A  are
replaced with bidiagonal  matrix  B.  Other  elements  contain  elementary
reflections which form MxM matrix Q and NxN matrix P, respectively.

If M>=N, B is the upper  bidiagonal  MxN  matrix  and  is  stored  in  the
corresponding  elements  of  matrix  A.  Matrix  Q  is  represented  as  a
product   of   elementary   reflections   Q = H(0)*H(1)*...*H(n-1),  where
H(i) = 1-tau*v*v'. Here tau is a scalar which is stored  in  TauQ[i],  and
vector v has the following  structure:  v(0:i-1)=0, v(i)=1, v(i+1:m-1)  is
stored   in   elements   A(i+1:m-1,i).   Matrix   P  is  as  follows:  P =
G(0)*G(1)*...*G(n-2), where G(i) = 1 - tau*u*u'. Tau is stored in TauP[i],
u(0:i)=0, u(i+1)=1, u(i+2:n-1) is stored in elements A(i,i+2:n-1).

If M<N, B is the  lower  bidiagonal  MxN  matrix  and  is  stored  in  the
corresponding   elements  of  matrix  A.  Q = H(0)*H(1)*...*H(m-2),  where
H(i) = 1 - tau*v*v', tau is stored in TauQ, v(0:i)=0, v(i+1)=1, v(i+2:m-1)
is    stored    in   elements   A(i+2:m-1,i).    P = G(0)*G(1)*...*G(m-1),
G(i) = 1-tau*u*u', tau is stored in  TauP,  u(0:i-1)=0, u(i)=1, u(i+1:n-1)
is stored in A(i,i+1:n-1).

EXAMPLE:

m=6, n=5 (m > n):               m=5, n=6 (m < n):

(  d   e   u1  u1  u1 )         (  d   u1  u1  u1  u1  u1 )
(  v1  d   e   u2  u2 )         (  e   d   u2  u2  u2  u2 )
(  v1  v2  d   e   u3 )         (  v1  e   d   u3  u3  u3 )
(  v1  v2  v3  d   e  )         (  v1  v2  e   d   u4  u4 )
(  v1  v2  v3  v4  d  )         (  v1  v2  v3  e   d   u5 )
(  v1  v2  v3  v4  v5 )

Here vi and ui are vectors which form H(i) and G(i), and d and e -
are the diagonal and off-diagonal elements of matrix B.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************)
procedure RMatrixBD(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var TauQ : TReal1DArray;
     var TauP : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    MinMN : AlglibInteger;
    MaxMN : AlglibInteger;
    I : AlglibInteger;
    LTau : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Prepare
    //
    if (N<=0) or (M<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    MaxMN := Max(M, N);
    SetLength(Work, MaxMN+1);
    SetLength(T, MaxMN+1);
    if M>=N then
    begin
        SetLength(TauQ, N);
        SetLength(TauP, N);
    end
    else
    begin
        SetLength(TauQ, M);
        SetLength(TauP, M);
    end;
    if M>=N then
    begin
        
        //
        // Reduce to upper bidiagonal form
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Generate elementary reflector H(i) to annihilate A(i+1:m-1,i)
            //
            i1_ := (I) - (1);
            for i_ := 1 to M-I do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            GenerateReflection(T, M-I, LTau);
            TauQ[I] := LTau;
            i1_ := (1) - (I);
            for i_ := I to M-1 do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            T[1] := 1;
            
            //
            // Apply H(i) to A(i:m-1,i+1:n-1) from the left
            //
            ApplyReflectionFromTheLeft(A, LTau, T, I, M-1, I+1, N-1, WORK);
            if I<N-1 then
            begin
                
                //
                // Generate elementary reflector G(i) to annihilate
                // A(i,i+2:n-1)
                //
                APVMove(@T[0], 1, N-I-1, @A[I][0], I+1, N-1);
                GenerateReflection(T, N-1-I, LTau);
                TauP[I] := LTau;
                APVMove(@A[I][0], I+1, N-1, @T[0], 1, N-1-I);
                T[1] := 1;
                
                //
                // Apply G(i) to A(i+1:m-1,i+1:n-1) from the right
                //
                ApplyReflectionFromTheRight(A, LTau, T, I+1, M-1, I+1, N-1, WORK);
            end
            else
            begin
                TAUP[I] := 0;
            end;
            Inc(I);
        end;
    end
    else
    begin
        
        //
        // Reduce to lower bidiagonal form
        //
        I:=0;
        while I<=M-1 do
        begin
            
            //
            // Generate elementary reflector G(i) to annihilate A(i,i+1:n-1)
            //
            APVMove(@T[0], 1, N-I, @A[I][0], I, N-1);
            GenerateReflection(T, N-I, LTau);
            TauP[I] := LTau;
            APVMove(@A[I][0], I, N-1, @T[0], 1, N-I);
            T[1] := 1;
            
            //
            // Apply G(i) to A(i+1:m-1,i:n-1) from the right
            //
            ApplyReflectionFromTheRight(A, LTau, T, I+1, M-1, I, N-1, WORK);
            if I<M-1 then
            begin
                
                //
                // Generate elementary reflector H(i) to annihilate
                // A(i+2:m-1,i)
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to M-1-I do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                GenerateReflection(T, M-1-I, LTau);
                TauQ[I] := LTau;
                i1_ := (1) - (I+1);
                for i_ := I+1 to M-1 do
                begin
                    A[i_,I] := T[i_+i1_];
                end;
                T[1] := 1;
                
                //
                // Apply H(i) to A(i+1:m-1,i+1:n-1) from the left
                //
                ApplyReflectionFromTheLeft(A, LTau, T, I+1, M-1, I+1, N-1, WORK);
            end
            else
            begin
                TAUQ[I] := 0;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces a matrix to bidiagonal form.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    QColumns    -   required number of columns in matrix Q.
                    M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array[0..M-1, 0..QColumns-1]
                    If QColumns=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDUnpackQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(QColumns<=M, 'RMatrixBDUnpackQ: QColumns>M!');
    Assert(QColumns>=0, 'RMatrixBDUnpackQ: QColumns<0!');
    if (M=0) or (N=0) or (QColumns=0) then
    begin
        Exit;
    end;
    
    //
    // prepare Q
    //
    SetLength(Q, M, QColumns);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=QColumns-1 do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Calculate
    //
    RMatrixBDMultiplyByQ(QP, M, N, TauQ, Q, M, QColumns, False, False);
end;


(*************************************************************************
Multiplication by matrix Q which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by Q or Q'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    Z           -   multiplied matrix.
                    array[0..ZRows-1,0..ZColumns-1]
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=M, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=M, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by Q or Q'.

Output parameters:
    Z           -   product of Z and Q.
                    Array[0..ZRows-1,0..ZColumns-1]
                    If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDMultiplyByQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
var
    I : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IStep : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    Mx : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) or (ZRows<=0) or (ZColumns<=0) then
    begin
        Exit;
    end;
    Assert(FromTheRight and (ZColumns=M) or  not FromTheRight and (ZRows=M), 'RMatrixBDMultiplyByQ: incorrect Z size!');
    
    //
    // init
    //
    Mx := Max(M, N);
    Mx := Max(Mx, ZRows);
    Mx := Max(Mx, ZColumns);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    if M>=N then
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := 0;
            I2 := N-1;
            IStep := +1;
        end
        else
        begin
            I1 := N-1;
            I2 := 0;
            IStep := -1;
        end;
        if DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        I := I1;
        repeat
            i1_ := (I) - (1);
            for i_ := 1 to M-I do
            begin
                V[i_] := QP[i_+i1_,I];
            end;
            V[1] := 1;
            if FromTheRight then
            begin
                ApplyReflectionFromTheRight(Z, TAUQ[I], V, 0, ZRows-1, I, M-1, WORK);
            end
            else
            begin
                ApplyReflectionFromTheLeft(Z, TAUQ[I], V, I, M-1, 0, ZColumns-1, WORK);
            end;
            I := I+IStep;
        until I=I2+IStep;
    end
    else
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := 0;
            I2 := M-2;
            IStep := +1;
        end
        else
        begin
            I1 := M-2;
            I2 := 0;
            IStep := -1;
        end;
        if DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        if M-1>0 then
        begin
            I := I1;
            repeat
                i1_ := (I+1) - (1);
                for i_ := 1 to M-I-1 do
                begin
                    V[i_] := QP[i_+i1_,I];
                end;
                V[1] := 1;
                if FromTheRight then
                begin
                    ApplyReflectionFromTheRight(Z, TAUQ[I], V, 0, ZRows-1, I+1, M-1, WORK);
                end
                else
                begin
                    ApplyReflectionFromTheLeft(Z, TAUQ[I], V, I+1, M-1, 0, ZColumns-1, WORK);
                end;
                I := I+IStep;
            until I=I2+IStep;
        end;
    end;
end;


(*************************************************************************
Unpacking matrix P which reduces matrix A to bidiagonal form.
The subroutine returns transposed matrix P.

Input parameters:
    QP      -   matrices Q and P in compact form.
                Output of ToBidiagonal subroutine.
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.
    TAUP    -   scalar factors which are used to form P.
                Output of ToBidiagonal subroutine.
    PTRows  -   required number of rows of matrix P^T. N >= PTRows >= 0.

Output parameters:
    PT      -   first PTRows columns of matrix P^T
                Array[0..PTRows-1, 0..N-1]
                If PTRows=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDUnpackPT(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     PTRows : AlglibInteger;
     var PT : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(PTRows<=N, 'RMatrixBDUnpackPT: PTRows>N!');
    Assert(PTRows>=0, 'RMatrixBDUnpackPT: PTRows<0!');
    if (M=0) or (N=0) or (PTRows=0) then
    begin
        Exit;
    end;
    
    //
    // prepare PT
    //
    SetLength(PT, PTRows, N);
    I:=0;
    while I<=PTRows-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                PT[I,J] := 1;
            end
            else
            begin
                PT[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Calculate
    //
    RMatrixBDMultiplyByP(QP, M, N, TauP, PT, PTRows, N, True, True);
end;


(*************************************************************************
Multiplication by matrix P which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by P or P'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of RMatrixBD subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUP        -   scalar factors which are used to form P.
                    Output of RMatrixBD subroutine.
    Z           -   multiplied matrix.
                    Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=N, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=N, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by P or P'.

Output parameters:
    Z - product of Z and P.
                Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
                If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDMultiplyByP(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
var
    I : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    Mx : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IStep : AlglibInteger;
begin
    if (M<=0) or (N<=0) or (ZRows<=0) or (ZColumns<=0) then
    begin
        Exit;
    end;
    Assert(FromTheRight and (ZColumns=N) or  not FromTheRight and (ZRows=N), 'RMatrixBDMultiplyByP: incorrect Z size!');
    
    //
    // init
    //
    Mx := Max(M, N);
    Mx := Max(Mx, ZRows);
    Mx := Max(Mx, ZColumns);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    if M>=N then
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := N-2;
            I2 := 0;
            IStep := -1;
        end
        else
        begin
            I1 := 0;
            I2 := N-2;
            IStep := +1;
        end;
        if  not DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        if N-1>0 then
        begin
            I := I1;
            repeat
                APVMove(@V[0], 1, N-1-I, @QP[I][0], I+1, N-1);
                V[1] := 1;
                if FromTheRight then
                begin
                    ApplyReflectionFromTheRight(Z, TAUP[I], V, 0, ZRows-1, I+1, N-1, WORK);
                end
                else
                begin
                    ApplyReflectionFromTheLeft(Z, TAUP[I], V, I+1, N-1, 0, ZColumns-1, WORK);
                end;
                I := I+IStep;
            until I=I2+IStep;
        end;
    end
    else
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := M-1;
            I2 := 0;
            IStep := -1;
        end
        else
        begin
            I1 := 0;
            I2 := M-1;
            IStep := +1;
        end;
        if  not DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        I := I1;
        repeat
            APVMove(@V[0], 1, N-I, @QP[I][0], I, N-1);
            V[1] := 1;
            if FromTheRight then
            begin
                ApplyReflectionFromTheRight(Z, TauP[I], V, 0, ZRows-1, I, N-1, WORK);
            end
            else
            begin
                ApplyReflectionFromTheLeft(Z, TauP[I], V, I, N-1, 0, ZColumns-1, WORK);
            end;
            I := I+IStep;
        until I=I2+IStep;
    end;
end;


(*************************************************************************
Unpacking of the main and secondary diagonals of bidiagonal decomposition
of matrix A.

Input parameters:
    B   -   output of RMatrixBD subroutine.
    M   -   number of rows in matrix B.
    N   -   number of columns in matrix B.

Output parameters:
    IsUpper -   True, if the matrix is upper bidiagonal.
                otherwise IsUpper is False.
    D       -   the main diagonal.
                Array whose index ranges within [0..Min(M,N)-1].
    E       -   the secondary diagonal (upper or lower, depending on
                the value of IsUpper).
                Array index ranges within [0..Min(M,N)-1], the last
                element is not used.

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDUnpackDiagonals(const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var IsUpper : Boolean;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
begin
    IsUpper := M>=N;
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        SetLength(D, N);
        SetLength(E, N);
        I:=0;
        while I<=N-2 do
        begin
            D[I] := B[I,I];
            E[I] := B[I,I+1];
            Inc(I);
        end;
        D[N-1] := B[N-1,N-1];
    end
    else
    begin
        SetLength(D, M);
        SetLength(E, M);
        I:=0;
        while I<=M-2 do
        begin
            D[I] := B[I,I];
            E[I] := B[I+1,I];
            Inc(I);
        end;
        D[M-1] := B[M-1,M-1];
    end;
end;


(*************************************************************************
Reduction of a square matrix to  upper Hessenberg form: Q'*A*Q = H,
where Q is an orthogonal matrix, H - Hessenberg matrix.

Input parameters:
    A       -   matrix A with elements [0..N-1, 0..N-1]
    N       -   size of matrix A.

Output parameters:
    A       -   matrices Q and P in  compact form (see below).
    Tau     -   array of scalar factors which are used to form matrix Q.
                Array whose index ranges within [0..N-2]

Matrix H is located on the main diagonal, on the lower secondary  diagonal
and above the main diagonal of matrix A. The elements which are used to
form matrix Q are situated in array Tau and below the lower secondary
diagonal of matrix A as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(0)*H(2)*...*H(n-2),

where each H(i) is given by

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - is a real vector,
so that v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) stored in A(i+2:n-1,i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure RMatrixHessenberg(var A : TReal2DArray;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    I : AlglibInteger;
    V : Double;
    T : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=0, 'RMatrixHessenberg: incorrect N!');
    
    //
    // Quick return if possible
    //
    if N<=1 then
    begin
        Exit;
    end;
    SetLength(Tau, N-2+1);
    SetLength(T, N+1);
    SetLength(WORK, N-1+1);
    I:=0;
    while I<=N-2 do
    begin
        
        //
        // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
        //
        i1_ := (I+1) - (1);
        for i_ := 1 to N-I-1 do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        GenerateReflection(T, N-I-1, V);
        i1_ := (1) - (I+1);
        for i_ := I+1 to N-1 do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        Tau[I] := V;
        T[1] := 1;
        
        //
        // Apply H(i) to A(1:ihi,i+1:ihi) from the right
        //
        ApplyReflectionFromTheRight(A, V, T, 0, N-1, I+1, N-1, WORK);
        
        //
        // Apply H(i) to A(i+1:ihi,i+1:n) from the left
        //
        ApplyReflectionFromTheLeft(A, V, T, I+1, N-1, I+1, N-1, WORK);
        Inc(I);
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces matrix A to upper Hessenberg form

Input parameters:
    A   -   output of RMatrixHessenberg subroutine.
    N   -   size of matrix A.
    Tau -   scalar factors which are used to form Q.
            Output of RMatrixHessenberg subroutine.

Output parameters:
    Q   -   matrix Q.
            Array whose indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixHessenbergUnpackQ(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(Q, N-1+1, N-1+1);
    SetLength(V, N-1+1);
    SetLength(WORK, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // unpack Q
    //
    I:=0;
    while I<=N-2 do
    begin
        
        //
        // Apply H(i)
        //
        i1_ := (I+1) - (1);
        for i_ := 1 to N-I-1 do
        begin
            V[i_] := A[i_+i1_,I];
        end;
        V[1] := 1;
        ApplyReflectionFromTheRight(Q, Tau[I], V, 0, N-1, I+1, N-1, WORK);
        Inc(I);
    end;
end;


(*************************************************************************
Unpacking matrix H (the result of matrix A reduction to upper Hessenberg form)

Input parameters:
    A   -   output of RMatrixHessenberg subroutine.
    N   -   size of matrix A.

Output parameters:
    H   -   matrix H. Array whose indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     2005-2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixHessenbergUnpackH(const A : TReal2DArray;
     N : AlglibInteger;
     var H : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
begin
    if N=0 then
    begin
        Exit;
    end;
    SetLength(H, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I-2 do
        begin
            H[I,J] := 0;
            Inc(J);
        end;
        J := Max(0, I-1);
        APVMove(@H[I][0], J, N-1, @A[I][0], J, N-1);
        Inc(I);
    end;
end;


(*************************************************************************
Reduction of a symmetric matrix which is given by its higher or lower
triangular part to a tridiagonal matrix using orthogonal similarity
transformation: Q'*A*Q=T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is given
                by its upper triangle, and the lower triangle is not used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  d   e   v1  v2  v3 )              (  d                  )
    (      d   e   v2  v3 )              (  e   d              )
    (          d   e   v3 )              (  v0  e   d          )
    (              d   e  )              (  v0  v1  e   d      )
    (                  d  )              (  v0  v1  v2  e   d  )

  where d and e denote diagonal and off-diagonal elements of T, and vi
  denotes an element of the vector defining H(i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure SMatrixTD(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TReal1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
    ALPHA : Double;
    TAUI : Double;
    V : Double;
    T : TReal1DArray;
    T2 : TReal1DArray;
    T3 : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(T, N+1);
    SetLength(T2, N+1);
    SetLength(T3, N+1);
    if N>1 then
    begin
        SetLength(Tau, N-2+1);
    end;
    SetLength(D, N-1+1);
    if N>1 then
    begin
        SetLength(E, N-2+1);
    end;
    if IsUpper then
    begin
        
        //
        // Reduce the upper triangle of A
        //
        I:=N-2;
        while I>=0 do
        begin
            
            //
            // Generate elementary reflector H() = E - tau * v * v'
            //
            if I>=1 then
            begin
                i1_ := (0) - (2);
                for i_ := 2 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
            end;
            T[1] := A[I,I+1];
            GenerateReflection(T, I+1, TauI);
            if I>=1 then
            begin
                i1_ := (2) - (0);
                for i_ := 0 to I-1 do
                begin
                    A[i_,I+1] := T[i_+i1_];
                end;
            end;
            A[I,I+1] := T[1];
            E[I] := A[I,I+1];
            if AP_FP_Neq(TAUI,0) then
            begin
                
                //
                // Apply H from both sides to A
                //
                A[I,I+1] := 1;
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                SymmetricMatrixVectorMultiply(A, IsUpper, 0, I, T, TauI, T3);
                APVMove(@Tau[0], 0, I, @T3[0], 1, I+1);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                V := 0.0;
                for i_ := 0 to I do
                begin
                    V := V + Tau[i_]*A[i_,I+1];
                end;
                ALPHA := -Double(0.5)*TAUI*V;
                for i_ := 0 to I do
                begin
                    Tau[i_] := Tau[i_] + Alpha*A[i_,I+1];
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                APVMove(@T3[0], 1, I+1, @Tau[0], 0, I);
                SymmetricRank2Update(A, IsUpper, 0, I, T, T3, T2, -1);
                A[I,I+1] := E[I];
            end;
            D[I+1] := A[I+1,I+1];
            TAU[I] := TAUI;
            Dec(I);
        end;
        D[0] := A[0,0];
    end
    else
    begin
        
        //
        // Reduce the lower triangle of A
        //
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Generate elementary reflector H = E - tau * v * v'
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I-1 do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            GenerateReflection(T, N-I-1, TauI);
            i1_ := (1) - (I+1);
            for i_ := I+1 to N-1 do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            E[I] := A[I+1,I];
            if AP_FP_Neq(TAUI,0) then
            begin
                
                //
                // Apply H from both sides to A
                //
                A[I+1,I] := 1;
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                SymmetricMatrixVectorMultiply(A, IsUpper, I+1, N-1, T, TauI, T2);
                APVMove(@Tau[0], I, N-2, @T2[0], 1, N-I-1);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                i1_ := (I+1)-(I);
                V := 0.0;
                for i_ := I to N-2 do
                begin
                    V := V + Tau[i_]*A[i_+i1_,I];
                end;
                ALPHA := -Double(0.5)*TAUI*V;
                i1_ := (I+1) - (I);
                for i_ := I to N-2 do
                begin
                    Tau[i_] := Tau[i_] + Alpha*A[i_+i1_,I];
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //     A := A - v * w' - w * v'
                //
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                APVMove(@T2[0], 1, N-I-1, @Tau[0], I, N-2);
                SymmetricRank2Update(A, IsUpper, I+1, N-1, T, T2, T3, -1);
                A[I+1,I] := E[I];
            end;
            D[I] := A[I,I];
            TAU[I] := TAUI;
            Inc(I);
        end;
        D[N-1] := A[N-1,N-1];
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces symmetric matrix to a tridiagonal
form.

Input parameters:
    A       -   the result of a SMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of SMatrixTD subroutine)
    Tau     -   the result of a SMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005-2010 by Bochkanov Sergey
*************************************************************************)
procedure SMatrixTDUnpackQ(const A : TReal2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(Q, N-1+1, N-1+1);
    SetLength(V, N+1);
    SetLength(WORK, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // unpack Q
    //
    if IsUpper then
    begin
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Apply H(i)
            //
            i1_ := (0) - (1);
            for i_ := 1 to I+1 do
            begin
                V[i_] := A[i_+i1_,I+1];
            end;
            V[I+1] := 1;
            ApplyReflectionFromTheLeft(Q, Tau[I], V, 0, I, 0, N-1, WORK);
            Inc(I);
        end;
    end
    else
    begin
        I:=N-2;
        while I>=0 do
        begin
            
            //
            // Apply H(i)
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I-1 do
            begin
                V[i_] := A[i_+i1_,I];
            end;
            V[1] := 1;
            ApplyReflectionFromTheLeft(Q, Tau[I], V, I+1, N-1, 0, N-1, WORK);
            Dec(I);
        end;
    end;
end;


(*************************************************************************
Reduction of a Hermitian matrix which is given  by  its  higher  or  lower
triangular part to a real  tridiagonal  matrix  using  unitary  similarity
transformation: Q'*A*Q = T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is  given
                by its upper triangle, and the lower triangle is not  used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of real symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of real symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a complex scalar, and v is a complex vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a complex scalar, and v is a complex vector with
  v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  d   e   v1  v2  v3 )              (  d                  )
    (      d   e   v2  v3 )              (  e   d              )
    (          d   e   v3 )              (  v0  e   d          )
    (              d   e  )              (  v0  v1  e   d      )
    (                  d  )              (  v0  v1  v2  e   d  )

where d and e denote diagonal and off-diagonal elements of T, and vi
denotes an element of the vector defining H(i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure HMatrixTD(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TComplex1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
    Alpha : Complex;
    TauI : Complex;
    V : Complex;
    T : TComplex1DArray;
    T2 : TComplex1DArray;
    T3 : TComplex1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        Assert(AP_FP_Eq(A[I,I].Y,0));
        Inc(I);
    end;
    if N>1 then
    begin
        SetLength(Tau, N-2+1);
        SetLength(E, N-2+1);
    end;
    SetLength(D, N-1+1);
    SetLength(T, N-1+1);
    SetLength(T2, N-1+1);
    SetLength(T3, N-1+1);
    if IsUpper then
    begin
        
        //
        // Reduce the upper triangle of A
        //
        A[N-1,N-1] := C_Complex(A[N-1,N-1].X);
        I:=N-2;
        while I>=0 do
        begin
            
            //
            // Generate elementary reflector H = I+1 - tau * v * v'
            //
            ALPHA := A[I,I+1];
            T[1] := ALPHA;
            if I>=1 then
            begin
                i1_ := (0) - (2);
                for i_ := 2 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
            end;
            ComplexGenerateReflection(T, I+1, TauI);
            if I>=1 then
            begin
                i1_ := (2) - (0);
                for i_ := 0 to I-1 do
                begin
                    A[i_,I+1] := T[i_+i1_];
                end;
            end;
            Alpha := T[1];
            E[I] := ALPHA.X;
            if C_NotEqualR(TAUI,0) then
            begin
                
                //
                // Apply H(I+1) from both sides to A
                //
                A[I,I+1] := C_Complex(1);
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                HermitianMatrixVectorMultiply(A, IsUpper, 0, I, T, TauI, T2);
                i1_ := (1) - (0);
                for i_ := 0 to I do
                begin
                    Tau[i_] := T2[i_+i1_];
                end;
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                V := C_Complex(0.0);
                for i_ := 0 to I do
                begin
                    V := C_Add(V,C_Mul(Conj(Tau[i_]),A[i_,I+1]));
                end;
                ALPHA := C_Opposite(C_Mul(C_MulR(TauI,Double(0.5)),V));
                for i_ := 0 to I do
                begin
                    Tau[i_] := C_Add(Tau[i_], C_Mul(Alpha, A[i_,I+1]));
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T3[i_] := Tau[i_+i1_];
                end;
                HermitianRank2Update(A, IsUpper, 0, I, T, T3, T2, C_Complex(-1));
            end
            else
            begin
                A[I,I] := C_Complex(A[I,I].X);
            end;
            A[I,I+1] := C_Complex(E[I]);
            D[I+1] := A[I+1,I+1].X;
            TAU[I] := TAUI;
            Dec(I);
        end;
        D[0] := A[0,0].X;
    end
    else
    begin
        
        //
        // Reduce the lower triangle of A
        //
        A[0,0] := C_Complex(A[0,0].X);
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Generate elementary reflector H = I - tau * v * v'
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I-1 do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            ComplexGenerateReflection(T, N-I-1, TauI);
            i1_ := (1) - (I+1);
            for i_ := I+1 to N-1 do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            E[I] := A[I+1,I].X;
            if C_NotEqualR(TauI,0) then
            begin
                
                //
                // Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                A[I+1,I] := C_Complex(1);
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                HermitianMatrixVectorMultiply(A, IsUpper, I+1, N-1, T, TauI, T2);
                i1_ := (1) - (I);
                for i_ := I to N-2 do
                begin
                    Tau[i_] := T2[i_+i1_];
                end;
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                i1_ := (I+1)-(I);
                V := C_Complex(0.0);
                for i_ := I to N-2 do
                begin
                    V := C_Add(V,C_Mul(Conj(Tau[i_]),A[i_+i1_,I]));
                end;
                ALPHA := C_Opposite(C_Mul(C_MulR(TauI,Double(0.5)),V));
                i1_ := (I+1) - (I);
                for i_ := I to N-2 do
                begin
                    Tau[i_] := C_Add(Tau[i_], C_Mul(Alpha, A[i_+i1_,I]));
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                // A := A - v * w' - w * v'
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                i1_ := (I) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T2[i_] := Tau[i_+i1_];
                end;
                HermitianRank2Update(A, IsUpper, I+1, N-1, T, T2, T3, C_Complex(-1));
            end
            else
            begin
                A[I+1,I+1] := C_Complex(A[I+1,I+1].X);
            end;
            A[I+1,I] := C_Complex(E[I]);
            D[I] := A[I,I].X;
            TAU[I] := TauI;
            Inc(I);
        end;
        D[N-1] := A[N-1,N-1].X;
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces a Hermitian matrix to a real  tridiagonal
form.

Input parameters:
    A       -   the result of a HMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of HMatrixTD subroutine)
    Tau     -   the result of a HMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005-2010 by Bochkanov Sergey
*************************************************************************)
procedure HMatrixTDUnpackQ(const A : TComplex2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TComplex1DArray;
     var Q : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TComplex1DArray;
    WORK : TComplex1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(Q, N-1+1, N-1+1);
    SetLength(V, N+1);
    SetLength(WORK, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                Q[I,J] := C_Complex(1);
            end
            else
            begin
                Q[I,J] := C_Complex(0);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // unpack Q
    //
    if IsUpper then
    begin
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Apply H(i)
            //
            i1_ := (0) - (1);
            for i_ := 1 to I+1 do
            begin
                V[i_] := A[i_+i1_,I+1];
            end;
            V[I+1] := C_Complex(1);
            ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, 0, I, 0, N-1, WORK);
            Inc(I);
        end;
    end
    else
    begin
        I:=N-2;
        while I>=0 do
        begin
            
            //
            // Apply H(i)
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I-1 do
            begin
                V[i_] := A[i_+i1_,I];
            end;
            V[1] := C_Complex(1);
            ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, I+1, N-1, 0, N-1, WORK);
            Dec(I);
        end;
    end;
end;


(*************************************************************************
Base case for real QR

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************)
procedure RMatrixQRBaseCase(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TReal1DArray;
     var T : TReal1DArray;
     var Tau : TReal1DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    
    //
    // Test the input arguments
    //
    K := MinMN;
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        i1_ := (I) - (1);
        for i_ := 1 to M-I do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        GenerateReflection(T, M-I, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to M-1 do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        T[1] := 1;
        if I<N then
        begin
            
            //
            // Apply H(i) to A(i:m-1,i+1:n-1) from the left
            //
            ApplyReflectionFromTheLeft(A, Tau[I], T, I, M-1, I+1, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Base case for real LQ

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************)
procedure RMatrixLQBaseCase(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TReal1DArray;
     var T : TReal1DArray;
     var Tau : TReal1DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Double;
begin
    MinMN := Min(M, N);
    K := Min(M, N);
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i,i+1:n-1)
        //
        APVMove(@T[0], 1, N-I, @A[I][0], I, N-1);
        GenerateReflection(T, N-I, Tmp);
        Tau[I] := Tmp;
        APVMove(@A[I][0], I, N-1, @T[0], 1, N-I);
        T[1] := 1;
        if I<N then
        begin
            
            //
            // Apply H(i) to A(i+1:m,i:n) from the right
            //
            ApplyReflectionFromTheRight(A, Tau[I], T, I+1, M-1, I, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Base case for complex QR

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************)
procedure CMatrixQRBaseCase(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TComplex1DArray;
     var T : TComplex1DArray;
     var Tau : TComplex1DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    MMI : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    if MinMN<=0 then
    begin
        Exit;
    end;
    
    //
    // Test the input arguments
    //
    K := Min(M, N);
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        MMI := M-I;
        i1_ := (I) - (1);
        for i_ := 1 to MMI do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        ComplexGenerateReflection(T, MMI, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to M-1 do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        T[1] := C_Complex(1);
        if I<N-1 then
        begin
            
            //
            // Apply H'(i) to A(i:m,i+1:n) from the left
            //
            ComplexApplyReflectionFromTheLeft(A, Conj(Tau[I]), T, I, M-1, I+1, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Base case for complex LQ

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994.
     Sergey Bochkanov, ALGLIB project, translation from FORTRAN to
     pseudocode, 2007-2010.
*************************************************************************)
procedure CMatrixLQBaseCase(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var WORK : TComplex1DArray;
     var T : TComplex1DArray;
     var Tau : TComplex1DArray);
var
    I : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    if MinMN<=0 then
    begin
        Exit;
    end;
    
    //
    // Test the input arguments
    //
    I:=0;
    while I<=MinMN-1 do
    begin
        
        //
        // Generate elementary reflector H(i)
        //
        // NOTE: ComplexGenerateReflection() generates left reflector,
        // i.e. H which reduces x by applyiong from the left, but we
        // need RIGHT reflector. So we replace H=E-tau*v*v' by H^H,
        // which changes v to conj(v).
        //
        i1_ := (I) - (1);
        for i_ := 1 to N-I do
        begin
            T[i_] := Conj(A[I,i_+i1_]);
        end;
        ComplexGenerateReflection(T, N-I, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to N-1 do
        begin
            A[I,i_] := Conj(T[i_+i1_]);
        end;
        T[1] := C_Complex(1);
        if I<M-1 then
        begin
            
            //
            // Apply H'(i)
            //
            ComplexApplyReflectionFromTheRight(A, Tau[I], T, I+1, M-1, I, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Generate block reflector:
* fill unused parts of reflectors matrix by zeros
* fill diagonal of reflectors matrix by ones
* generate triangular factor T

PARAMETERS:
    A           -   either LengthA*BlockSize (if ColumnwiseA) or
                    BlockSize*LengthA (if not ColumnwiseA) matrix of
                    elementary reflectors.
                    Modified on exit.
    Tau         -   scalar factors
    ColumnwiseA -   reflectors are stored in rows or in columns
    LengthA     -   length of largest reflector
    BlockSize   -   number of reflectors
    T           -   array[BlockSize,2*BlockSize]. Left BlockSize*BlockSize
                    submatrix stores triangular factor on exit.
    WORK        -   array[BlockSize]
    
  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixBlockReflector(var A : TReal2DArray;
     var Tau : TReal1DArray;
     ColumnwiseA : Boolean;
     LengthA : AlglibInteger;
     BlockSize : AlglibInteger;
     var T : TReal2DArray;
     var WORK : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // fill beginning of new column with zeros,
    // load 1.0 in the first non-zero element
    //
    K:=0;
    while K<=BlockSize-1 do
    begin
        if ColumnwiseA then
        begin
            I:=0;
            while I<=K-1 do
            begin
                A[I,K] := 0;
                Inc(I);
            end;
        end
        else
        begin
            I:=0;
            while I<=K-1 do
            begin
                A[K,I] := 0;
                Inc(I);
            end;
        end;
        A[K,K] := 1;
        Inc(K);
    end;
    
    //
    // Calculate Gram matrix of A
    //
    I:=0;
    while I<=BlockSize-1 do
    begin
        J:=0;
        while J<=BlockSize-1 do
        begin
            T[I,BlockSize+J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    K:=0;
    while K<=LengthA-1 do
    begin
        J:=1;
        while J<=BlockSize-1 do
        begin
            if ColumnwiseA then
            begin
                V := A[K,J];
                if AP_FP_Neq(V,0) then
                begin
                    APVAdd(@T[J][0], BlockSize, BlockSize+J-1, @A[K][0], 0, J-1, V);
                end;
            end
            else
            begin
                V := A[J,K];
                if AP_FP_Neq(V,0) then
                begin
                    i1_ := (0) - (BlockSize);
                    for i_ := BlockSize to BlockSize+J-1 do
                    begin
                        T[J,i_] := T[J,i_] + V*A[i_+i1_,K];
                    end;
                end;
            end;
            Inc(J);
        end;
        Inc(K);
    end;
    
    //
    // Prepare Y (stored in TmpA) and T (stored in TmpT)
    //
    K:=0;
    while K<=BlockSize-1 do
    begin
        
        //
        // fill non-zero part of T, use pre-calculated Gram matrix
        //
        APVMove(@WORK[0], 0, K-1, @T[K][0], BlockSize, BlockSize+K-1);
        I:=0;
        while I<=K-1 do
        begin
            V := APVDotProduct(@T[I][0], I, K-1, @WORK[0], I, K-1);
            T[I,K] := -Tau[K]*V;
            Inc(I);
        end;
        T[K,K] := -Tau[K];
        
        //
        // Rest of T is filled by zeros
        //
        I:=K+1;
        while I<=BlockSize-1 do
        begin
            T[I,K] := 0;
            Inc(I);
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Generate block reflector (complex):
* fill unused parts of reflectors matrix by zeros
* fill diagonal of reflectors matrix by ones
* generate triangular factor T


  -- ALGLIB routine --
     17.02.2010
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixBlockReflector(var A : TComplex2DArray;
     var Tau : TComplex1DArray;
     ColumnwiseA : Boolean;
     LengthA : AlglibInteger;
     BlockSize : AlglibInteger;
     var T : TComplex2DArray;
     var WORK : TComplex1DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    
    //
    // Prepare Y (stored in TmpA) and T (stored in TmpT)
    //
    K:=0;
    while K<=BlockSize-1 do
    begin
        
        //
        // fill beginning of new column with zeros,
        // load 1.0 in the first non-zero element
        //
        if ColumnwiseA then
        begin
            I:=0;
            while I<=K-1 do
            begin
                A[I,K] := C_Complex(0);
                Inc(I);
            end;
        end
        else
        begin
            I:=0;
            while I<=K-1 do
            begin
                A[K,I] := C_Complex(0);
                Inc(I);
            end;
        end;
        A[K,K] := C_Complex(1);
        
        //
        // fill non-zero part of T,
        //
        I:=0;
        while I<=K-1 do
        begin
            if ColumnwiseA then
            begin
                V := C_Complex(0.0);
                for i_ := K to LengthA-1 do
                begin
                    V := C_Add(V,C_Mul(Conj(A[i_,I]),A[i_,K]));
                end;
            end
            else
            begin
                V := C_Complex(0.0);
                for i_ := K to LengthA-1 do
                begin
                    V := C_Add(V,C_Mul(A[I,i_],Conj(A[K,i_])));
                end;
            end;
            WORK[I] := V;
            Inc(I);
        end;
        I:=0;
        while I<=K-1 do
        begin
            V := C_Complex(0.0);
            for i_ := I to K-1 do
            begin
                V := C_Add(V,C_Mul(T[I,i_],WORK[i_]));
            end;
            T[I,K] := C_Opposite(C_Mul(Tau[K],V));
            Inc(I);
        end;
        T[K,K] := C_Opposite(Tau[K]);
        
        //
        // Rest of T is filled by zeros
        //
        I:=K+1;
        while I<=BlockSize-1 do
        begin
            T[I,K] := C_Complex(0);
            Inc(I);
        end;
        Inc(K);
    end;
end;


end.