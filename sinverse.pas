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
unit sinverse;
interface
uses Math, Sysutils, Ap, sblas, ldlt;

function SMatrixLDLTInverse(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function SMatrixInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function InverseLDLT(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function InverseSymmetricIndefinite(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;

implementation

(*************************************************************************
Inversion of a symmetric indefinite matrix

The algorithm gets an LDLT-decomposition as an input, generates matrix A^-1
and saves the lower or upper triangle of an inverse matrix depending on the
input (U*D*U' or L*D*L').

Input parameters:
    A       -   LDLT-decomposition of the matrix,
                Output of subroutine SMatrixLDLT.
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then the symmetric matrix
                is given as decomposition A = U*D*U' and this decomposition
                is stored in the upper triangle of matrix A and on the main
                diagonal, and the lower triangle of matrix A is not used.
    Pivots  -   a table of permutations, output of subroutine SMatrixLDLT.

Output parameters:
    A       -   inverse of the matrix, whose LDLT-decomposition was stored
                in matrix A as a subroutine input.
                Array with elements [0..N-1, 0..N-1].
                If IsUpper = True, then A contains the upper triangle of
                matrix A^-1, and the elements below the main diagonal are
                not used nor changed. The same applies if IsUpper = False.

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular and could not be inverted.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************)
function SMatrixLDLTInverse(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    WORK : TReal1DArray;
    WORK2 : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    KP : AlglibInteger;
    KSTEP : AlglibInteger;
    AK : Double;
    AKKP1 : Double;
    AKP1 : Double;
    D : Double;
    T : Double;
    TEMP : Double;
    KM1 : AlglibInteger;
    KP1 : AlglibInteger;
    L : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    SetLength(Work, N+1);
    SetLength(Work2, N+1);
    Result := True;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // Check that the diagonal matrix D is nonsingular.
    //
    I:=0;
    while I<=N-1 do
    begin
        if (Pivots[I]>=0) and AP_FP_Eq(A[I,I],0) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    if IsUPPER then
    begin
        
        //
        // Compute inv(A) from the factorization A = U*D*U'.
        //
        // K+1 is the main loop index, increasing from 1 to N in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        K := 0;
        while K<=N-1 do
        begin
            if Pivots[K]>=0 then
            begin
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                A[K,K] := 1/A[K,K];
                
                //
                // Compute column K+1 of the inverse.
                //
                if K>0 then
                begin
                    i1_ := (0) - (1);
                    for i_ := 1 to K do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, 1-1, K+1-1-1, WORK, -1, WORK2);
                    i1_ := (1) - (0);
                    for i_ := 0 to K-1 do
                    begin
                        A[i_,K] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK2[0], 1, K, @WORK[0], 1, K);
                    A[K,K] := A[K,K]-V;
                end;
                KSTEP := 1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                T := ABSReal(A[K,K+1]);
                AK := A[K,K]/T;
                AKP1 := A[K+1,K+1]/T;
                AKKP1 := A[K,K+1]/T;
                D := T*(AK*AKP1-1);
                A[K,K] := AKP1/D;
                A[K+1,K+1] := AK/D;
                A[K,K+1] := -AKKP1/D;
                
                //
                // Compute columns K+1 and K+1+1 of the inverse.
                //
                if K>0 then
                begin
                    i1_ := (0) - (1);
                    for i_ := 1 to K do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, 0, K-1, WORK, -1, WORK2);
                    i1_ := (1) - (0);
                    for i_ := 0 to K-1 do
                    begin
                        A[i_,K] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, K, @WORK2[0], 1, K);
                    A[K,K] := A[K,K]-V;
                    V := 0.0;
                    for i_ := 0 to K-1 do
                    begin
                        V := V + A[i_,K]*A[i_,K+1];
                    end;
                    A[K,K+1] := A[K,K+1]-V;
                    i1_ := (0) - (1);
                    for i_ := 1 to K do
                    begin
                        WORK[i_] := A[i_+i1_,K+1];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, 0, K-1, WORK, -1, WORK2);
                    i1_ := (1) - (0);
                    for i_ := 0 to K-1 do
                    begin
                        A[i_,K+1] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, K, @WORK2[0], 1, K);
                    A[K+1,K+1] := A[K+1,K+1]-V;
                end;
                KSTEP := 2;
            end;
            if Pivots[K]>=0 then
            begin
                KP := Pivots[K];
            end
            else
            begin
                KP := N+Pivots[K];
            end;
            if KP<>K then
            begin
                
                //
                // Interchange rows and columns K and KP in the leading
                // submatrix
                //
                i1_ := (0) - (1);
                for i_ := 1 to KP do
                begin
                    WORK[i_] := A[i_+i1_,K];
                end;
                for i_ := 0 to KP-1 do
                begin
                    A[i_,K] := A[i_,KP];
                end;
                i1_ := (1) - (0);
                for i_ := 0 to KP-1 do
                begin
                    A[i_,KP] := WORK[i_+i1_];
                end;
                i1_ := (KP+1) - (1);
                for i_ := 1 to K-1-KP do
                begin
                    WORK[i_] := A[i_+i1_,K];
                end;
                for i_ := KP+1 to K-1 do
                begin
                    A[i_,K] := A[KP,i_];
                end;
                APVMove(@A[KP][0], KP+1, K-1, @WORK[0], 1, K-1-KP);
                TEMP := A[K,K];
                A[K,K] := A[KP,KP];
                A[KP,KP] := TEMP;
                if KSTEP=2 then
                begin
                    TEMP := A[K,K+1];
                    A[K,K+1] := A[KP,K+1];
                    A[KP,K+1] := TEMP;
                end;
            end;
            K := K+KSTEP;
        end;
    end
    else
    begin
        
        //
        // Compute inv(A) from the factorization A = L*D*L'.
        //
        // K is the main loop index, increasing from 0 to N-1 in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        K := N-1;
        while K>=0 do
        begin
            if Pivots[K]>=0 then
            begin
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                A[K,K] := 1/A[K,K];
                
                //
                // Compute column K+1 of the inverse.
                //
                if K<N-1 then
                begin
                    i1_ := (K+1) - (1);
                    for i_ := 1 to N-K-1 do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, K+1, N-1, WORK, -1, WORK2);
                    i1_ := (1) - (K+1);
                    for i_ := K+1 to N-1 do
                    begin
                        A[i_,K] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, N-K-1, @WORK2[0], 1, N-K-1);
                    A[K,K] := A[K,K]-V;
                end;
                KSTEP := 1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                T := ABSReal(A[K,K-1]);
                AK := A[K-1,K-1]/T;
                AKP1 := A[K,K]/T;
                AKKP1 := A[K,K-1]/T;
                D := T*(AK*AKP1-1);
                A[K-1,K-1] := AKP1/D;
                A[K,K] := AK/D;
                A[K,K-1] := -AKKP1/D;
                
                //
                // Compute columns K+1-1 and K+1 of the inverse.
                //
                if K<N-1 then
                begin
                    i1_ := (K+1) - (1);
                    for i_ := 1 to N-K-1 do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, K+1, N-1, WORK, -1, WORK2);
                    i1_ := (1) - (K+1);
                    for i_ := K+1 to N-1 do
                    begin
                        A[i_,K] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, N-K-1, @WORK2[0], 1, N-K-1);
                    A[K,K] := A[K,K]-V;
                    V := 0.0;
                    for i_ := K+1 to N-1 do
                    begin
                        V := V + A[i_,K]*A[i_,K-1];
                    end;
                    A[K,K-1] := A[K,K-1]-V;
                    i1_ := (K+1) - (1);
                    for i_ := 1 to N-K-1 do
                    begin
                        WORK[i_] := A[i_+i1_,K-1];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, K+1, N-1, WORK, -1, WORK2);
                    i1_ := (1) - (K+1);
                    for i_ := K+1 to N-1 do
                    begin
                        A[i_,K-1] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, N-K-1, @WORK2[0], 1, N-K-1);
                    A[K-1,K-1] := A[K-1,K-1]-V;
                end;
                KSTEP := 2;
            end;
            if Pivots[K]>=0 then
            begin
                KP := Pivots[K];
            end
            else
            begin
                KP := Pivots[K]+N;
            end;
            if KP<>K then
            begin
                
                //
                // Interchange rows and columns K and KP
                //
                if KP<N-1 then
                begin
                    i1_ := (KP+1) - (1);
                    for i_ := 1 to N-KP-1 do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    for i_ := KP+1 to N-1 do
                    begin
                        A[i_,K] := A[i_,KP];
                    end;
                    i1_ := (1) - (KP+1);
                    for i_ := KP+1 to N-1 do
                    begin
                        A[i_,KP] := WORK[i_+i1_];
                    end;
                end;
                i1_ := (K+1) - (1);
                for i_ := 1 to KP-K-1 do
                begin
                    WORK[i_] := A[i_+i1_,K];
                end;
                for i_ := K+1 to KP-1 do
                begin
                    A[i_,K] := A[KP,i_];
                end;
                APVMove(@A[KP][0], K+1, KP-1, @WORK[0], 1, KP-K-1);
                TEMP := A[K,K];
                A[K,K] := A[KP,KP];
                A[KP,KP] := TEMP;
                if KSTEP=2 then
                begin
                    TEMP := A[K,K-1];
                    A[K,K-1] := A[KP,K-1];
                    A[KP,K-1] := TEMP;
                end;
            end;
            K := K-KSTEP;
        end;
    end;
end;


(*************************************************************************
Inversion of a symmetric indefinite matrix

Given a lower or upper triangle of matrix A, the algorithm generates
matrix A^-1 and saves the lower or upper triangle depending on the input.

Input parameters:
    A       -   matrix to be inverted (upper or lower triangle).
                Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then the upper
                triangle of matrix A is given, otherwise the lower
                triangle is given.

Output parameters:
    A       -   inverse of matrix A.
                Array with elements [0..N-1, 0..N-1].
                If IsUpper = True, then A contains the upper triangle of
                matrix A^-1, and the elements below the main diagonal are
                not used nor changed.
                The same applies if IsUpper = False.

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular and could not be inverted.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************)
function SMatrixInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    Pivots : TInteger1DArray;
begin
    SMatrixLDLT(A, N, IsUpper, Pivots);
    Result := SMatrixLDLTInverse(A, Pivots, N, IsUpper);
end;


function InverseLDLT(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    WORK : TReal1DArray;
    WORK2 : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    KP : AlglibInteger;
    KSTEP : AlglibInteger;
    AK : Double;
    AKKP1 : Double;
    AKP1 : Double;
    D : Double;
    T : Double;
    TEMP : Double;
    KM1 : AlglibInteger;
    KP1 : AlglibInteger;
    L : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    SetLength(Work, N+1);
    SetLength(Work2, N+1);
    Result := True;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // Check that the diagonal matrix D is nonsingular.
    //
    I:=1;
    while I<=N do
    begin
        if (Pivots[I]>0) and AP_FP_Eq(A[I,I],0) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    if IsUPPER then
    begin
        
        //
        // Compute inv(A) from the factorization A = U*D*U'.
        //
        // K is the main loop index, increasing from 1 to N in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        K := 1;
        while K<=N do
        begin
            if Pivots[K]>0 then
            begin
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                A[K,K] := 1/A[K,K];
                
                //
                // Compute column K of the inverse.
                //
                if K>1 then
                begin
                    KM1 := K-1;
                    for i_ := 1 to KM1 do
                    begin
                        WORK[i_] := A[i_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, 1, K-1, WORK, -1, WORK2);
                    for i_ := 1 to KM1 do
                    begin
                        A[i_,K] := WORK2[i_];
                    end;
                    V := APVDotProduct(@WORK2[0], 1, KM1, @WORK[0], 1, KM1);
                    A[K,K] := A[K,K]-V;
                end;
                KSTEP := 1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                T := ABSReal(A[K,K+1]);
                AK := A[K,K]/T;
                AKP1 := A[K+1,K+1]/T;
                AKKP1 := A[K,K+1]/T;
                D := T*(AK*AKP1-1);
                A[K,K] := AKP1/D;
                A[K+1,K+1] := AK/D;
                A[K,K+1] := -AKKP1/D;
                
                //
                // Compute columns K and K+1 of the inverse.
                //
                if K>1 then
                begin
                    KM1 := K-1;
                    KP1 := K+1;
                    for i_ := 1 to KM1 do
                    begin
                        WORK[i_] := A[i_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, 1, K-1, WORK, -1, WORK2);
                    for i_ := 1 to KM1 do
                    begin
                        A[i_,K] := WORK2[i_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, KM1, @WORK2[0], 1, KM1);
                    A[K,K] := A[K,K]-V;
                    V := 0.0;
                    for i_ := 1 to KM1 do
                    begin
                        V := V + A[i_,K]*A[i_,KP1];
                    end;
                    A[K,K+1] := A[K,K+1]-V;
                    for i_ := 1 to KM1 do
                    begin
                        WORK[i_] := A[i_,KP1];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, 1, K-1, WORK, -1, WORK2);
                    for i_ := 1 to KM1 do
                    begin
                        A[i_,KP1] := WORK2[i_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, KM1, @WORK2[0], 1, KM1);
                    A[K+1,K+1] := A[K+1,K+1]-V;
                end;
                KSTEP := 2;
            end;
            KP := ABSInt(Pivots[K]);
            if KP<>K then
            begin
                
                //
                // Interchange rows and columns K and KP in the leading
                // submatrix A(1:k+1,1:k+1)
                //
                L := KP-1;
                for i_ := 1 to L do
                begin
                    WORK[i_] := A[i_,K];
                end;
                for i_ := 1 to L do
                begin
                    A[i_,K] := A[i_,KP];
                end;
                for i_ := 1 to L do
                begin
                    A[i_,KP] := WORK[i_];
                end;
                L := K-KP-1;
                I1 := KP+1;
                I2 := K-1;
                i1_ := (I1) - (1);
                for i_ := 1 to L do
                begin
                    WORK[i_] := A[i_+i1_,K];
                end;
                for i_ := I1 to I2 do
                begin
                    A[i_,K] := A[KP,i_];
                end;
                APVMove(@A[KP][0], I1, I2, @WORK[0], 1, L);
                TEMP := A[K,K];
                A[K,K] := A[KP,KP];
                A[KP,KP] := TEMP;
                if KSTEP=2 then
                begin
                    TEMP := A[K,K+1];
                    A[K,K+1] := A[KP,K+1];
                    A[KP,K+1] := TEMP;
                end;
            end;
            K := K+KSTEP;
        end;
    end
    else
    begin
        
        //
        // Compute inv(A) from the factorization A = L*D*L'.
        //
        // K is the main loop index, increasing from 1 to N in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        K := N;
        while K>=1 do
        begin
            if Pivots[K]>0 then
            begin
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                A[K,K] := 1/A[K,K];
                
                //
                // Compute column K of the inverse.
                //
                if K<N then
                begin
                    KP1 := K+1;
                    KM1 := K-1;
                    L := N-K;
                    i1_ := (KP1) - (1);
                    for i_ := 1 to L do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, K+1, N, WORK, -1, WORK2);
                    i1_ := (1) - (KP1);
                    for i_ := KP1 to N do
                    begin
                        A[i_,K] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, L, @WORK2[0], 1, L);
                    A[K,K] := A[K,K]-V;
                end;
                KSTEP := 1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                T := ABSReal(A[K,K-1]);
                AK := A[K-1,K-1]/T;
                AKP1 := A[K,K]/T;
                AKKP1 := A[K,K-1]/T;
                D := T*(AK*AKP1-1);
                A[K-1,K-1] := AKP1/D;
                A[K,K] := AK/D;
                A[K,K-1] := -AKKP1/D;
                
                //
                // Compute columns K-1 and K of the inverse.
                //
                if K<N then
                begin
                    KP1 := K+1;
                    KM1 := K-1;
                    L := N-K;
                    i1_ := (KP1) - (1);
                    for i_ := 1 to L do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, K+1, N, WORK, -1, WORK2);
                    i1_ := (1) - (KP1);
                    for i_ := KP1 to N do
                    begin
                        A[i_,K] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, L, @WORK2[0], 1, L);
                    A[K,K] := A[K,K]-V;
                    V := 0.0;
                    for i_ := KP1 to N do
                    begin
                        V := V + A[i_,K]*A[i_,KM1];
                    end;
                    A[K,K-1] := A[K,K-1]-V;
                    i1_ := (KP1) - (1);
                    for i_ := 1 to L do
                    begin
                        WORK[i_] := A[i_+i1_,KM1];
                    end;
                    SymmetricMatrixVectorMultiply(A, IsUpper, K+1, N, WORK, -1, WORK2);
                    i1_ := (1) - (KP1);
                    for i_ := KP1 to N do
                    begin
                        A[i_,KM1] := WORK2[i_+i1_];
                    end;
                    V := APVDotProduct(@WORK[0], 1, L, @WORK2[0], 1, L);
                    A[K-1,K-1] := A[K-1,K-1]-V;
                end;
                KSTEP := 2;
            end;
            KP := ABSInt(Pivots[K]);
            if KP<>K then
            begin
                
                //
                // Interchange rows and columns K and KP in the trailing
                // submatrix A(k-1:n,k-1:n)
                //
                if KP<N then
                begin
                    L := N-KP;
                    KP1 := KP+1;
                    i1_ := (KP1) - (1);
                    for i_ := 1 to L do
                    begin
                        WORK[i_] := A[i_+i1_,K];
                    end;
                    for i_ := KP1 to N do
                    begin
                        A[i_,K] := A[i_,KP];
                    end;
                    i1_ := (1) - (KP1);
                    for i_ := KP1 to N do
                    begin
                        A[i_,KP] := WORK[i_+i1_];
                    end;
                end;
                L := KP-K-1;
                I1 := K+1;
                I2 := KP-1;
                i1_ := (I1) - (1);
                for i_ := 1 to L do
                begin
                    WORK[i_] := A[i_+i1_,K];
                end;
                for i_ := I1 to I2 do
                begin
                    A[i_,K] := A[KP,i_];
                end;
                APVMove(@A[KP][0], I1, I2, @WORK[0], 1, L);
                TEMP := A[K,K];
                A[K,K] := A[KP,KP];
                A[KP,KP] := TEMP;
                if KSTEP=2 then
                begin
                    TEMP := A[K,K-1];
                    A[K,K-1] := A[KP,K-1];
                    A[KP,K-1] := TEMP;
                end;
            end;
            K := K-KSTEP;
        end;
    end;
end;


function InverseSymmetricIndefinite(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    Pivots : TInteger1DArray;
begin
    LDLTDecomposition(A, N, IsUpper, Pivots);
    Result := InverseLDLT(A, Pivots, N, IsUpper);
end;


end.