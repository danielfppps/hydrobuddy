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
unit ldlt;
interface
uses Math, Sysutils, Ap;

procedure SMatrixLDLT(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Pivots : TInteger1DArray);
procedure LDLTDecomposition(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Pivots : TInteger1DArray);

implementation

(*************************************************************************
LDLTDecomposition of a symmetric matrix

The algorithm represents a symmetric matrix (which is not necessarily
positive definite) as A=L*D*L' or A = U*D*U', where D is a block-diagonal
matrix with blocks 1x1 or 2x2, matrix L (matrix U) is a product of lower
(upper) triangular matrices with unit diagonal and permutation matrices.

Input parameters:
    A       -   factorized matrix, array with elements [0..N-1, 0..N-1].
                If IsUpper – True, then the upper triangle contains
                elements of symmetric matrix A, and the lower triangle is
                not used.
                The same applies if IsUpper = False.
    N       -   size of factorized matrix.
    IsUpper -   parameter which shows a method of matrix definition (lower
                or upper triangle).

Output parameters:
    A       -   matrices D and U, if IsUpper = True, or L, if IsUpper = False,
                in compact form, replacing the upper (lower) triangle of
                matrix A. In that case, the elements under (over) the main
                diagonal are not used nor modified.
    Pivots  -   tables of performed permutations (see below).

If IsUpper = True, then A = U*D*U', U = P(n)*U(n)*...*P(k)*U(k), where
P(k) is the permutation matrix, U(k) - upper triangular matrix with its
unit main diagonal and k decreases from n with step s which is equal to
1 or 2 (according to the size of the blocks of matrix D).

        (   I    v    0   )   k-s+1
U(k) =  (   0    I    0   )   s
        (   0    0    I   )   n-k-1
           k-s+1 s   n-k-1

If Pivots[k]>=0, then s=1, P(k) - permutation of rows k and Pivots[k], the
vectorv forming matrix U(k) is stored in elements A(0:k-1,k), D(k) replaces
A(k,k). If Pivots[k]=Pivots[k-1]<0 then s=2, P(k) - permutation of rows k-1
and N+Pivots[k-1], the vector v forming matrix U(k) is stored in elements
A(0:k-1,k:k+1), the upper triangle of block D(k) is stored in A(k,k),
A(k,k+1) and A(k+1,k+1).

If IsUpper = False, then A = L*D*L', L=P(0)*L(0)*...*P(k)*L(k), where P(k)
is the permutation matrix, L(k) – lower triangular matrix with unit main
diagonal and k decreases from 1 with step s which is equal to 1 or 2
(according to the size of the blocks of matrix D).

        (   I    0     0   )  k-1
L(k) =  (   0    I     0   )  s
        (   0    v     I   )  n-k-s+1
           k-1   s  n-k-s+1

If Pivots[k]>=0 then s=1, P(k) – permutation of rows k and Pivots[k], the
vector v forming matrix L(k) is stored in elements A(k+1:n-1,k), D(k)
replaces A(k,k). If Pivots[k]=Pivots[k+1]<0 then s=2, P(k) - permutation
of rows k+1 and N+Pivots[k+1], the vector v forming matrix L(k) is stored
in elements A(k+2:n-1,k:k+1), the lower triangle of block D(k) is stored in
A(k,k), A(k+1,k) and A(k+1,k+1).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1999
*************************************************************************)
procedure SMatrixLDLT(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    IMAX : AlglibInteger;
    J : AlglibInteger;
    JMAX : AlglibInteger;
    K : AlglibInteger;
    KK : AlglibInteger;
    KP : AlglibInteger;
    KSTEP : AlglibInteger;
    ABSAKK : Double;
    ALPHA : Double;
    COLMAX : Double;
    D11 : Double;
    D12 : Double;
    D21 : Double;
    D22 : Double;
    R1 : Double;
    ROWMAX : Double;
    T : Double;
    WK : Double;
    WKM1 : Double;
    WKP1 : Double;
    II : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    VV : Double;
    Temp : TReal1DArray;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, N-1+1);
    SetLength(Temp, N-1+1);
    
    //
    // Initialize ALPHA for use in choosing pivot block size.
    //
    ALPHA := (1+SQRT(17))/8;
    if IsUpper then
    begin
        
        //
        // Factorize A as U*D*U' using the upper triangle of A
        //
        //
        // K is the main loop index, decreasing from N to 1 in steps of
        // 1 or 2
        //
        K := N-1;
        while K>=0 do
        begin
            KSTEP := 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            ABSAKK := ABSReal(A[K,K]);
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K+1, and COLMAX is its absolute value
            //
            if K>0 then
            begin
                IMAX := 1;
                II:=2;
                while II<=K do
                begin
                    if AP_FP_Greater(AbsReal(A[II-1,K]),AbsReal(A[IMAX-1,K])) then
                    begin
                        IMAX := II;
                    end;
                    Inc(II);
                end;
                COLMAX := ABSReal(A[IMAX-1,K]);
            end
            else
            begin
                COLMAX := 0;
            end;
            if AP_FP_Eq(Max(ABSAKK, COLMAX),0) then
            begin
                
                //
                // Column K is zero
                //
                KP := K;
            end
            else
            begin
                if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX) then
                begin
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    KP := K;
                end
                else
                begin
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    JMAX := IMAX+1;
                    II:=IMAX+2;
                    while II<=K+1 do
                    begin
                        if AP_FP_Greater(AbsReal(A[IMAX-1,II-1]),AbsReal(A[IMAX-1,JMAX-1])) then
                        begin
                            JMAX := II;
                        end;
                        Inc(II);
                    end;
                    ROWMAX := ABSReal(A[IMAX-1,JMAX-1]);
                    if IMAX>1 then
                    begin
                        JMAX := 1;
                        II:=2;
                        while II<=IMAX-1 do
                        begin
                            if AP_FP_Greater(AbsReal(A[II-1,IMAX-1]),AbsReal(A[JMAX-1,IMAX-1])) then
                            begin
                                JMAX := II;
                            end;
                            Inc(II);
                        end;
                        ROWMAX := Max(ROWMAX, ABSReal(A[JMAX-1,IMAX-1]));
                    end;
                    VV := COLMAX/ROWMAX;
                    if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX*VV) then
                    begin
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        KP := K;
                    end
                    else
                    begin
                        if AP_FP_Greater_Eq(ABSReal(A[IMAX-1,IMAX-1]),ALPHA*ROWMAX) then
                        begin
                            
                            //
                            // interchange rows and columns K and IMAX, use 1-by-1
                            // pivot block
                            //
                            KP := IMAX-1;
                        end
                        else
                        begin
                            
                            //
                            // interchange rows and columns K-1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            KP := IMAX-1;
                            KSTEP := 2;
                        end;
                    end;
                end;
                KK := K+1-KSTEP;
                if KP+1<>KK+1 then
                begin
                    
                    //
                    // Interchange rows and columns KK and KP+1 in the leading
                    // submatrix A(0:K,0:K)
                    //
                    for i_ := 0 to KP-1 do
                    begin
                        Temp[i_] := A[i_,KK];
                    end;
                    for i_ := 0 to KP-1 do
                    begin
                        A[i_,KK] := A[i_,KP];
                    end;
                    for i_ := 0 to KP-1 do
                    begin
                        A[i_,KP] := Temp[i_];
                    end;
                    for i_ := KP+1 to KK-1 do
                    begin
                        Temp[i_] := A[i_,KK];
                    end;
                    for i_ := KP+1 to KK-1 do
                    begin
                        A[i_,KK] := A[KP,i_];
                    end;
                    APVMove(@A[KP][0], KP+1, KK-1, @Temp[0], KP+1, KK-1);
                    T := A[KK,KK];
                    A[KK,KK] := A[KP,KP];
                    A[KP,KP] := T;
                    if KSTEP=2 then
                    begin
                        T := A[K-1,K];
                        A[K-1,K] := A[KP,K];
                        A[KP,K] := T;
                    end;
                end;
                
                //
                // Update the leading submatrix
                //
                if KSTEP=1 then
                begin
                    
                    //
                    // 1-by-1 pivot block D(k): column k now holds
                    //
                    // W(k) = U(k)*D(k)
                    //
                    // where U(k) is the k-th column of U
                    //
                    // Perform a rank-1 update of A(1:k-1,1:k-1) as
                    //
                    // A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
                    //
                    R1 := 1/A[K,K];
                    I:=0;
                    while I<=K-1 do
                    begin
                        VV := -R1*A[I,K];
                        for i_ := I to K-1 do
                        begin
                            A[I,i_] := A[I,i_] + VV*A[i_,K];
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Store U(K+1) in column K+1
                    //
                    for i_ := 0 to K-1 do
                    begin
                        A[i_,K] := R1*A[i_,K];
                    end;
                end
                else
                begin
                    
                    //
                    // 2-by-2 pivot block D(k): columns k and k-1 now hold
                    //
                    // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                    //
                    // where U(k) and U(k-1) are the k-th and (k-1)-th columns
                    // of U
                    //
                    // Perform a rank-2 update of A(1:k-2,1:k-2) as
                    //
                    // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
                    //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
                    //
                    if K>1 then
                    begin
                        D12 := A[K-1,K];
                        D22 := A[K-1,K-1]/D12;
                        D11 := A[K,K]/D12;
                        T := 1/(D11*D22-1);
                        D12 := T/D12;
                        J:=K-2;
                        while J>=0 do
                        begin
                            WKM1 := D12*(D11*A[J,K-1]-A[J,K]);
                            WK := D12*(D22*A[J,K]-A[J,K-1]);
                            for i_ := 0 to J do
                            begin
                                A[i_,J] := A[i_,J] - WK*A[i_,K];
                            end;
                            for i_ := 0 to J do
                            begin
                                A[i_,J] := A[i_,J] - WKM1*A[i_,K-1];
                            end;
                            A[J,K] := WK;
                            A[J,K-1] := WKM1;
                            Dec(J);
                        end;
                    end;
                end;
            end;
            
            //
            // Store details of the interchanges in IPIV
            //
            if KSTEP=1 then
            begin
                Pivots[K] := KP;
            end
            else
            begin
                Pivots[K] := KP-N;
                Pivots[K-1] := KP-N;
            end;
            
            //
            // Decrease K+1 and return to the start of the main loop
            //
            K := K-KSTEP;
        end;
    end
    else
    begin
        
        //
        // Factorize A as L*D*L' using the lower triangle of A
        //
        // K+1 is the main loop index, increasing from 1 to N in steps of
        // 1 or 2
        //
        K := 0;
        while K<=N-1 do
        begin
            KSTEP := 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            ABSAKK := ABSReal(A[K,K]);
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K+1, and COLMAX is its absolute value
            //
            if K<N-1 then
            begin
                IMAX := K+1+1;
                II:=K+1+2;
                while II<=N do
                begin
                    if AP_FP_Greater(AbsReal(A[II-1,K]),AbsReal(A[IMAX-1,K])) then
                    begin
                        IMAX := II;
                    end;
                    Inc(II);
                end;
                COLMAX := ABSReal(A[IMAX-1,K]);
            end
            else
            begin
                COLMAX := 0;
            end;
            if AP_FP_Eq(Max(ABSAKK, COLMAX),0) then
            begin
                
                //
                // Column K+1 is zero
                //
                KP := K;
            end
            else
            begin
                if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX) then
                begin
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    KP := K;
                end
                else
                begin
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    JMAX := K+1;
                    II:=K+1+1;
                    while II<=IMAX-1 do
                    begin
                        if AP_FP_Greater(AbsReal(A[IMAX-1,II-1]),AbsReal(A[IMAX-1,JMAX-1])) then
                        begin
                            JMAX := II;
                        end;
                        Inc(II);
                    end;
                    ROWMAX := ABSReal(A[IMAX-1,JMAX-1]);
                    if IMAX<N then
                    begin
                        JMAX := IMAX+1;
                        II:=IMAX+2;
                        while II<=N do
                        begin
                            if AP_FP_Greater(AbsReal(A[II-1,IMAX-1]),AbsReal(A[JMAX-1,IMAX-1])) then
                            begin
                                JMAX := II;
                            end;
                            Inc(II);
                        end;
                        ROWMAX := Max(ROWMAX, ABSReal(A[JMAX-1,IMAX-1]));
                    end;
                    VV := COLMAX/ROWMAX;
                    if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX*VV) then
                    begin
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        KP := K;
                    end
                    else
                    begin
                        if AP_FP_Greater_Eq(ABSReal(A[IMAX-1,IMAX-1]),ALPHA*ROWMAX) then
                        begin
                            
                            //
                            // interchange rows and columns K+1 and IMAX, use 1-by-1
                            // pivot block
                            //
                            KP := IMAX-1;
                        end
                        else
                        begin
                            
                            //
                            // interchange rows and columns K+1+1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            KP := IMAX-1;
                            KSTEP := 2;
                        end;
                    end;
                end;
                KK := K+KSTEP-1;
                if KP<>KK then
                begin
                    
                    //
                    //              Interchange rows and columns KK+1 and KP+1 in the trailing
                    //              submatrix A(K+1:n,K+1:n)
                    //
                    if KP+1<N then
                    begin
                        for i_ := KP+1 to N-1 do
                        begin
                            Temp[i_] := A[i_,KK];
                        end;
                        for i_ := KP+1 to N-1 do
                        begin
                            A[i_,KK] := A[i_,KP];
                        end;
                        for i_ := KP+1 to N-1 do
                        begin
                            A[i_,KP] := Temp[i_];
                        end;
                    end;
                    for i_ := KK+1 to KP-1 do
                    begin
                        Temp[i_] := A[i_,KK];
                    end;
                    for i_ := KK+1 to KP-1 do
                    begin
                        A[i_,KK] := A[KP,i_];
                    end;
                    APVMove(@A[KP][0], KK+1, KP-1, @Temp[0], KK+1, KP-1);
                    T := A[KK,KK];
                    A[KK,KK] := A[KP,KP];
                    A[KP,KP] := T;
                    if KSTEP=2 then
                    begin
                        T := A[K+1,K];
                        A[K+1,K] := A[KP,K];
                        A[KP,K] := T;
                    end;
                end;
                
                //
                // Update the trailing submatrix
                //
                if KSTEP=1 then
                begin
                    
                    //
                    // 1-by-1 pivot block D(K+1): column K+1 now holds
                    //
                    // W(K+1) = L(K+1)*D(K+1)
                    //
                    // where L(K+1) is the K+1-th column of L
                    //
                    if K+1<N then
                    begin
                        
                        //
                        // Perform a rank-1 update of A(K+1+1:n,K+1+1:n) as
                        //
                        // A := A - L(K+1)*D(K+1)*L(K+1)' = A - W(K+1)*(1/D(K+1))*W(K+1)'
                        //
                        D11 := 1/A[K+1-1,K+1-1];
                        II:=K+1;
                        while II<=N-1 do
                        begin
                            VV := -D11*A[II,K];
                            for i_ := K+1 to II do
                            begin
                                A[II,i_] := A[II,i_] + VV*A[i_,K];
                            end;
                            Inc(II);
                        end;
                        
                        //
                        // Store L(K+1) in column K+1
                        //
                        for i_ := K+1 to N-1 do
                        begin
                            A[i_,K] := D11*A[i_,K];
                        end;
                    end;
                end
                else
                begin
                    
                    //
                    // 2-by-2 pivot block D(K+1)
                    //
                    if K<N-2 then
                    begin
                        
                        //
                        // Perform a rank-2 update of A(K+1+2:n,K+1+2:n) as
                        //
                        // A := A - ( (A(K+1) A(K+1+1))*D(K+1)**(-1) ) * (A(K+1) A(K+1+1))'
                        //
                        // where L(K+1) and L(K+1+1) are the K+1-th and (K+1+1)-th
                        // columns of L
                        //
                        D21 := A[K+1,K];
                        D11 := A[K+1,K+1]/D21;
                        D22 := A[K,K]/D21;
                        T := 1/(D11*D22-1);
                        D21 := T/D21;
                        J:=K+2;
                        while J<=N-1 do
                        begin
                            WK := D21*(D11*A[J,K]-A[J,K+1]);
                            WKP1 := D21*(D22*A[J,K+1]-A[J,K]);
                            for i_ := J to N-1 do
                            begin
                                A[i_,J] := A[i_,J] - WK*A[i_,K];
                            end;
                            for i_ := J to N-1 do
                            begin
                                A[i_,J] := A[i_,J] - WKP1*A[i_,K+1];
                            end;
                            A[J,K] := WK;
                            A[J,K+1] := WKP1;
                            Inc(J);
                        end;
                    end;
                end;
            end;
            
            //
            // Store details of the interchanges in IPIV
            //
            if KSTEP=1 then
            begin
                Pivots[K+1-1] := KP+1-1;
            end
            else
            begin
                Pivots[K+1-1] := KP+1-1-N;
                Pivots[K+1+1-1] := KP+1-1-N;
            end;
            
            //
            // Increase K+1 and return to the start of the main loop
            //
            K := K+KSTEP;
        end;
    end;
end;


procedure LDLTDecomposition(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    IMAX : AlglibInteger;
    J : AlglibInteger;
    JMAX : AlglibInteger;
    K : AlglibInteger;
    KK : AlglibInteger;
    KP : AlglibInteger;
    KSTEP : AlglibInteger;
    ABSAKK : Double;
    ALPHA : Double;
    COLMAX : Double;
    D11 : Double;
    D12 : Double;
    D21 : Double;
    D22 : Double;
    R1 : Double;
    ROWMAX : Double;
    T : Double;
    WK : Double;
    WKM1 : Double;
    WKP1 : Double;
    II : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    VV : Double;
    Temp : TReal1DArray;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, N+1);
    SetLength(Temp, N+1);
    
    //
    // Initialize ALPHA for use in choosing pivot block size.
    //
    ALPHA := (1+SQRT(17))/8;
    if IsUpper then
    begin
        
        //
        // Factorize A as U*D*U' using the upper triangle of A
        //
        //
        // K is the main loop index, decreasing from N to 1 in steps of
        // 1 or 2
        //
        K := N;
        while K>=1 do
        begin
            KSTEP := 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            ABSAKK := ABSReal(A[K,K]);
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K, and COLMAX is its absolute value
            //
            if K>1 then
            begin
                IMAX := 1;
                II:=2;
                while II<=K-1 do
                begin
                    if AP_FP_Greater(AbsReal(A[II,K]),AbsReal(A[IMAX,K])) then
                    begin
                        IMAX := II;
                    end;
                    Inc(II);
                end;
                COLMAX := ABSReal(A[IMAX,K]);
            end
            else
            begin
                COLMAX := 0;
            end;
            if AP_FP_Eq(Max(ABSAKK, COLMAX),0) then
            begin
                
                //
                // Column K is zero
                //
                KP := K;
            end
            else
            begin
                if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX) then
                begin
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    KP := K;
                end
                else
                begin
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    JMAX := IMAX+1;
                    II:=IMAX+2;
                    while II<=K do
                    begin
                        if AP_FP_Greater(AbsReal(A[IMAX,II]),AbsReal(A[IMAX,JMAX])) then
                        begin
                            JMAX := II;
                        end;
                        Inc(II);
                    end;
                    ROWMAX := ABSReal(A[IMAX,JMAX]);
                    if IMAX>1 then
                    begin
                        JMAX := 1;
                        II:=2;
                        while II<=IMAX-1 do
                        begin
                            if AP_FP_Greater(AbsReal(A[II,IMAX]),AbsReal(A[JMAX,IMAX])) then
                            begin
                                JMAX := II;
                            end;
                            Inc(II);
                        end;
                        ROWMAX := Max(ROWMAX, ABSReal(A[JMAX,IMAX]));
                    end;
                    VV := COLMAX/ROWMAX;
                    if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX*VV) then
                    begin
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        KP := K;
                    end
                    else
                    begin
                        if AP_FP_Greater_Eq(ABSReal(A[IMAX,IMAX]),ALPHA*ROWMAX) then
                        begin
                            
                            //
                            // interchange rows and columns K and IMAX, use 1-by-1
                            // pivot block
                            //
                            KP := IMAX;
                        end
                        else
                        begin
                            
                            //
                            // interchange rows and columns K-1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            KP := IMAX;
                            KSTEP := 2;
                        end;
                    end;
                end;
                KK := K-KSTEP+1;
                if KP<>KK then
                begin
                    
                    //
                    // Interchange rows and columns KK and KP in the leading
                    // submatrix A(1:k,1:k)
                    //
                    I1 := KP-1;
                    for i_ := 1 to I1 do
                    begin
                        Temp[i_] := A[i_,KK];
                    end;
                    for i_ := 1 to I1 do
                    begin
                        A[i_,KK] := A[i_,KP];
                    end;
                    for i_ := 1 to I1 do
                    begin
                        A[i_,KP] := Temp[i_];
                    end;
                    I1 := KP+1;
                    I2 := KK-1;
                    for i_ := I1 to I2 do
                    begin
                        Temp[i_] := A[i_,KK];
                    end;
                    for i_ := I1 to I2 do
                    begin
                        A[i_,KK] := A[KP,i_];
                    end;
                    APVMove(@A[KP][0], I1, I2, @Temp[0], I1, I2);
                    T := A[KK,KK];
                    A[KK,KK] := A[KP,KP];
                    A[KP,KP] := T;
                    if KSTEP=2 then
                    begin
                        T := A[K-1,K];
                        A[K-1,K] := A[KP,K];
                        A[KP,K] := T;
                    end;
                end;
                
                //
                // Update the leading submatrix
                //
                if KSTEP=1 then
                begin
                    
                    //
                    // 1-by-1 pivot block D(k): column k now holds
                    //
                    // W(k) = U(k)*D(k)
                    //
                    // where U(k) is the k-th column of U
                    //
                    // Perform a rank-1 update of A(1:k-1,1:k-1) as
                    //
                    // A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
                    //
                    R1 := 1/A[K,K];
                    I:=1;
                    while I<=K-1 do
                    begin
                        I2 := K-1;
                        VV := -R1*A[I,K];
                        for i_ := I to I2 do
                        begin
                            A[I,i_] := A[I,i_] + VV*A[i_,K];
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Store U(k) in column k
                    //
                    I2 := K-1;
                    for i_ := 1 to I2 do
                    begin
                        A[i_,K] := R1*A[i_,K];
                    end;
                end
                else
                begin
                    
                    //
                    // 2-by-2 pivot block D(k): columns k and k-1 now hold
                    //
                    // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                    //
                    // where U(k) and U(k-1) are the k-th and (k-1)-th columns
                    // of U
                    //
                    // Perform a rank-2 update of A(1:k-2,1:k-2) as
                    //
                    // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
                    //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
                    //
                    if K>2 then
                    begin
                        D12 := A[K-1,K];
                        D22 := A[K-1,K-1]/D12;
                        D11 := A[K,K]/D12;
                        T := 1/(D11*D22-1);
                        D12 := T/D12;
                        J:=K-2;
                        while J>=1 do
                        begin
                            WKM1 := D12*(D11*A[J,K-1]-A[J,K]);
                            WK := D12*(D22*A[J,K]-A[J,K-1]);
                            for i_ := 1 to J do
                            begin
                                A[i_,J] := A[i_,J] - WK*A[i_,K];
                            end;
                            I1 := K-1;
                            for i_ := 1 to J do
                            begin
                                A[i_,J] := A[i_,J] - WKM1*A[i_,I1];
                            end;
                            A[J,K] := WK;
                            A[J,K-1] := WKM1;
                            Dec(J);
                        end;
                    end;
                end;
            end;
            
            //
            // Store details of the interchanges in IPIV
            //
            if KSTEP=1 then
            begin
                Pivots[K] := KP;
            end
            else
            begin
                Pivots[K] := -KP;
                Pivots[K-1] := -KP;
            end;
            
            //
            // Decrease K and return to the start of the main loop
            //
            K := K-KSTEP;
        end;
    end
    else
    begin
        
        //
        // Factorize A as L*D*L' using the lower triangle of A
        //
        // K is the main loop index, increasing from 1 to N in steps of
        // 1 or 2
        //
        K := 1;
        while K<=N do
        begin
            KSTEP := 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            ABSAKK := ABSReal(A[K,K]);
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K, and COLMAX is its absolute value
            //
            if K<N then
            begin
                IMAX := K+1;
                II:=K+2;
                while II<=N do
                begin
                    if AP_FP_Greater(AbsReal(A[II,K]),AbsReal(A[IMAX,K])) then
                    begin
                        IMAX := II;
                    end;
                    Inc(II);
                end;
                COLMAX := ABSReal(A[IMAX,K]);
            end
            else
            begin
                COLMAX := 0;
            end;
            if AP_FP_Eq(Max(ABSAKK, COLMAX),0) then
            begin
                
                //
                // Column K is zero
                //
                KP := K;
            end
            else
            begin
                if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX) then
                begin
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    KP := K;
                end
                else
                begin
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    JMAX := K;
                    II:=K+1;
                    while II<=IMAX-1 do
                    begin
                        if AP_FP_Greater(AbsReal(A[IMAX,II]),AbsReal(A[IMAX,JMAX])) then
                        begin
                            JMAX := II;
                        end;
                        Inc(II);
                    end;
                    ROWMAX := ABSReal(A[IMAX,JMAX]);
                    if IMAX<N then
                    begin
                        JMAX := IMAX+1;
                        II:=IMAX+2;
                        while II<=N do
                        begin
                            if AP_FP_Greater(AbsReal(A[II,IMAX]),AbsReal(A[JMAX,IMAX])) then
                            begin
                                JMAX := II;
                            end;
                            Inc(II);
                        end;
                        ROWMAX := Max(ROWMAX, ABSReal(A[JMAX,IMAX]));
                    end;
                    VV := COLMAX/ROWMAX;
                    if AP_FP_Greater_Eq(ABSAKK,ALPHA*COLMAX*VV) then
                    begin
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        KP := K;
                    end
                    else
                    begin
                        if AP_FP_Greater_Eq(ABSReal(A[IMAX,IMAX]),ALPHA*ROWMAX) then
                        begin
                            
                            //
                            // interchange rows and columns K and IMAX, use 1-by-1
                            // pivot block
                            //
                            KP := IMAX;
                        end
                        else
                        begin
                            
                            //
                            // interchange rows and columns K+1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            KP := IMAX;
                            KSTEP := 2;
                        end;
                    end;
                end;
                KK := K+KSTEP-1;
                if KP<>KK then
                begin
                    
                    //
                    //              Interchange rows and columns KK and KP in the trailing
                    //              submatrix A(k:n,k:n)
                    //
                    if KP<N then
                    begin
                        I1 := KP+1;
                        for i_ := I1 to N do
                        begin
                            Temp[i_] := A[i_,KK];
                        end;
                        for i_ := I1 to N do
                        begin
                            A[i_,KK] := A[i_,KP];
                        end;
                        for i_ := I1 to N do
                        begin
                            A[i_,KP] := Temp[i_];
                        end;
                    end;
                    I1 := KK+1;
                    I2 := KP-1;
                    for i_ := I1 to I2 do
                    begin
                        Temp[i_] := A[i_,KK];
                    end;
                    for i_ := I1 to I2 do
                    begin
                        A[i_,KK] := A[KP,i_];
                    end;
                    APVMove(@A[KP][0], I1, I2, @Temp[0], I1, I2);
                    T := A[KK,KK];
                    A[KK,KK] := A[KP,KP];
                    A[KP,KP] := T;
                    if KSTEP=2 then
                    begin
                        T := A[K+1,K];
                        A[K+1,K] := A[KP,K];
                        A[KP,K] := T;
                    end;
                end;
                
                //
                // Update the trailing submatrix
                //
                if KSTEP=1 then
                begin
                    
                    //
                    // 1-by-1 pivot block D(k): column k now holds
                    //
                    // W(k) = L(k)*D(k)
                    //
                    // where L(k) is the k-th column of L
                    //
                    if K<N then
                    begin
                        
                        //
                        // Perform a rank-1 update of A(k+1:n,k+1:n) as
                        //
                        // A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
                        //
                        D11 := 1/A[K,K];
                        II:=K+1;
                        while II<=N do
                        begin
                            I1 := K+1;
                            I2 := II;
                            VV := -D11*A[II,K];
                            for i_ := I1 to I2 do
                            begin
                                A[II,i_] := A[II,i_] + VV*A[i_,K];
                            end;
                            Inc(II);
                        end;
                        
                        //
                        // Store L(k) in column K
                        //
                        I1 := K+1;
                        for i_ := I1 to N do
                        begin
                            A[i_,K] := D11*A[i_,K];
                        end;
                    end;
                end
                else
                begin
                    
                    //
                    // 2-by-2 pivot block D(k)
                    //
                    if K<N-1 then
                    begin
                        
                        //
                        // Perform a rank-2 update of A(k+2:n,k+2:n) as
                        //
                        // A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
                        //
                        // where L(k) and L(k+1) are the k-th and (k+1)-th
                        // columns of L
                        //
                        D21 := A[K+1,K];
                        D11 := A[K+1,K+1]/D21;
                        D22 := A[K,K]/D21;
                        T := 1/(D11*D22-1);
                        D21 := T/D21;
                        J:=K+2;
                        while J<=N do
                        begin
                            WK := D21*(D11*A[J,K]-A[J,K+1]);
                            WKP1 := D21*(D22*A[J,K+1]-A[J,K]);
                            II := K+1;
                            for i_ := J to N do
                            begin
                                A[i_,J] := A[i_,J] - WK*A[i_,K];
                            end;
                            for i_ := J to N do
                            begin
                                A[i_,J] := A[i_,J] - WKP1*A[i_,II];
                            end;
                            A[J,K] := WK;
                            A[J,K+1] := WKP1;
                            Inc(J);
                        end;
                    end;
                end;
            end;
            
            //
            // Store details of the interchanges in IPIV
            //
            if KSTEP=1 then
            begin
                Pivots[K] := KP;
            end
            else
            begin
                Pivots[K] := -KP;
                Pivots[K+1] := -KP;
            end;
            
            //
            // Increase K and return to the start of the main loop
            //
            K := K+KSTEP;
        end;
    end;
end;


end.