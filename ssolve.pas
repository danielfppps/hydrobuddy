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
unit ssolve;
interface
uses Math, Sysutils, Ap, ldlt;

function SMatrixLDLTSolve(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
function SMatrixSolve(A : TReal2DArray;
     const B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
function SolveSystemLDLT(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
function SolveSymmetricSystem(A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;

implementation

(*************************************************************************
Solving  a system  of linear equations  with a system matrix  given by its
LDLT decomposition

The algorithm solves systems with a square matrix only.

Input parameters:
    A       -   LDLT decomposition of the matrix (the result of the
                SMatrixLDLT subroutine).
    Pivots  -   row permutation table (the result of the SMatrixLDLT subroutine).
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   points to the triangle of matrix A in which the LDLT
                decomposition is stored.
                If IsUpper=True, the decomposition has the form of U*D*U',
                matrix U is stored in the upper triangle of  matrix A  (in
                that case, the lower triangle isn't used and isn't changed
                by the subroutine).
                Similarly, if IsUpper=False, the decomposition has the form
                of L*D*L' and the lower triangle stores matrix L.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular. X contains the solution.
    False, if the matrix is singular (the determinant of matrix D is equal
to 0). In this case, X doesn't contain a solution.
*************************************************************************)
function SMatrixLDLTSolve(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    K : AlglibInteger;
    KP : AlglibInteger;
    AK : Double;
    AKM1 : Double;
    AKM1K : Double;
    BK : Double;
    BKM1 : Double;
    DENOM : Double;
    V : Double;
    i_ : AlglibInteger;
begin
    B := DynamicArrayCopy(B);
    
    //
    // Quick return if possible
    //
    Result := True;
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // Check that the diagonal matrix D is nonsingular
    //
    if IsUpper then
    begin
        
        //
        // Upper triangular storage: examine D from bottom to top
        //
        I:=N-1;
        while I>=0 do
        begin
            if (Pivots[I]>=0) and AP_FP_Eq(A[I,I],0) then
            begin
                Result := False;
                Exit;
            end;
            Dec(I);
        end;
    end
    else
    begin
        
        //
        // Lower triangular storage: examine D from top to bottom.
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
    end;
    
    //
    // Solve Ax = b
    //
    if IsUpper then
    begin
        
        //
        // Solve A*X = B, where A = U*D*U'.
        //
        // First solve U*D*X = B, overwriting B with X.
        //
        // K+1 is the main loop index, decreasing from N to 1 in steps of
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
                // Interchange rows K+1 and IPIV(K+1).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(U(K+1)), where U(K+1) is the transformation
                // stored in column K+1 of A.
                //
                V := B[K];
                for i_ := 0 to K-1 do
                begin
                    B[i_] := B[i_] - V*A[i_,K];
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                B[K] := B[K]/A[K,K];
                K := K-1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K+1-1 and -IPIV(K+1).
                //
                KP := Pivots[K]+N;
                if KP<>K-1 then
                begin
                    V := B[K-1];
                    B[K-1] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(U(K+1)), where U(K+1) is the transformation
                // stored in columns K+1-1 and K+1 of A.
                //
                V := B[K];
                for i_ := 0 to K-2 do
                begin
                    B[i_] := B[i_] - V*A[i_,K];
                end;
                V := B[K-1];
                for i_ := 0 to K-2 do
                begin
                    B[i_] := B[i_] - V*A[i_,K-1];
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                AKM1K := A[K-1,K];
                AKM1 := A[K-1,K-1]/AKM1K;
                AK := A[K,K]/AKM1K;
                DENOM := AKM1*AK-1;
                BKM1 := B[K-1]/AKM1K;
                BK := B[K]/AKM1K;
                B[K-1] := (AK*BKM1-BK)/DENOM;
                B[K] := (AKM1*BK-BKM1)/DENOM;
                K := K-2;
            end;
        end;
        
        //
        // Next solve U'*X = B, overwriting B with X.
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
                // Multiply by inv(U'(K+1)), where U(K+1) is the transformation
                // stored in column K+1 of A.
                //
                V := 0.0;
                for i_ := 0 to K-1 do
                begin
                    V := V + B[i_]*A[i_,K];
                end;
                B[K] := B[K]-V;
                
                //
                // Interchange rows K+1 and IPIV(K+1).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K+1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(U'(K+1+1)), where U(K+1+1) is the transformation
                // stored in columns K+1 and K+1+1 of A.
                //
                V := 0.0;
                for i_ := 0 to K-1 do
                begin
                    V := V + B[i_]*A[i_,K];
                end;
                B[K] := B[K]-V;
                V := 0.0;
                for i_ := 0 to K-1 do
                begin
                    V := V + B[i_]*A[i_,K+1];
                end;
                B[K+1] := B[K+1]-V;
                
                //
                // Interchange rows K+1 and -IPIV(K+1).
                //
                KP := Pivots[K]+N;
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K+2;
            end;
        end;
    end
    else
    begin
        
        //
        // Solve A*X = B, where A = L*D*L'.
        //
        // First solve L*D*X = B, overwriting B with X.
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
                // Interchange rows K+1 and IPIV(K+1).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(L(K+1)), where L(K+1) is the transformation
                // stored in column K+1 of A.
                //
                if K+1<N then
                begin
                    V := B[K];
                    for i_ := K+1 to N-1 do
                    begin
                        B[i_] := B[i_] - V*A[i_,K];
                    end;
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                B[K] := B[K]/A[K,K];
                K := K+1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K+1+1 and -IPIV(K+1).
                //
                KP := Pivots[K]+N;
                if KP<>K+1 then
                begin
                    V := B[K+1];
                    B[K+1] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(L(K+1)), where L(K+1) is the transformation
                // stored in columns K+1 and K+1+1 of A.
                //
                if K+1<N-1 then
                begin
                    V := B[K];
                    for i_ := K+2 to N-1 do
                    begin
                        B[i_] := B[i_] - V*A[i_,K];
                    end;
                    V := B[K+1];
                    for i_ := K+2 to N-1 do
                    begin
                        B[i_] := B[i_] - V*A[i_,K+1];
                    end;
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                AKM1K := A[K+1,K];
                AKM1 := A[K,K]/AKM1K;
                AK := A[K+1,K+1]/AKM1K;
                DENOM := AKM1*AK-1;
                BKM1 := B[K]/AKM1K;
                BK := B[K+1]/AKM1K;
                B[K] := (AK*BKM1-BK)/DENOM;
                B[K+1] := (AKM1*BK-BKM1)/DENOM;
                K := K+2;
            end;
        end;
        
        //
        // Next solve L'*X = B, overwriting B with X.
        //
        // K+1 is the main loop index, decreasing from N to 1 in steps of
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
                // Multiply by inv(L'(K+1)), where L(K+1) is the transformation
                // stored in column K+1 of A.
                //
                if K+1<N then
                begin
                    V := 0.0;
                    for i_ := K+1 to N-1 do
                    begin
                        V := V + B[i_]*A[i_,K];
                    end;
                    B[K] := B[K]-V;
                end;
                
                //
                // Interchange rows K+1 and IPIV(K+1).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K-1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(L'(K+1-1)), where L(K+1-1) is the transformation
                // stored in columns K+1-1 and K+1 of A.
                //
                if K+1<N then
                begin
                    V := 0.0;
                    for i_ := K+1 to N-1 do
                    begin
                        V := V + B[i_]*A[i_,K];
                    end;
                    B[K] := B[K]-V;
                    V := 0.0;
                    for i_ := K+1 to N-1 do
                    begin
                        V := V + B[i_]*A[i_,K-1];
                    end;
                    B[K-1] := B[K-1]-V;
                end;
                
                //
                // Interchange rows K+1 and -IPIV(K+1).
                //
                KP := Pivots[K]+N;
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K-2;
            end;
        end;
    end;
    SetLength(X, N-1+1);
    APVMove(@X[0], 0, N-1, @B[0], 0, N-1);
end;


(*************************************************************************
Solving a system of linear equations with a symmetric system matrix

Input parameters:
    A       -   system matrix (upper or lower triangle).
                Array whose indexes range within [0..N-1, 0..N-1].
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   If IsUpper = True, A contains the upper triangle,
                otherwise A contains the lower triangle.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular. X contains the solution.
    False, if the matrix is singular (the determinant of the matrix is equal
to 0). In this case, X doesn't contain a solution.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
function SMatrixSolve(A : TReal2DArray;
     const B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
var
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    SMatrixLDLT(A, N, IsUpper, Pivots);
    Result := SMatrixLDLTSolve(A, Pivots, B, N, IsUpper, X);
end;


function SolveSystemLDLT(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    K : AlglibInteger;
    KP : AlglibInteger;
    KM1 : AlglibInteger;
    KM2 : AlglibInteger;
    KP1 : AlglibInteger;
    KP2 : AlglibInteger;
    AK : Double;
    AKM1 : Double;
    AKM1K : Double;
    BK : Double;
    BKM1 : Double;
    DENOM : Double;
    V : Double;
    i_ : AlglibInteger;
begin
    B := DynamicArrayCopy(B);
    
    //
    // Quick return if possible
    //
    Result := True;
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // Check that the diagonal matrix D is nonsingular
    //
    if IsUpper then
    begin
        
        //
        // Upper triangular storage: examine D from bottom to top
        //
        I:=N;
        while I>=1 do
        begin
            if (Pivots[I]>0) and AP_FP_Eq(A[I,I],0) then
            begin
                Result := False;
                Exit;
            end;
            Dec(I);
        end;
    end
    else
    begin
        
        //
        // Lower triangular storage: examine D from top to bottom.
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
    end;
    
    //
    // Solve Ax = b
    //
    if IsUpper then
    begin
        
        //
        // Solve A*X = B, where A = U*D*U'.
        //
        // First solve U*D*X = B, overwriting B with X.
        //
        // K is the main loop index, decreasing from N to 1 in steps of
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
                // Interchange rows K and IPIV(K).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(U(K)), where U(K) is the transformation
                // stored in column K of A.
                //
                KM1 := K-1;
                V := B[K];
                for i_ := 1 to KM1 do
                begin
                    B[i_] := B[i_] - V*A[i_,K];
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                B[K] := B[K]/A[K,K];
                K := K-1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K-1 and -IPIV(K).
                //
                KP := -Pivots[K];
                if KP<>K-1 then
                begin
                    V := B[K-1];
                    B[K-1] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(U(K)), where U(K) is the transformation
                // stored in columns K-1 and K of A.
                //
                KM2 := K-2;
                KM1 := K-1;
                V := B[K];
                for i_ := 1 to KM2 do
                begin
                    B[i_] := B[i_] - V*A[i_,K];
                end;
                V := B[K-1];
                for i_ := 1 to KM2 do
                begin
                    B[i_] := B[i_] - V*A[i_,KM1];
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                AKM1K := A[K-1,K];
                AKM1 := A[K-1,K-1]/AKM1K;
                AK := A[K,K]/AKM1K;
                DENOM := AKM1*AK-1;
                BKM1 := B[K-1]/AKM1K;
                BK := B[K]/AKM1K;
                B[K-1] := (AK*BKM1-BK)/DENOM;
                B[K] := (AKM1*BK-BKM1)/DENOM;
                K := K-2;
            end;
        end;
        
        //
        // Next solve U'*X = B, overwriting B with X.
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
                // Multiply by inv(U'(K)), where U(K) is the transformation
                // stored in column K of A.
                //
                KM1 := K-1;
                V := 0.0;
                for i_ := 1 to KM1 do
                begin
                    V := V + B[i_]*A[i_,K];
                end;
                B[K] := B[K]-V;
                
                //
                // Interchange rows K and IPIV(K).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K+1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(U'(K+1)), where U(K+1) is the transformation
                // stored in columns K and K+1 of A.
                //
                KM1 := K-1;
                KP1 := K+1;
                V := 0.0;
                for i_ := 1 to KM1 do
                begin
                    V := V + B[i_]*A[i_,K];
                end;
                B[K] := B[K]-V;
                V := 0.0;
                for i_ := 1 to KM1 do
                begin
                    V := V + B[i_]*A[i_,KP1];
                end;
                B[K+1] := B[K+1]-V;
                
                //
                // Interchange rows K and -IPIV(K).
                //
                KP := -Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K+2;
            end;
        end;
    end
    else
    begin
        
        //
        // Solve A*X = B, where A = L*D*L'.
        //
        // First solve L*D*X = B, overwriting B with X.
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
                // Interchange rows K and IPIV(K).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(L(K)), where L(K) is the transformation
                // stored in column K of A.
                //
                if K<N then
                begin
                    KP1 := K+1;
                    V := B[K];
                    for i_ := KP1 to N do
                    begin
                        B[i_] := B[i_] - V*A[i_,K];
                    end;
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                B[K] := B[K]/A[K,K];
                K := K+1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K+1 and -IPIV(K).
                //
                KP := -Pivots[K];
                if KP<>K+1 then
                begin
                    V := B[K+1];
                    B[K+1] := B[KP];
                    B[KP] := V;
                end;
                
                //
                // Multiply by inv(L(K)), where L(K) is the transformation
                // stored in columns K and K+1 of A.
                //
                if K<N-1 then
                begin
                    KP1 := K+1;
                    KP2 := K+2;
                    V := B[K];
                    for i_ := KP2 to N do
                    begin
                        B[i_] := B[i_] - V*A[i_,K];
                    end;
                    V := B[K+1];
                    for i_ := KP2 to N do
                    begin
                        B[i_] := B[i_] - V*A[i_,KP1];
                    end;
                end;
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                AKM1K := A[K+1,K];
                AKM1 := A[K,K]/AKM1K;
                AK := A[K+1,K+1]/AKM1K;
                DENOM := AKM1*AK-1;
                BKM1 := B[K]/AKM1K;
                BK := B[K+1]/AKM1K;
                B[K] := (AK*BKM1-BK)/DENOM;
                B[K+1] := (AKM1*BK-BKM1)/DENOM;
                K := K+2;
            end;
        end;
        
        //
        // Next solve L'*X = B, overwriting B with X.
        //
        // K is the main loop index, decreasing from N to 1 in steps of
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
                // Multiply by inv(L'(K)), where L(K) is the transformation
                // stored in column K of A.
                //
                if K<N then
                begin
                    KP1 := K+1;
                    V := 0.0;
                    for i_ := KP1 to N do
                    begin
                        V := V + B[i_]*A[i_,K];
                    end;
                    B[K] := B[K]-V;
                end;
                
                //
                // Interchange rows K and IPIV(K).
                //
                KP := Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K-1;
            end
            else
            begin
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(L'(K-1)), where L(K-1) is the transformation
                // stored in columns K-1 and K of A.
                //
                if K<N then
                begin
                    KP1 := K+1;
                    KM1 := K-1;
                    V := 0.0;
                    for i_ := KP1 to N do
                    begin
                        V := V + B[i_]*A[i_,K];
                    end;
                    B[K] := B[K]-V;
                    V := 0.0;
                    for i_ := KP1 to N do
                    begin
                        V := V + B[i_]*A[i_,KM1];
                    end;
                    B[K-1] := B[K-1]-V;
                end;
                
                //
                // Interchange rows K and -IPIV(K).
                //
                KP := -Pivots[K];
                if KP<>K then
                begin
                    V := B[K];
                    B[K] := B[KP];
                    B[KP] := V;
                end;
                K := K-2;
            end;
        end;
    end;
    SetLength(X, N+1);
    APVMove(@X[0], 1, N, @B[0], 1, N);
end;


function SolveSymmetricSystem(A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
var
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    LDLTDecomposition(A, N, IsUpper, Pivots);
    Result := SolveSystemLDLT(A, Pivots, B, N, IsUpper, X);
end;


end.
