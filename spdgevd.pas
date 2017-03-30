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
unit spdgevd;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, sblas, blas, trlinsolve, safesolve, rcond, matinv, hblas, ortfac, rotations, hsschur, evd;

function SMatrixGEVD(A : TReal2DArray;
     N : AlglibInteger;
     IsUpperA : Boolean;
     const B : TReal2DArray;
     IsUpperB : Boolean;
     ZNeeded : AlglibInteger;
     ProblemType : AlglibInteger;
     var D : TReal1DArray;
     var Z : TReal2DArray):Boolean;
function SMatrixGEVDReduce(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpperA : Boolean;
     const B : TReal2DArray;
     IsUpperB : Boolean;
     ProblemType : AlglibInteger;
     var R : TReal2DArray;
     var IsUpperR : Boolean):Boolean;

implementation

(*************************************************************************
Algorithm for solving the following generalized symmetric positive-definite
eigenproblem:
    A*x = lambda*B*x (1) or
    A*B*x = lambda*x (2) or
    B*A*x = lambda*x (3).
where A is a symmetric matrix, B - symmetric positive-definite matrix.
The problem is solved by reducing it to an ordinary  symmetric  eigenvalue
problem.

Input parameters:
    A           -   symmetric matrix which is given by its upper or lower
                    triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
    N           -   size of matrices A and B.
    IsUpperA    -   storage format of matrix A.
    B           -   symmetric positive-definite matrix which is given by
                    its upper or lower triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
    IsUpperB    -   storage format of matrix B.
    ZNeeded     -   if ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.
    ProblemType -   if ProblemType is equal to:
                     * 1, the following problem is solved: A*x = lambda*B*x;
                     * 2, the following problem is solved: A*B*x = lambda*x;
                     * 3, the following problem is solved: B*A*x = lambda*x.

Output parameters:
    D           -   eigenvalues in ascending order.
                    Array whose index ranges within [0..N-1].
    Z           -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains eigenvectors.
                    Array whose indexes range within [0..N-1, 0..N-1].
                    The eigenvectors are stored in matrix columns. It should
                    be noted that the eigenvectors in such problems do not
                    form an orthogonal system.

Result:
    True, if the problem was solved successfully.
    False, if the error occurred during the Cholesky decomposition of matrix
    B (the matrix isn’t positive-definite) or during the work of the iterative
    algorithm for solving the symmetric eigenproblem.

See also the GeneralizedSymmetricDefiniteEVDReduce subroutine.

  -- ALGLIB --
     Copyright 1.28.2006 by Bochkanov Sergey
*************************************************************************)
function SMatrixGEVD(A : TReal2DArray;
     N : AlglibInteger;
     IsUpperA : Boolean;
     const B : TReal2DArray;
     IsUpperB : Boolean;
     ZNeeded : AlglibInteger;
     ProblemType : AlglibInteger;
     var D : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    R : TReal2DArray;
    T : TReal2DArray;
    IsUpperR : Boolean;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    J1INC : AlglibInteger;
    J2INC : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    A := DynamicArrayCopy(A);
    
    //
    // Reduce and solve
    //
    Result := SMatrixGEVDReduce(A, N, IsUpperA, B, IsUpperB, ProblemType, R, IsUpperR);
    if  not Result then
    begin
        Exit;
    end;
    Result := SMatrixEVD(A, N, ZNeeded, IsUpperA, D, T);
    if  not Result then
    begin
        Exit;
    end;
    
    //
    // Transform eigenvectors if needed
    //
    if ZNeeded<>0 then
    begin
        
        //
        // fill Z with zeros
        //
        SetLength(Z, N-1+1, N-1+1);
        J:=0;
        while J<=N-1 do
        begin
            Z[0,J] := Double(0.0);
            Inc(J);
        end;
        I:=1;
        while I<=N-1 do
        begin
            APVMove(@Z[I][0], 0, N-1, @Z[0][0], 0, N-1);
            Inc(I);
        end;
        
        //
        // Setup R properties
        //
        if IsUpperR then
        begin
            J1 := 0;
            J2 := N-1;
            J1INC := +1;
            J2INC := 0;
        end
        else
        begin
            J1 := 0;
            J2 := 0;
            J1INC := 0;
            J2INC := +1;
        end;
        
        //
        // Calculate R*Z
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=J1;
            while J<=J2 do
            begin
                V := R[I,J];
                APVAdd(@Z[I][0], 0, N-1, @T[J][0], 0, N-1, V);
                Inc(J);
            end;
            J1 := J1+J1INC;
            J2 := J2+J2INC;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Algorithm for reduction of the following generalized symmetric positive-
definite eigenvalue problem:
    A*x = lambda*B*x (1) or
    A*B*x = lambda*x (2) or
    B*A*x = lambda*x (3)
to the symmetric eigenvalues problem C*y = lambda*y (eigenvalues of this and
the given problems are the same, and the eigenvectors of the given problem
could be obtained by multiplying the obtained eigenvectors by the
transformation matrix x = R*y).

Here A is a symmetric matrix, B - symmetric positive-definite matrix.

Input parameters:
    A           -   symmetric matrix which is given by its upper or lower
                    triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
    N           -   size of matrices A and B.
    IsUpperA    -   storage format of matrix A.
    B           -   symmetric positive-definite matrix which is given by
                    its upper or lower triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
    IsUpperB    -   storage format of matrix B.
    ProblemType -   if ProblemType is equal to:
                     * 1, the following problem is solved: A*x = lambda*B*x;
                     * 2, the following problem is solved: A*B*x = lambda*x;
                     * 3, the following problem is solved: B*A*x = lambda*x.

Output parameters:
    A           -   symmetric matrix which is given by its upper or lower
                    triangle depending on IsUpperA. Contains matrix C.
                    Array whose indexes range within [0..N-1, 0..N-1].
    R           -   upper triangular or low triangular transformation matrix
                    which is used to obtain the eigenvectors of a given problem
                    as the product of eigenvectors of C (from the right) and
                    matrix R (from the left). If the matrix is upper
                    triangular, the elements below the main diagonal
                    are equal to 0 (and vice versa). Thus, we can perform
                    the multiplication without taking into account the
                    internal structure (which is an easier though less
                    effective way).
                    Array whose indexes range within [0..N-1, 0..N-1].
    IsUpperR    -   type of matrix R (upper or lower triangular).

Result:
    True, if the problem was reduced successfully.
    False, if the error occurred during the Cholesky decomposition of
        matrix B (the matrix is not positive-definite).

  -- ALGLIB --
     Copyright 1.28.2006 by Bochkanov Sergey
*************************************************************************)
function SMatrixGEVDReduce(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpperA : Boolean;
     const B : TReal2DArray;
     IsUpperB : Boolean;
     ProblemType : AlglibInteger;
     var R : TReal2DArray;
     var IsUpperR : Boolean):Boolean;
var
    T : TReal2DArray;
    W1 : TReal1DArray;
    W2 : TReal1DArray;
    W3 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Rep : MatInvReport;
    Info : AlglibInteger;
    i_ : AlglibInteger;
begin
    Assert(N>0, 'SMatrixGEVDReduce: N<=0!');
    Assert((ProblemType=1) or (ProblemType=2) or (ProblemType=3), 'SMatrixGEVDReduce: incorrect ProblemType!');
    Result := True;
    
    //
    // Problem 1:  A*x = lambda*B*x
    //
    // Reducing to:
    //     C*y = lambda*y
    //     C = L^(-1) * A * L^(-T)
    //     x = L^(-T) * y
    //
    if ProblemType=1 then
    begin
        
        //
        // Factorize B in T: B = LL'
        //
        SetLength(T, N-1+1, N-1+1);
        if IsUpperB then
        begin
            I:=0;
            while I<=N-1 do
            begin
                for i_ := I to N-1 do
                begin
                    T[i_,I] := B[I,i_];
                end;
                Inc(I);
            end;
        end
        else
        begin
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@T[I][0], 0, I, @B[I][0], 0, I);
                Inc(I);
            end;
        end;
        if  not SPDMatrixCholesky(T, N, False) then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Invert L in T
        //
        RMatrixTRInverse(T, N, False, False, Info, Rep);
        if Info<=0 then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Build L^(-1) * A * L^(-T) in R
        //
        SetLength(W1, N+1);
        SetLength(W2, N+1);
        SetLength(R, N-1+1, N-1+1);
        J:=1;
        while J<=N do
        begin
            
            //
            // Form w2 = A * l'(j) (here l'(j) is j-th column of L^(-T))
            //
            APVMove(@W1[0], 1, J, @T[J-1][0], 0, J-1);
            SymmetricMatrixVectorMultiply(A, IsUpperA, 0, J-1, W1, Double(1.0), W2);
            if IsUpperA then
            begin
                MatrixVectorMultiply(A, 0, J-1, J, N-1, True, W1, 1, J, Double(1.0), W2, J+1, N, Double(0.0));
            end
            else
            begin
                MatrixVectorMultiply(A, J, N-1, 0, J-1, False, W1, 1, J, Double(1.0), W2, J+1, N, Double(0.0));
            end;
            
            //
            // Form l(i)*w2 (here l(i) is i-th row of L^(-1))
            //
            I:=1;
            while I<=N do
            begin
                V := APVDotProduct(@T[I-1][0], 0, I-1, @W2[0], 1, I);
                R[I-1,J-1] := V;
                Inc(I);
            end;
            Inc(J);
        end;
        
        //
        // Copy R to A
        //
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@A[I][0], 0, N-1, @R[I][0], 0, N-1);
            Inc(I);
        end;
        
        //
        // Copy L^(-1) from T to R and transpose
        //
        IsUpperR := True;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=I-1 do
            begin
                R[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            for i_ := I to N-1 do
            begin
                R[I,i_] := T[i_,I];
            end;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Problem 2:  A*B*x = lambda*x
    // or
    // problem 3:  B*A*x = lambda*x
    //
    // Reducing to:
    //     C*y = lambda*y
    //     C = U * A * U'
    //     B = U'* U
    //
    if (ProblemType=2) or (ProblemType=3) then
    begin
        
        //
        // Factorize B in T: B = U'*U
        //
        SetLength(T, N-1+1, N-1+1);
        if IsUpperB then
        begin
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@T[I][0], I, N-1, @B[I][0], I, N-1);
                Inc(I);
            end;
        end
        else
        begin
            I:=0;
            while I<=N-1 do
            begin
                for i_ := I to N-1 do
                begin
                    T[I,i_] := B[i_,I];
                end;
                Inc(I);
            end;
        end;
        if  not SPDMatrixCholesky(T, N, True) then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Build U * A * U' in R
        //
        SetLength(W1, N+1);
        SetLength(W2, N+1);
        SetLength(W3, N+1);
        SetLength(R, N-1+1, N-1+1);
        J:=1;
        while J<=N do
        begin
            
            //
            // Form w2 = A * u'(j) (here u'(j) is j-th column of U')
            //
            APVMove(@W1[0], 1, N-J+1, @T[J-1][0], J-1, N-1);
            SymmetricMatrixVectorMultiply(A, IsUpperA, J-1, N-1, W1, Double(1.0), W3);
            APVMove(@W2[0], J, N, @W3[0], 1, N-J+1);
            APVMove(@W1[0], J, N, @T[J-1][0], J-1, N-1);
            if IsUpperA then
            begin
                MatrixVectorMultiply(A, 0, J-2, J-1, N-1, False, W1, J, N, Double(1.0), W2, 1, J-1, Double(0.0));
            end
            else
            begin
                MatrixVectorMultiply(A, J-1, N-1, 0, J-2, True, W1, J, N, Double(1.0), W2, 1, J-1, Double(0.0));
            end;
            
            //
            // Form u(i)*w2 (here u(i) is i-th row of U)
            //
            I:=1;
            while I<=N do
            begin
                V := APVDotProduct(@T[I-1][0], I-1, N-1, @W2[0], I, N);
                R[I-1,J-1] := V;
                Inc(I);
            end;
            Inc(J);
        end;
        
        //
        // Copy R to A
        //
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@A[I][0], 0, N-1, @R[I][0], 0, N-1);
            Inc(I);
        end;
        if ProblemType=2 then
        begin
            
            //
            // Invert U in T
            //
            RMatrixTRInverse(T, N, True, False, Info, Rep);
            if Info<=0 then
            begin
                Result := False;
                Exit;
            end;
            
            //
            // Copy U^-1 from T to R
            //
            IsUpperR := True;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=I-1 do
                begin
                    R[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@R[I][0], I, N-1, @T[I][0], I, N-1);
                Inc(I);
            end;
        end
        else
        begin
            
            //
            // Copy U from T to R and transpose
            //
            IsUpperR := False;
            I:=0;
            while I<=N-1 do
            begin
                J:=I+1;
                while J<=N-1 do
                begin
                    R[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                for i_ := I to N-1 do
                begin
                    R[i_,I] := T[I,i_];
                end;
                Inc(I);
            end;
        end;
    end;
end;


end.