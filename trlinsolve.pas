{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
This file is a part of ALGLIB project.

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
unit trlinsolve;
interface
uses Math, Sysutils, Ap;

procedure RMatrixTRSafeSolve(const A : TReal2DArray;
     N : AlglibInteger;
     var X : TReal1DArray;
     var S : Double;
     IsUpper : Boolean;
     IsTrans : Boolean;
     Isunit : Boolean);
procedure SafeSolveTriangular(const A : TReal2DArray;
     N : AlglibInteger;
     var X : TReal1DArray;
     var S : Double;
     IsUpper : Boolean;
     IsTrans : Boolean;
     Isunit : Boolean;
     NORMIN : Boolean;
     var CNORM : TReal1DArray);

implementation

(*************************************************************************
Utility subroutine performing the "safe" solution of system of linear
equations with triangular coefficient matrices.

The subroutine uses scaling and solves the scaled system A*x=s*b (where  s
is  a  scalar  value)  instead  of  A*x=b,  choosing  s  so  that x can be
represented by a floating-point number. The closer the system  gets  to  a
singular, the less s is. If the system is singular, s=0 and x contains the
non-trivial solution of equation A*x=0.

The feature of an algorithm is that it could not cause an  overflow  or  a
division by zero regardless of the matrix used as the input.

The algorithm can solve systems of equations with  upper/lower  triangular
matrices,  with/without unit diagonal, and systems of type A*x=b or A'*x=b
(where A' is a transposed matrix A).

Input parameters:
    A       -   system matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    X       -   right-hand member of a system.
                Array whose index ranges within [0..N-1].
    IsUpper -   matrix type. If it is True, the system matrix is the upper
                triangular and is located in  the  corresponding  part  of
                matrix A.
    Trans   -   problem type. If it is True, the problem to be  solved  is
                A'*x=b, otherwise it is A*x=b.
    Isunit  -   matrix type. If it is True, the system matrix has  a  unit
                diagonal (the elements on the main diagonal are  not  used
                in the calculation process), otherwise the matrix is considered
                to be a general triangular matrix.

Output parameters:
    X       -   solution. Array whose index ranges within [0..N-1].
    S       -   scaling factor.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************)
procedure RMatrixTRSafeSolve(const A : TReal2DArray;
     N : AlglibInteger;
     var X : TReal1DArray;
     var S : Double;
     IsUpper : Boolean;
     IsTrans : Boolean;
     Isunit : Boolean);
var
    NORMIN : Boolean;
    CNORM : TReal1DArray;
    A1 : TReal2DArray;
    X1 : TReal1DArray;
    I : AlglibInteger;
begin
    
    //
    // From 0-based to 1-based
    //
    NORMIN := False;
    SetLength(A1, N+1, N+1);
    SetLength(X1, N+1);
    I:=1;
    while I<=N do
    begin
        APVMove(@A1[I][0], 1, N, @A[I-1][0], 0, N-1);
        Inc(I);
    end;
    APVMove(@X1[0], 1, N, @X[0], 0, N-1);
    
    //
    // Solve 1-based
    //
    SafeSolveTriangular(A1, N, X1, S, IsUpper, IsTrans, Isunit, NORMIN, CNORM);
    
    //
    // From 1-based to 0-based
    //
    APVMove(@X[0], 0, N-1, @X1[0], 1, N);
end;


(*************************************************************************
Obsolete 1-based subroutine.
See RMatrixTRSafeSolve for 0-based replacement.
*************************************************************************)
procedure SafeSolveTriangular(const A : TReal2DArray;
     N : AlglibInteger;
     var X : TReal1DArray;
     var S : Double;
     IsUpper : Boolean;
     IsTrans : Boolean;
     Isunit : Boolean;
     NORMIN : Boolean;
     var CNORM : TReal1DArray);
var
    I : AlglibInteger;
    IMAX : AlglibInteger;
    J : AlglibInteger;
    JFIRST : AlglibInteger;
    JINC : AlglibInteger;
    JLAST : AlglibInteger;
    JM1 : AlglibInteger;
    JP1 : AlglibInteger;
    IP1 : AlglibInteger;
    IM1 : AlglibInteger;
    K : AlglibInteger;
    Flg : AlglibInteger;
    V : Double;
    VD : Double;
    BIGNUM : Double;
    GROW : Double;
    REC : Double;
    SMLNUM : Double;
    SUMJ : Double;
    TJJ : Double;
    TJJS : Double;
    TMAX : Double;
    TSCAL : Double;
    USCAL : Double;
    XBND : Double;
    XJ : Double;
    XMAX : Double;
    NOTRAN : Boolean;
    UPPER : Boolean;
    NOunit : Boolean;
    i_ : AlglibInteger;
begin
    UPPER := IsUpper;
    NOTRAN :=  not IsTrans;
    NOunit :=  not Isunit;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // Determine machine dependent parameters to control overflow.
    //
    SMLNUM := MinRealNumber/(MachineEpsilon*2);
    BIGNUM := 1/SMLNUM;
    S := 1;
    if  not NORMIN then
    begin
        SetLength(CNORM, N+1);
        
        //
        // Compute the 1-norm of each column, not including the diagonal.
        //
        if UPPER then
        begin
            
            //
            // A is upper triangular.
            //
            J:=1;
            while J<=N do
            begin
                V := 0;
                K:=1;
                while K<=J-1 do
                begin
                    V := V+AbsReal(A[K,J]);
                    Inc(K);
                end;
                CNORM[J] := V;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // A is lower triangular.
            //
            J:=1;
            while J<=N-1 do
            begin
                V := 0;
                K:=J+1;
                while K<=N do
                begin
                    V := V+AbsReal(A[K,J]);
                    Inc(K);
                end;
                CNORM[J] := V;
                Inc(J);
            end;
            CNORM[N] := 0;
        end;
    end;
    
    //
    // Scale the column norms by TSCAL if the maximum element in CNORM is
    // greater than BIGNUM.
    //
    IMAX := 1;
    K:=2;
    while K<=N do
    begin
        if AP_FP_Greater(CNORM[K],CNORM[IMAX]) then
        begin
            IMAX := K;
        end;
        Inc(K);
    end;
    TMAX := CNORM[IMAX];
    if AP_FP_Less_Eq(TMAX,BIGNUM) then
    begin
        TSCAL := 1;
    end
    else
    begin
        TSCAL := 1/(SMLNUM*TMAX);
        APVMul(@CNORM[0], 1, N, TSCAL);
    end;
    
    //
    // Compute a bound on the computed solution vector to see if the
    // Level 2 BLAS routine DTRSV can be used.
    //
    J := 1;
    K:=2;
    while K<=N do
    begin
        if AP_FP_Greater(AbsReal(X[K]),AbsReal(X[J])) then
        begin
            J := K;
        end;
        Inc(K);
    end;
    XMAX := ABSReal(X[J]);
    XBND := XMAX;
    if NOTRAN then
    begin
        
        //
        // Compute the growth in A * x = b.
        //
        if UPPER then
        begin
            JFIRST := N;
            JLAST := 1;
            JINC := -1;
        end
        else
        begin
            JFIRST := 1;
            JLAST := N;
            JINC := 1;
        end;
        if AP_FP_Neq(TSCAL,1) then
        begin
            GROW := 0;
        end
        else
        begin
            if NOunit then
            begin
                
                //
                // A is non-unit triangular.
                //
                // Compute GROW = 1/G(j) and XBND = 1/M(j).
                // Initially, G(0) = max{x(i), i=1,...,n}.
                //
                GROW := 1/Max(XBND, SMLNUM);
                XBND := GROW;
                J := JFIRST;
                while (JINC>0) and (J<=JLAST) or (JINC<0) and (J>=JLAST) do
                begin
                    
                    //
                    // Exit the loop if the growth factor is too small.
                    //
                    if AP_FP_Less_Eq(GROW,SMLNUM) then
                    begin
                        Break;
                    end;
                    
                    //
                    // M(j) = G(j-1) / abs(A(j,j))
                    //
                    TJJ := ABSReal(A[J,J]);
                    XBND := Min(XBND, Min(1, TJJ)*GROW);
                    if AP_FP_Greater_Eq(TJJ+CNORM[J],SMLNUM) then
                    begin
                        
                        //
                        // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                        //
                        GROW := GROW*(TJJ/(TJJ+CNORM[J]));
                    end
                    else
                    begin
                        
                        //
                        // G(j) could overflow, set GROW to 0.
                        //
                        GROW := 0;
                    end;
                    if J=JLAST then
                    begin
                        GROW := XBND;
                    end;
                    J := J+JINC;
                end;
            end
            else
            begin
                
                //
                // A is unit triangular.
                //
                // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                //
                GROW := Min(1, 1/Max(XBND, SMLNUM));
                J := JFIRST;
                while (JINC>0) and (J<=JLAST) or (JINC<0) and (J>=JLAST) do
                begin
                    
                    //
                    // Exit the loop if the growth factor is too small.
                    //
                    if AP_FP_Less_Eq(GROW,SMLNUM) then
                    begin
                        Break;
                    end;
                    
                    //
                    // G(j) = G(j-1)*( 1 + CNORM(j) )
                    //
                    GROW := GROW*(1/(1+CNORM[J]));
                    J := J+JINC;
                end;
            end;
        end;
    end
    else
    begin
        
        //
        // Compute the growth in A' * x = b.
        //
        if UPPER then
        begin
            JFIRST := 1;
            JLAST := N;
            JINC := 1;
        end
        else
        begin
            JFIRST := N;
            JLAST := 1;
            JINC := -1;
        end;
        if AP_FP_Neq(TSCAL,1) then
        begin
            GROW := 0;
        end
        else
        begin
            if NOunit then
            begin
                
                //
                // A is non-unit triangular.
                //
                // Compute GROW = 1/G(j) and XBND = 1/M(j).
                // Initially, M(0) = max{x(i), i=1,...,n}.
                //
                GROW := 1/Max(XBND, SMLNUM);
                XBND := GROW;
                J := JFIRST;
                while (JINC>0) and (J<=JLAST) or (JINC<0) and (J>=JLAST) do
                begin
                    
                    //
                    // Exit the loop if the growth factor is too small.
                    //
                    if AP_FP_Less_Eq(GROW,SMLNUM) then
                    begin
                        Break;
                    end;
                    
                    //
                    // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
                    //
                    XJ := 1+CNORM[J];
                    GROW := Min(GROW, XBND/XJ);
                    
                    //
                    // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
                    //
                    TJJ := ABSReal(A[J,J]);
                    if AP_FP_Greater(XJ,TJJ) then
                    begin
                        XBND := XBND*(TJJ/XJ);
                    end;
                    if J=JLAST then
                    begin
                        GROW := Min(GROW, XBND);
                    end;
                    J := J+JINC;
                end;
            end
            else
            begin
                
                //
                // A is unit triangular.
                //
                // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                //
                GROW := Min(1, 1/Max(XBND, SMLNUM));
                J := JFIRST;
                while (JINC>0) and (J<=JLAST) or (JINC<0) and (J>=JLAST) do
                begin
                    
                    //
                    // Exit the loop if the growth factor is too small.
                    //
                    if AP_FP_Less_Eq(GROW,SMLNUM) then
                    begin
                        Break;
                    end;
                    
                    //
                    // G(j) = ( 1 + CNORM(j) )*G(j-1)
                    //
                    XJ := 1+CNORM[J];
                    GROW := GROW/XJ;
                    J := J+JINC;
                end;
            end;
        end;
    end;
    if AP_FP_Greater(GROW*TSCAL,SMLNUM) then
    begin
        
        //
        // Use the Level 2 BLAS solve if the reciprocal of the bound on
        // elements of X is not too small.
        //
        if UPPER and NOTRAN or  not UPPER and  not NOTRAN then
        begin
            if NOunit then
            begin
                VD := A[N,N];
            end
            else
            begin
                VD := 1;
            end;
            X[N] := X[N]/VD;
            I:=N-1;
            while I>=1 do
            begin
                IP1 := I+1;
                if Upper then
                begin
                    V := APVDotProduct(@A[I][0], IP1, N, @X[0], IP1, N);
                end
                else
                begin
                    V := 0.0;
                    for i_ := IP1 to N do
                    begin
                        V := V + A[i_,I]*X[i_];
                    end;
                end;
                if NOunit then
                begin
                    VD := A[I,I];
                end
                else
                begin
                    VD := 1;
                end;
                X[I] := (X[I]-V)/VD;
                Dec(I);
            end;
        end
        else
        begin
            if NOunit then
            begin
                VD := A[1,1];
            end
            else
            begin
                VD := 1;
            end;
            X[1] := X[1]/VD;
            I:=2;
            while I<=N do
            begin
                IM1 := I-1;
                if Upper then
                begin
                    V := 0.0;
                    for i_ := 1 to IM1 do
                    begin
                        V := V + A[i_,I]*X[i_];
                    end;
                end
                else
                begin
                    V := APVDotProduct(@A[I][0], 1, IM1, @X[0], 1, IM1);
                end;
                if NOunit then
                begin
                    VD := A[I,I];
                end
                else
                begin
                    VD := 1;
                end;
                X[I] := (X[I]-V)/VD;
                Inc(I);
            end;
        end;
    end
    else
    begin
        
        //
        // Use a Level 1 BLAS solve, scaling intermediate results.
        //
        if AP_FP_Greater(XMAX,BIGNUM) then
        begin
            
            //
            // Scale X so that its components are less than or equal to
            // BIGNUM in absolute value.
            //
            S := BIGNUM/XMAX;
            APVMul(@X[0], 1, N, S);
            XMAX := BIGNUM;
        end;
        if NOTRAN then
        begin
            
            //
            // Solve A * x = b
            //
            J := JFIRST;
            while (JINC>0) and (J<=JLAST) or (JINC<0) and (J>=JLAST) do
            begin
                
                //
                // Compute x(j) = b(j) / A(j,j), scaling x if necessary.
                //
                XJ := ABSReal(X[J]);
                Flg := 0;
                if NOunit then
                begin
                    TJJS := A[J,J]*TSCAL;
                end
                else
                begin
                    TJJS := TSCAL;
                    if AP_FP_Eq(TSCAL,1) then
                    begin
                        Flg := 100;
                    end;
                end;
                if Flg<>100 then
                begin
                    TJJ := ABSReal(TJJS);
                    if AP_FP_Greater(TJJ,SMLNUM) then
                    begin
                        
                        //
                        // abs(A(j,j)) > SMLNUM:
                        //
                        if AP_FP_Less(TJJ,1) then
                        begin
                            if AP_FP_Greater(XJ,TJJ*BIGNUM) then
                            begin
                                
                                //
                                // Scale x by 1/b(j).
                                //
                                REC := 1/XJ;
                                APVMul(@X[0], 1, N, REC);
                                S := S*REC;
                                XMAX := XMAX*REC;
                            end;
                        end;
                        X[J] := X[J]/TJJS;
                        XJ := ABSReal(X[J]);
                    end
                    else
                    begin
                        if AP_FP_Greater(TJJ,0) then
                        begin
                            
                            //
                            // 0 < abs(A(j,j)) <= SMLNUM:
                            //
                            if AP_FP_Greater(XJ,TJJ*BIGNUM) then
                            begin
                                
                                //
                                // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                                // to avoid overflow when dividing by A(j,j).
                                //
                                REC := TJJ*BIGNUM/XJ;
                                if AP_FP_Greater(CNORM[J],1) then
                                begin
                                    
                                    //
                                    // Scale by 1/CNORM(j) to avoid overflow when
                                    // multiplying x(j) times column j.
                                    //
                                    REC := REC/CNORM[J];
                                end;
                                APVMul(@X[0], 1, N, REC);
                                S := S*REC;
                                XMAX := XMAX*REC;
                            end;
                            X[J] := X[J]/TJJS;
                            XJ := ABSReal(X[J]);
                        end
                        else
                        begin
                            
                            //
                            // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                            // scale = 0, and compute a solution to A*x = 0.
                            //
                            I:=1;
                            while I<=N do
                            begin
                                X[I] := 0;
                                Inc(I);
                            end;
                            X[J] := 1;
                            XJ := 1;
                            S := 0;
                            XMAX := 0;
                        end;
                    end;
                end;
                
                //
                // Scale x if necessary to avoid overflow when adding a
                // multiple of column j of A.
                //
                if AP_FP_Greater(XJ,1) then
                begin
                    REC := 1/XJ;
                    if AP_FP_Greater(CNORM[J],(BIGNUM-XMAX)*REC) then
                    begin
                        
                        //
                        // Scale x by 1/(2*abs(x(j))).
                        //
                        REC := REC*Double(0.5);
                        APVMul(@X[0], 1, N, REC);
                        S := S*REC;
                    end;
                end
                else
                begin
                    if AP_FP_Greater(XJ*CNORM[J],BIGNUM-XMAX) then
                    begin
                        
                        //
                        // Scale x by 1/2.
                        //
                        APVMul(@X[0], 1, N, 0.5);
                        S := S*Double(0.5);
                    end;
                end;
                if UPPER then
                begin
                    if J>1 then
                    begin
                        
                        //
                        // Compute the update
                        // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                        //
                        V := X[J]*TSCAL;
                        JM1 := J-1;
                        for i_ := 1 to JM1 do
                        begin
                            X[i_] := X[i_] - V*A[i_,J];
                        end;
                        I := 1;
                        K:=2;
                        while K<=J-1 do
                        begin
                            if AP_FP_Greater(AbsReal(X[K]),AbsReal(X[I])) then
                            begin
                                I := K;
                            end;
                            Inc(K);
                        end;
                        XMAX := ABSReal(X[I]);
                    end;
                end
                else
                begin
                    if J<N then
                    begin
                        
                        //
                        // Compute the update
                        // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                        //
                        JP1 := J+1;
                        V := X[J]*TSCAL;
                        for i_ := JP1 to N do
                        begin
                            X[i_] := X[i_] - V*A[i_,J];
                        end;
                        I := J+1;
                        K:=J+2;
                        while K<=N do
                        begin
                            if AP_FP_Greater(AbsReal(X[K]),AbsReal(X[I])) then
                            begin
                                I := K;
                            end;
                            Inc(K);
                        end;
                        XMAX := ABSReal(X[I]);
                    end;
                end;
                J := J+JINC;
            end;
        end
        else
        begin
            
            //
            // Solve A' * x = b
            //
            J := JFIRST;
            while (JINC>0) and (J<=JLAST) or (JINC<0) and (J>=JLAST) do
            begin
                
                //
                // Compute x(j) = b(j) - sum A(k,j)*x(k).
                //   k<>j
                //
                XJ := ABSReal(X[J]);
                USCAL := TSCAL;
                REC := 1/Max(XMAX, 1);
                if AP_FP_Greater(CNORM[J],(BIGNUM-XJ)*REC) then
                begin
                    
                    //
                    // If x(j) could overflow, scale x by 1/(2*XMAX).
                    //
                    REC := REC*Double(0.5);
                    if NOunit then
                    begin
                        TJJS := A[J,J]*TSCAL;
                    end
                    else
                    begin
                        TJJS := TSCAL;
                    end;
                    TJJ := ABSReal(TJJS);
                    if AP_FP_Greater(TJJ,1) then
                    begin
                        
                        //
                        // Divide by A(j,j) when scaling x if A(j,j) > 1.
                        //
                        REC := Min(1, REC*TJJ);
                        USCAL := USCAL/TJJS;
                    end;
                    if AP_FP_Less(REC,1) then
                    begin
                        APVMul(@X[0], 1, N, REC);
                        S := S*REC;
                        XMAX := XMAX*REC;
                    end;
                end;
                SUMJ := 0;
                if AP_FP_Eq(USCAL,1) then
                begin
                    
                    //
                    // If the scaling needed for A in the dot product is 1,
                    // call DDOT to perform the dot product.
                    //
                    if UPPER then
                    begin
                        if J>1 then
                        begin
                            JM1 := J-1;
                            SUMJ := 0.0;
                            for i_ := 1 to JM1 do
                            begin
                                SUMJ := SUMJ + A[i_,J]*X[i_];
                            end;
                        end
                        else
                        begin
                            SUMJ := 0;
                        end;
                    end
                    else
                    begin
                        if J<N then
                        begin
                            JP1 := J+1;
                            SUMJ := 0.0;
                            for i_ := JP1 to N do
                            begin
                                SUMJ := SUMJ + A[i_,J]*X[i_];
                            end;
                        end;
                    end;
                end
                else
                begin
                    
                    //
                    // Otherwise, use in-line code for the dot product.
                    //
                    if UPPER then
                    begin
                        I:=1;
                        while I<=J-1 do
                        begin
                            V := A[I,J]*USCAL;
                            SUMJ := SUMJ+V*X[I];
                            Inc(I);
                        end;
                    end
                    else
                    begin
                        if J<N then
                        begin
                            I:=J+1;
                            while I<=N do
                            begin
                                V := A[I,J]*USCAL;
                                SUMJ := SUMJ+V*X[I];
                                Inc(I);
                            end;
                        end;
                    end;
                end;
                if AP_FP_Eq(USCAL,TSCAL) then
                begin
                    
                    //
                    // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                    // was not used to scale the dotproduct.
                    //
                    X[J] := X[J]-SUMJ;
                    XJ := ABSReal(X[J]);
                    Flg := 0;
                    if NOunit then
                    begin
                        TJJS := A[J,J]*TSCAL;
                    end
                    else
                    begin
                        TJJS := TSCAL;
                        if AP_FP_Eq(TSCAL,1) then
                        begin
                            Flg := 150;
                        end;
                    end;
                    
                    //
                    // Compute x(j) = x(j) / A(j,j), scaling if necessary.
                    //
                    if Flg<>150 then
                    begin
                        TJJ := ABSReal(TJJS);
                        if AP_FP_Greater(TJJ,SMLNUM) then
                        begin
                            
                            //
                            // abs(A(j,j)) > SMLNUM:
                            //
                            if AP_FP_Less(TJJ,1) then
                            begin
                                if AP_FP_Greater(XJ,TJJ*BIGNUM) then
                                begin
                                    
                                    //
                                    // Scale X by 1/abs(x(j)).
                                    //
                                    REC := 1/XJ;
                                    APVMul(@X[0], 1, N, REC);
                                    S := S*REC;
                                    XMAX := XMAX*REC;
                                end;
                            end;
                            X[J] := X[J]/TJJS;
                        end
                        else
                        begin
                            if AP_FP_Greater(TJJ,0) then
                            begin
                                
                                //
                                // 0 < abs(A(j,j)) <= SMLNUM:
                                //
                                if AP_FP_Greater(XJ,TJJ*BIGNUM) then
                                begin
                                    
                                    //
                                    // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                                    //
                                    REC := TJJ*BIGNUM/XJ;
                                    APVMul(@X[0], 1, N, REC);
                                    S := S*REC;
                                    XMAX := XMAX*REC;
                                end;
                                X[J] := X[J]/TJJS;
                            end
                            else
                            begin
                                
                                //
                                // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                                // scale = 0, and compute a solution to A'*x = 0.
                                //
                                I:=1;
                                while I<=N do
                                begin
                                    X[I] := 0;
                                    Inc(I);
                                end;
                                X[J] := 1;
                                S := 0;
                                XMAX := 0;
                            end;
                        end;
                    end;
                end
                else
                begin
                    
                    //
                    // Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                    // product has already been divided by 1/A(j,j).
                    //
                    X[J] := X[J]/TJJS-SUMJ;
                end;
                XMAX := Max(XMAX, ABSReal(X[J]));
                J := J+JINC;
            end;
        end;
        S := S/TSCAL;
    end;
    
    //
    // Scale the column norms by 1/TSCAL for return.
    //
    if AP_FP_Neq(TSCAL,1) then
    begin
        V := 1/TSCAL;
        APVMul(@CNORM[0], 1, N, V);
    end;
end;


end.