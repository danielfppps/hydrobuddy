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
unit estnorm;
interface
uses Math, Sysutils, Ap;

procedure IterativeEstimate1Norm(N : AlglibInteger;
     var V : TReal1DArray;
     var X : TReal1DArray;
     var ISGN : TInteger1DArray;
     var EST : Double;
     var KASE : AlglibInteger);
function DemoIterativeEstimate1Norm(const A : TReal2DArray;
     N : AlglibInteger):Double;

implementation

(*************************************************************************
Matrix norm estimation

The algorithm estimates the 1-norm of square matrix A  on  the  assumption
that the multiplication of matrix  A  by  the  vector  is  available  (the
iterative method is used). It is recommended to use this algorithm  if  it
is hard  to  calculate  matrix  elements  explicitly  (for  example,  when
estimating the inverse matrix norm).

The algorithm uses back communication for multiplying the  vector  by  the
matrix.  If  KASE=0  after  returning from a subroutine, its execution was
completed successfully, otherwise it is required to multiply the  returned
vector by matrix A and call the subroutine again.

The DemoIterativeEstimateNorm subroutine shows a simple example.

Parameters:
    N       -   size of matrix A.
    V       -   vector.   It is initialized by the subroutine on the first
                call. It is then passed into it on repeated calls.
    X       -   if KASE<>0, it contains the vector to be replaced by:
                    A * X,      if KASE=1
                    A^T * X,    if KASE=2
                Array whose index ranges within [1..N].
    ISGN    -   vector. It is initialized by the subroutine on  the  first
                call. It is then passed into it on repeated calls.
    EST     -   if KASE=0, it contains the lower boundary of the matrix
                norm estimate.
    KASE    -   on the first call, it should be equal to 0. After the last
                return, it is equal to 0 (EST contains the  matrix  norm),
                on intermediate returns it can be equal to 1 or 2 depending
                on the operation to be performed on vector X.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
procedure IterativeEstimate1Norm(N : AlglibInteger;
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
        SetLength(V, N+3+1);
        SetLength(X, N+1);
        SetLength(ISGN, N+4+1);
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


(*************************************************************************
Example of usage of an IterativeEstimateNorm subroutine

Input parameters:
    A   -   matrix.
            Array whose indexes range within [1..N, 1..N].

Return:
    Matrix norm estimated by the subroutine.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
function DemoIterativeEstimate1Norm(const A : TReal2DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    S : Double;
    X : TReal1DArray;
    T : TReal1DArray;
    V : TReal1DArray;
    IV : TInteger1DArray;
    KASE : AlglibInteger;
    i_ : AlglibInteger;
begin
    KASE := 0;
    SetLength(T, N+1);
    IterativeEstimate1Norm(N, V, X, IV, Result, KASE);
    while KASE<>0 do
    begin
        if KASE=1 then
        begin
            I:=1;
            while I<=N do
            begin
                S := APVDotProduct(@A[I][0], 1, N, @X[0], 1, N);
                T[I] := S;
                Inc(I);
            end;
        end
        else
        begin
            I:=1;
            while I<=N do
            begin
                S := 0.0;
                for i_ := 1 to N do
                begin
                    S := S + A[i_,I]*X[i_];
                end;
                T[I] := S;
                Inc(I);
            end;
        end;
        APVMove(@X[0], 1, N, @T[0], 1, N);
        IterativeEstimate1Norm(N, V, X, IV, Result, KASE);
    end;
end;


end.