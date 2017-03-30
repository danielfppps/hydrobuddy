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
unit srcond;
interface
uses Math, Sysutils, Ap, ldlt, ssolve, estnorm;

function SMatrixRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function SMatrixLDLTRCond(const L : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function RCondSymmetric(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function RCondLDLT(const L : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
procedure InternalLDLTRCond(const L : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RCOND : Double);

implementation

(*************************************************************************
Condition number estimate of a symmetric matrix

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   symmetric definite matrix which is given by its upper or
                lower triangle depending on IsUpper.
                Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A))
*************************************************************************)
function SMatrixRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    A1 : TReal2DArray;
begin
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        if IsUpper then
        begin
            J:=I;
            while J<=N do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
        end
        else
        begin
            J:=1;
            while J<=I do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
        end;
        Inc(I);
    end;
    Result := RCondSymmetric(A1, N, IsUpper);
end;


(*************************************************************************
Condition number estimate of a matrix given by LDLT-decomposition

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    L       -   LDLT-decomposition of matrix A given by the upper or lower
                triangle depending on IsUpper.
                Output of SMatrixLDLT subroutine.
    Pivots  -   table of permutations which were made during LDLT-decomposition,
                Output of SMatrixLDLT subroutine.
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A))
*************************************************************************)
function SMatrixLDLTRCond(const L : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : TReal2DArray;
    P1 : TInteger1DArray;
begin
    SetLength(L1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        if IsUpper then
        begin
            J:=I;
            while J<=N do
            begin
                L1[I,J] := L[I-1,J-1];
                Inc(J);
            end;
        end
        else
        begin
            J:=1;
            while J<=I do
            begin
                L1[I,J] := L[I-1,J-1];
                Inc(J);
            end;
        end;
        Inc(I);
    end;
    SetLength(P1, N+1);
    I:=1;
    while I<=N do
    begin
        if Pivots[I-1]>=0 then
        begin
            P1[I] := Pivots[I-1]+1;
        end
        else
        begin
            P1[I] := -(Pivots[I-1]+N+1);
        end;
        Inc(I);
    end;
    Result := RCondLDLT(L1, P1, N, IsUpper);
end;


function RCondSymmetric(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    IM : AlglibInteger;
    JM : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Nrm := 0;
    J:=1;
    while J<=N do
    begin
        V := 0;
        I:=1;
        while I<=N do
        begin
            IM := I;
            JM := J;
            if IsUpper and (J<I) then
            begin
                IM := J;
                JM := I;
            end;
            if  not IsUpper and (J>I) then
            begin
                IM := J;
                JM := I;
            end;
            V := V+AbsReal(A[IM,JM]);
            Inc(I);
        end;
        Nrm := Max(Nrm, V);
        Inc(J);
    end;
    LDLTDecomposition(A, N, IsUpper, Pivots);
    InternalLDLTRCond(A, Pivots, N, IsUpper, True, Nrm, V);
    Result := V;
end;


function RCondLDLT(const L : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    V : Double;
begin
    InternalLDLTRCond(L, Pivots, N, IsUpper, False, 0, V);
    Result := V;
end;


procedure InternalLDLTRCond(const L : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RCOND : Double);
var
    I : AlglibInteger;
    KASE : AlglibInteger;
    K : AlglibInteger;
    KM1 : AlglibInteger;
    KM2 : AlglibInteger;
    KP1 : AlglibInteger;
    KP2 : AlglibInteger;
    AINVNM : Double;
    WORK0 : TReal1DArray;
    WORK1 : TReal1DArray;
    WORK2 : TReal1DArray;
    IWORK : TInteger1DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    Assert(N>=0);
    
    //
    // Check that the diagonal matrix D is nonsingular.
    //
    RCond := 0;
    if IsUpper then
    begin
        I:=N;
        while I>=1 do
        begin
            if (Pivots[I]>0) and AP_FP_Eq(L[I,I],0) then
            begin
                Exit;
            end;
            Dec(I);
        end;
    end
    else
    begin
        I:=1;
        while I<=N do
        begin
            if (Pivots[I]>0) and AP_FP_Eq(L[I,I],0) then
            begin
                Exit;
            end;
            Inc(I);
        end;
    end;
    
    //
    // Estimate the norm of A.
    //
    if  not IsNormProvided then
    begin
        KASE := 0;
        ANORM := 0;
        while True do
        begin
            IterativeEstimate1Norm(N, WORK1, WORK0, IWORK, ANORM, KASE);
            if KASE=0 then
            begin
                Break;
            end;
            if IsUpper then
            begin
                
                //
                // Multiply by U'
                //
                K := N;
                while K>=1 do
                begin
                    if Pivots[K]>0 then
                    begin
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K];
                        WORK0[K] := WORK0[Pivots[K]];
                        WORK0[Pivots[K]] := V;
                        
                        //
                        // U(k)
                        //
                        KM1 := K-1;
                        V := 0.0;
                        for i_ := 1 to KM1 do
                        begin
                            V := V + WORK0[i_]*L[i_,K];
                        end;
                        WORK0[K] := WORK0[K]+V;
                        
                        //
                        // Next k
                        //
                        K := K-1;
                    end
                    else
                    begin
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K-1];
                        WORK0[K-1] := WORK0[-Pivots[K-1]];
                        WORK0[-Pivots[K-1]] := V;
                        
                        //
                        // U(k)
                        //
                        KM1 := K-1;
                        KM2 := K-2;
                        V := 0.0;
                        for i_ := 1 to KM2 do
                        begin
                            V := V + WORK0[i_]*L[i_,KM1];
                        end;
                        WORK0[KM1] := WORK0[KM1]+V;
                        V := 0.0;
                        for i_ := 1 to KM2 do
                        begin
                            V := V + WORK0[i_]*L[i_,K];
                        end;
                        WORK0[K] := WORK0[K]+V;
                        
                        //
                        // Next k
                        //
                        K := K-2;
                    end;
                end;
                
                //
                // Multiply by D
                //
                K := N;
                while K>=1 do
                begin
                    if Pivots[K]>0 then
                    begin
                        WORK0[K] := WORK0[K]*L[K,K];
                        K := K-1;
                    end
                    else
                    begin
                        V := WORK0[K-1];
                        WORK0[K-1] := L[K-1,K-1]*WORK0[K-1]+L[K-1,K]*WORK0[K];
                        WORK0[K] := L[K-1,K]*V+L[K,K]*WORK0[K];
                        K := K-2;
                    end;
                end;
                
                //
                // Multiply by U
                //
                K := 1;
                while K<=N do
                begin
                    if Pivots[K]>0 then
                    begin
                        
                        //
                        // U(k)
                        //
                        KM1 := K-1;
                        V := WORK0[K];
                        for i_ := 1 to KM1 do
                        begin
                            WORK0[i_] := WORK0[i_] + V*L[i_,K];
                        end;
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K];
                        WORK0[K] := WORK0[Pivots[K]];
                        WORK0[Pivots[K]] := V;
                        
                        //
                        // Next k
                        //
                        K := K+1;
                    end
                    else
                    begin
                        
                        //
                        // U(k)
                        //
                        KM1 := K-1;
                        KP1 := K+1;
                        V := WORK0[K];
                        for i_ := 1 to KM1 do
                        begin
                            WORK0[i_] := WORK0[i_] + V*L[i_,K];
                        end;
                        V := WORK0[KP1];
                        for i_ := 1 to KM1 do
                        begin
                            WORK0[i_] := WORK0[i_] + V*L[i_,KP1];
                        end;
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K];
                        WORK0[K] := WORK0[-Pivots[K]];
                        WORK0[-Pivots[K]] := V;
                        
                        //
                        // Next k
                        //
                        K := K+2;
                    end;
                end;
            end
            else
            begin
                
                //
                // Multiply by L'
                //
                K := 1;
                while K<=N do
                begin
                    if Pivots[K]>0 then
                    begin
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K];
                        WORK0[K] := WORK0[Pivots[K]];
                        WORK0[Pivots[K]] := V;
                        
                        //
                        // L(k)
                        //
                        KP1 := K+1;
                        V := 0.0;
                        for i_ := KP1 to N do
                        begin
                            V := V + WORK0[i_]*L[i_,K];
                        end;
                        WORK0[K] := WORK0[K]+V;
                        
                        //
                        // Next k
                        //
                        K := K+1;
                    end
                    else
                    begin
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K+1];
                        WORK0[K+1] := WORK0[-Pivots[K+1]];
                        WORK0[-Pivots[K+1]] := V;
                        
                        //
                        // L(k)
                        //
                        KP1 := K+1;
                        KP2 := K+2;
                        V := 0.0;
                        for i_ := KP2 to N do
                        begin
                            V := V + WORK0[i_]*L[i_,K];
                        end;
                        WORK0[K] := WORK0[K]+V;
                        V := 0.0;
                        for i_ := KP2 to N do
                        begin
                            V := V + WORK0[i_]*L[i_,KP1];
                        end;
                        WORK0[KP1] := WORK0[KP1]+V;
                        
                        //
                        // Next k
                        //
                        K := K+2;
                    end;
                end;
                
                //
                // Multiply by D
                //
                K := N;
                while K>=1 do
                begin
                    if Pivots[K]>0 then
                    begin
                        WORK0[K] := WORK0[K]*L[K,K];
                        K := K-1;
                    end
                    else
                    begin
                        V := WORK0[K-1];
                        WORK0[K-1] := L[K-1,K-1]*WORK0[K-1]+L[K,K-1]*WORK0[K];
                        WORK0[K] := L[K,K-1]*V+L[K,K]*WORK0[K];
                        K := K-2;
                    end;
                end;
                
                //
                // Multiply by L
                //
                K := N;
                while K>=1 do
                begin
                    if Pivots[K]>0 then
                    begin
                        
                        //
                        // L(k)
                        //
                        KP1 := K+1;
                        V := WORK0[K];
                        for i_ := KP1 to N do
                        begin
                            WORK0[i_] := WORK0[i_] + V*L[i_,K];
                        end;
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K];
                        WORK0[K] := WORK0[Pivots[K]];
                        WORK0[Pivots[K]] := V;
                        
                        //
                        // Next k
                        //
                        K := K-1;
                    end
                    else
                    begin
                        
                        //
                        // L(k)
                        //
                        KP1 := K+1;
                        KM1 := K-1;
                        V := WORK0[K];
                        for i_ := KP1 to N do
                        begin
                            WORK0[i_] := WORK0[i_] + V*L[i_,K];
                        end;
                        V := WORK0[KM1];
                        for i_ := KP1 to N do
                        begin
                            WORK0[i_] := WORK0[i_] + V*L[i_,KM1];
                        end;
                        
                        //
                        // P(k)
                        //
                        V := WORK0[K];
                        WORK0[K] := WORK0[-Pivots[K]];
                        WORK0[-Pivots[K]] := V;
                        
                        //
                        // Next k
                        //
                        K := K-2;
                    end;
                end;
            end;
        end;
    end;
    
    //
    // Quick return if possible
    //
    RCOND := 0;
    if N=0 then
    begin
        RCOND := 1;
        Exit;
    end;
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    
    //
    // Estimate the 1-norm of inv(A).
    //
    KASE := 0;
    while True do
    begin
        IterativeEstimate1Norm(N, WORK1, WORK0, IWORK, AINVNM, KASE);
        if KASE=0 then
        begin
            Break;
        end;
        SolveSystemLDLT(L, Pivots, WORK0, N, IsUpper, WORK2);
        APVMove(@WORK0[0], 1, N, @WORK2[0], 1, N);
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        V := 1/AINVNM;
        RCOND := V/ANORM;
    end;
end;


end.