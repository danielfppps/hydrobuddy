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
unit creflections;
interface
uses Math, Sysutils, Ap;

procedure ComplexGenerateReflection(var X : TComplex1DArray;
     N : AlglibInteger;
     var Tau : Complex);
procedure ComplexApplyReflectionFromTheLeft(var C : TComplex2DArray;
     Tau : Complex;
     const V : TComplex1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TComplex1DArray);
procedure ComplexApplyReflectionFromTheRight(var C : TComplex2DArray;
     Tau : Complex;
     var V : TComplex1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TComplex1DArray);

implementation

(*************************************************************************
Generation of an elementary complex reflection transformation

The subroutine generates elementary complex reflection H of  order  N,  so
that, for a given X, the following equality holds true:

     ( X(1) )   ( Beta )
H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number
     ( X(n) )   (  0   )

where

              ( V(1) )
H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array with elements [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced by vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau.

This subroutine is the modification of CLARFG subroutines  from the LAPACK
library. It has similar functionality except for the fact that it  doesn’t
handle errors when intermediate results cause an overflow.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ComplexGenerateReflection(var X : TComplex1DArray;
     N : AlglibInteger;
     var Tau : Complex);
var
    J : AlglibInteger;
    ALPHA : Complex;
    ALPHI : Double;
    ALPHR : Double;
    BETA : Double;
    XNORM : Double;
    MX : Double;
    T : Complex;
    S : Double;
    V : Complex;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        TAU := C_Complex(0);
        Exit;
    end;
    
    //
    // Scale if needed (to avoid overflow/underflow during intermediate
    // calculations).
    //
    MX := 0;
    J:=1;
    while J<=N do
    begin
        MX := Max(AbsComplex(X[J]), MX);
        Inc(J);
    end;
    S := 1;
    if AP_FP_Neq(MX,0) then
    begin
        if AP_FP_Less(MX,1) then
        begin
            S := Sqrt(MinRealNumber);
            V := C_Complex(1/S);
            for i_ := 1 to N do
            begin
                X[i_] := C_Mul(V, X[i_]);
            end;
        end
        else
        begin
            S := Sqrt(MaxRealNumber);
            V := C_Complex(1/S);
            for i_ := 1 to N do
            begin
                X[i_] := C_Mul(V, X[i_]);
            end;
        end;
    end;
    
    //
    // calculate
    //
    ALPHA := X[1];
    MX := 0;
    J:=2;
    while J<=N do
    begin
        MX := Max(AbsComplex(X[J]), MX);
        Inc(J);
    end;
    XNORM := 0;
    if AP_FP_Neq(MX,0) then
    begin
        J:=2;
        while J<=N do
        begin
            T := C_DivR(X[J],MX);
            XNORM := XNORM+C_Mul(T,Conj(T)).X;
            Inc(J);
        end;
        XNORM := Sqrt(XNORM)*MX;
    end;
    ALPHR := ALPHA.X;
    ALPHI := ALPHA.Y;
    if AP_FP_Eq(XNORM,0) and AP_FP_Eq(ALPHI,0) then
    begin
        TAU := C_Complex(0);
        X[1] := C_MulR(X[1],S);
        Exit;
    end;
    MX := Max(AbsReal(ALPHR), AbsReal(ALPHI));
    MX := Max(MX, AbsReal(XNORM));
    BETA := -MX*Sqrt(AP_Sqr(ALPHR/MX)+AP_Sqr(ALPHI/MX)+AP_Sqr(XNORM/MX));
    if AP_FP_Less(ALPHR,0) then
    begin
        BETA := -BETA;
    end;
    TAU.X := (BETA-ALPHR)/BETA;
    TAU.Y := -ALPHI/BETA;
    ALPHA := C_RDiv(1,C_SubR(ALPHA,BETA));
    if N>1 then
    begin
        for i_ := 2 to N do
        begin
            X[i_] := C_Mul(ALPHA, X[i_]);
        end;
    end;
    ALPHA := C_Complex(BETA);
    X[1] := ALPHA;
    
    //
    // Scale back
    //
    X[1] := C_MulR(X[1],S);
end;


(*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  pre-multiplies  the  matrix  by  an  elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description of the GenerateReflection). Not the whole matrix  but  only  a
part of it is transformed (rows from M1 to M2, columns from N1 to N2). Only
the elements of this submatrix are changed.

Note: the matrix is multiplied by H, not by H'.   If  it  is  required  to
multiply the matrix by H', it is necessary to pass Conj(Tau) instead of Tau.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..M2-M1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ComplexApplyReflectionFromTheLeft(var C : TComplex2DArray;
     Tau : Complex;
     const V : TComplex1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TComplex1DArray);
var
    T : Complex;
    I : AlglibInteger;
    VM : AlglibInteger;
    i_ : AlglibInteger;
begin
    if C_EqualR(Tau,0) or (N1>N2) or (M1>M2) then
    begin
        Exit;
    end;
    
    //
    // w := C^T * conj(v)
    //
    VM := M2-M1+1;
    I:=N1;
    while I<=N2 do
    begin
        WORK[I] := C_Complex(0);
        Inc(I);
    end;
    I:=M1;
    while I<=M2 do
    begin
        T := Conj(V[I+1-M1]);
        for i_ := N1 to N2 do
        begin
            WORK[i_] := C_Add(WORK[i_], C_Mul(T, C[I,i_]));
        end;
        Inc(I);
    end;
    
    //
    // C := C - tau * v * w^T
    //
    I:=M1;
    while I<=M2 do
    begin
        T := C_Mul(V[I-M1+1],TAU);
        for i_ := N1 to N2 do
        begin
            C[I,i_] := C_Sub(C[I,i_], C_Mul(T, WORK[i_]));
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  post-multiplies  the  matrix  by  an elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description  of  the  GenerateReflection). Not the whole matrix but only a
part  of  it  is  transformed (rows from M1 to M2, columns from N1 to N2).
Only the elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..N2-N1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ComplexApplyReflectionFromTheRight(var C : TComplex2DArray;
     Tau : Complex;
     var V : TComplex1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TComplex1DArray);
var
    T : Complex;
    I : AlglibInteger;
    VM : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if C_EqualR(Tau,0) or (N1>N2) or (M1>M2) then
    begin
        Exit;
    end;
    
    //
    // w := C * v
    //
    VM := N2-N1+1;
    I:=M1;
    while I<=M2 do
    begin
        i1_ := (1)-(N1);
        T := C_Complex(0.0);
        for i_ := N1 to N2 do
        begin
            T := C_Add(T,C_Mul(C[I,i_],V[i_+i1_]));
        end;
        WORK[I] := T;
        Inc(I);
    end;
    
    //
    // C := C - w * conj(v^T)
    //
    for i_ := 1 to VM do
    begin
        V[i_] := Conj(V[i_]);
    end;
    I:=M1;
    while I<=M2 do
    begin
        T := C_Mul(WORK[I],TAU);
        i1_ := (1) - (N1);
        for i_ := N1 to N2 do
        begin
            C[I,i_] := C_Sub(C[I,i_], C_Mul(T, V[i_+i1_]));
        end;
        Inc(I);
    end;
    for i_ := 1 to VM do
    begin
        V[i_] := Conj(V[i_]);
    end;
end;


end.