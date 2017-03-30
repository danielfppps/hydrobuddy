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
unit reflections;
interface
uses Math, Sysutils, Ap;

procedure GenerateReflection(var X : TReal1DArray;
     N : AlglibInteger;
     var Tau : Double);
procedure ApplyReflectionFromTheLeft(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TReal1DArray);
procedure ApplyReflectionFromTheRight(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TReal1DArray);

implementation

(*************************************************************************
Generation of an elementary reflection transformation

The subroutine generates elementary reflection H of order N, so that, for
a given X, the following equality holds true:

    ( X(1) )   ( Beta )
H * (  ..  ) = (  0   )
    ( X(n) )   (  0   )

where
              ( V(1) )
H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array whose index ranges within [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced with vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau. If X is a null vector, Tau equals 0,
            otherwise 1 <= Tau <= 2.

This subroutine is the modification of the DLARFG subroutines from
the LAPACK library.

MODIFICATIONS:
    24.12.2005 sign(Alpha) was replaced with an analogous to the Fortran SIGN code.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure GenerateReflection(var X : TReal1DArray;
     N : AlglibInteger;
     var Tau : Double);
var
    J : AlglibInteger;
    Alpha : Double;
    XNORM : Double;
    V : Double;
    Beta : Double;
    MX : Double;
    S : Double;
begin
    if N<=1 then
    begin
        Tau := 0;
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
        MX := Max(AbsReal(X[J]), MX);
        Inc(J);
    end;
    S := 1;
    if AP_FP_Neq(MX,0) then
    begin
        if AP_FP_Less_Eq(MX,MinRealNumber/MachineEpsilon) then
        begin
            S := MinRealNumber/MachineEpsilon;
            V := 1/S;
            APVMul(@X[0], 1, N, V);
            MX := MX*V;
        end
        else
        begin
            if AP_FP_Greater_Eq(MX,MaxRealNumber*MachineEpsilon) then
            begin
                S := MaxRealNumber*MachineEpsilon;
                V := 1/S;
                APVMul(@X[0], 1, N, V);
                MX := MX*V;
            end;
        end;
    end;
    
    //
    // XNORM = DNRM2( N-1, X, INCX )
    //
    Alpha := X[1];
    XNORM := 0;
    if AP_FP_Neq(MX,0) then
    begin
        J:=2;
        while J<=N do
        begin
            XNORM := XNORM+AP_Sqr(X[J]/MX);
            Inc(J);
        end;
        XNORM := Sqrt(XNORM)*MX;
    end;
    if AP_FP_Eq(XNORM,0) then
    begin
        
        //
        // H  =  I
        //
        TAU := 0;
        X[1] := X[1]*S;
        Exit;
    end;
    
    //
    // general case
    //
    MX := Max(AbsReal(Alpha), AbsReal(XNORM));
    Beta := -MX*Sqrt(AP_Sqr(Alpha/MX)+AP_Sqr(XNORM/MX));
    if AP_FP_Less(Alpha,0) then
    begin
        Beta := -Beta;
    end;
    TAU := (BETA-ALPHA)/BETA;
    V := 1/(Alpha-Beta);
    APVMul(@X[0], 2, N, V);
    X[1] := Beta;
    
    //
    // Scale back outputs
    //
    X[1] := X[1]*S;
end;


(*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm pre-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the elements
of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..M2-M1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ApplyReflectionFromTheLeft(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TReal1DArray);
var
    T : Double;
    I : AlglibInteger;
    VM : AlglibInteger;
begin
    if AP_FP_Eq(Tau,0) or (N1>N2) or (M1>M2) then
    begin
        Exit;
    end;
    
    //
    // w := C' * v
    //
    VM := M2-M1+1;
    I:=N1;
    while I<=N2 do
    begin
        WORK[I] := 0;
        Inc(I);
    end;
    I:=M1;
    while I<=M2 do
    begin
        T := V[I+1-M1];
        APVAdd(@WORK[0], N1, N2, @C[I][0], N1, N2, T);
        Inc(I);
    end;
    
    //
    // C := C - tau * v * w'
    //
    I:=M1;
    while I<=M2 do
    begin
        T := V[I-M1+1]*TAU;
        APVSub(@C[I][0], N1, N2, @WORK[0], N1, N2, T);
        Inc(I);
    end;
end;


(*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm post-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the
elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..N2-N1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ApplyReflectionFromTheRight(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     var WORK : TReal1DArray);
var
    T : Double;
    I : AlglibInteger;
    VM : AlglibInteger;
begin
    if AP_FP_Eq(Tau,0) or (N1>N2) or (M1>M2) then
    begin
        Exit;
    end;
    VM := N2-N1+1;
    I:=M1;
    while I<=M2 do
    begin
        T := APVDotProduct(@C[I][0], N1, N2, @V[0], 1, VM);
        T := T*TAU;
        APVSub(@C[I][0], N1, N2, @V[0], 1, VM, T);
        Inc(I);
    end;
end;


end.