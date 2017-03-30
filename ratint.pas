{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007-2009, Sergey Bochkanov (ALGLIB project).

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
unit ratint;
interface
uses Math, Sysutils, Ap, tsort, ratinterpolation, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit;

type
(*************************************************************************
Barycentric interpolant.
*************************************************************************)
BarycentricInterpolant = record
    N : AlglibInteger;
    SY : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
end;


(*************************************************************************
Barycentric fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************)
BarycentricFitReport = record
    TaskRCond : Double;
    DBest : AlglibInteger;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    MaxError : Double;
end;



function BarycentricCalc(const B : BarycentricInterpolant; T : Double):Double;
procedure BarycentricDiff1(const B : BarycentricInterpolant;
     T : Double;
     var F : Double;
     var DF : Double);
procedure BarycentricDiff2(const B : BarycentricInterpolant;
     T : Double;
     var F : Double;
     var DF : Double;
     var D2F : Double);
procedure BarycentricLinTransX(var B : BarycentricInterpolant;
     CA : Double;
     CB : Double);
procedure BarycentricLinTransY(var B : BarycentricInterpolant;
     CA : Double;
     CB : Double);
procedure BarycentricUnpack(const B : BarycentricInterpolant;
     var N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray;
     var W : TReal1DArray);
procedure BarycentricSerialize(const B : BarycentricInterpolant;
     var RA : TReal1DArray;
     var RALen : AlglibInteger);
procedure BarycentricUnserialize(const RA : TReal1DArray;
     var B : BarycentricInterpolant);
procedure BarycentricCopy(const B : BarycentricInterpolant;
     var B2 : BarycentricInterpolant);
procedure BarycentricBuildXYW(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     var B : BarycentricInterpolant);
procedure BarycentricBuildFloaterHormann(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     D : AlglibInteger;
     var B : BarycentricInterpolant);
procedure BarycentricFitFloaterHormannWC(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var B : BarycentricInterpolant;
     var Rep : BarycentricFitReport);
procedure BarycentricFitFloaterHormann(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var B : BarycentricInterpolant;
     var Rep : BarycentricFitReport);

implementation

const
    BRCVNum = 10;

procedure BarycentricNormalize(var B : BarycentricInterpolant);forward;
procedure BarycentricCalcBasis(const B : BarycentricInterpolant;
     T : Double;
     var Y : TReal1DArray);forward;
procedure BarycentricFitWCFixedD(X : TReal1DArray;
     Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     XC : TReal1DArray;
     YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     D : AlglibInteger;
     var Info : AlglibInteger;
     var B : BarycentricInterpolant;
     var Rep : BarycentricFitReport);forward;


(*************************************************************************
Rational interpolation using barycentric formula

F(t) = SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))

Input parameters:
    B   -   barycentric interpolant built with one of model building
            subroutines.
    T   -   interpolation point

Result:
    barycentric interpolant F(t)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
function BarycentricCalc(const B : BarycentricInterpolant; T : Double):Double;
var
    S1 : Double;
    S2 : Double;
    S : Double;
    V : Double;
    I : AlglibInteger;
begin
    
    //
    // special case: N=1
    //
    if B.N=1 then
    begin
        Result := B.SY*B.Y[0];
        Exit;
    end;
    
    //
    // Here we assume that task is normalized, i.e.:
    // 1. abs(Y[i])<=1
    // 2. abs(W[i])<=1
    // 3. X[] is ordered
    //
    S := AbsReal(T-B.X[0]);
    I:=0;
    while I<=B.N-1 do
    begin
        V := B.X[I];
        if AP_FP_Eq(V,T) then
        begin
            Result := B.SY*B.Y[I];
            Exit;
        end;
        V := AbsReal(T-V);
        if AP_FP_Less(V,S) then
        begin
            S := V;
        end;
        Inc(I);
    end;
    S1 := 0;
    S2 := 0;
    I:=0;
    while I<=B.N-1 do
    begin
        V := S/(T-B.X[I]);
        V := V*B.W[I];
        S1 := S1+V*B.Y[I];
        S2 := S2+V;
        Inc(I);
    end;
    Result := B.SY*S1/S2;
end;


(*************************************************************************
Differentiation of barycentric interpolant: first derivative.

Algorithm used in this subroutine is very robust and should not fail until
provided with values too close to MaxRealNumber  (usually  MaxRealNumber/N
or greater will overflow).

INPUT PARAMETERS:
    B   -   barycentric interpolant built with one of model building
            subroutines.
    T   -   interpolation point

OUTPUT PARAMETERS:
    F   -   barycentric interpolant at T
    DF  -   first derivative
    
NOTE


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricDiff1(const B : BarycentricInterpolant;
     T : Double;
     var F : Double;
     var DF : Double);
var
    V : Double;
    VV : Double;
    I : AlglibInteger;
    K : AlglibInteger;
    N0 : Double;
    N1 : Double;
    D0 : Double;
    D1 : Double;
    S0 : Double;
    S1 : Double;
    XK : Double;
    XI : Double;
    XMin : Double;
    XMax : Double;
    XScale1 : Double;
    XOffs1 : Double;
    XScale2 : Double;
    XOffs2 : Double;
    XPrev : Double;
begin
    
    //
    // special case: N=1
    //
    if B.N=1 then
    begin
        F := B.SY*B.Y[0];
        DF := 0;
        Exit;
    end;
    if AP_FP_Eq(B.SY,0) then
    begin
        F := 0;
        DF := 0;
        Exit;
    end;
    Assert(AP_FP_Greater(B.SY,0), 'BarycentricDiff1: internal error');
    
    //
    // We assume than N>1 and B.SY>0. Find:
    // 1. pivot point (X[i] closest to T)
    // 2. width of interval containing X[i]
    //
    V := AbsReal(B.X[0]-T);
    K := 0;
    XMin := B.X[0];
    XMax := B.X[0];
    I:=1;
    while I<=B.N-1 do
    begin
        VV := B.X[I];
        if AP_FP_Less(AbsReal(VV-T),V) then
        begin
            V := AbsReal(VV-T);
            K := I;
        end;
        XMin := Min(XMin, VV);
        XMax := Max(XMax, VV);
        Inc(I);
    end;
    
    //
    // pivot point found, calculate dNumerator and dDenominator
    //
    XScale1 := 1/(XMax-XMin);
    XOffs1 := -XMin/(XMax-XMin)+1;
    XScale2 := 2;
    XOffs2 := -3;
    T := T*XScale1+XOffs1;
    T := T*XScale2+XOffs2;
    XK := B.X[K];
    XK := XK*XScale1+XOffs1;
    XK := XK*XScale2+XOffs2;
    V := T-XK;
    N0 := 0;
    N1 := 0;
    D0 := 0;
    D1 := 0;
    XPrev := -2;
    I:=0;
    while I<=B.N-1 do
    begin
        XI := B.X[I];
        XI := XI*XScale1+XOffs1;
        XI := XI*XScale2+XOffs2;
        Assert(AP_FP_Greater(XI,XPrev), 'BarycentricDiff1: points are too close!');
        XPrev := XI;
        if I<>K then
        begin
            VV := AP_Sqr(T-XI);
            S0 := (T-XK)/(T-XI);
            S1 := (XK-XI)/VV;
        end
        else
        begin
            S0 := 1;
            S1 := 0;
        end;
        VV := B.W[I]*B.Y[I];
        N0 := N0+S0*VV;
        N1 := N1+S1*VV;
        VV := B.W[I];
        D0 := D0+S0*VV;
        D1 := D1+S1*VV;
        Inc(I);
    end;
    F := B.SY*N0/D0;
    DF := (N1*D0-N0*D1)/AP_Sqr(D0);
    if AP_FP_Neq(DF,0) then
    begin
        DF := Sign(DF)*Exp(Ln(AbsReal(DF))+Ln(B.SY)+Ln(XScale1)+Ln(XScale2));
    end;
end;


(*************************************************************************
Differentiation of barycentric interpolant: first/second derivatives.

INPUT PARAMETERS:
    B   -   barycentric interpolant built with one of model building
            subroutines.
    T   -   interpolation point

OUTPUT PARAMETERS:
    F   -   barycentric interpolant at T
    DF  -   first derivative
    D2F -   second derivative

NOTE: this algorithm may fail due to overflow/underflor if  used  on  data
whose values are close to MaxRealNumber or MinRealNumber.  Use more robust
BarycentricDiff1() subroutine in such cases.


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricDiff2(const B : BarycentricInterpolant;
     T : Double;
     var F : Double;
     var DF : Double;
     var D2F : Double);
var
    V : Double;
    VV : Double;
    I : AlglibInteger;
    K : AlglibInteger;
    N0 : Double;
    N1 : Double;
    N2 : Double;
    D0 : Double;
    D1 : Double;
    D2 : Double;
    S0 : Double;
    S1 : Double;
    S2 : Double;
    XK : Double;
    XI : Double;
begin
    F := 0;
    DF := 0;
    D2F := 0;
    
    //
    // special case: N=1
    //
    if B.N=1 then
    begin
        F := B.SY*B.Y[0];
        DF := 0;
        D2F := 0;
        Exit;
    end;
    if AP_FP_Eq(B.SY,0) then
    begin
        F := 0;
        DF := 0;
        D2F := 0;
        Exit;
    end;
    Assert(AP_FP_Greater(B.SY,0), 'BarycentricDiff: internal error');
    
    //
    // We assume than N>1 and B.SY>0. Find:
    // 1. pivot point (X[i] closest to T)
    // 2. width of interval containing X[i]
    //
    V := AbsReal(B.X[0]-T);
    K := 0;
    I:=1;
    while I<=B.N-1 do
    begin
        VV := B.X[I];
        if AP_FP_Less(AbsReal(VV-T),V) then
        begin
            V := AbsReal(VV-T);
            K := I;
        end;
        Inc(I);
    end;
    
    //
    // pivot point found, calculate dNumerator and dDenominator
    //
    XK := B.X[K];
    V := T-XK;
    N0 := 0;
    N1 := 0;
    N2 := 0;
    D0 := 0;
    D1 := 0;
    D2 := 0;
    I:=0;
    while I<=B.N-1 do
    begin
        if I<>K then
        begin
            XI := B.X[I];
            VV := AP_Sqr(T-XI);
            S0 := (T-XK)/(T-XI);
            S1 := (XK-XI)/VV;
            S2 := -2*(XK-XI)/(VV*(T-XI));
        end
        else
        begin
            S0 := 1;
            S1 := 0;
            S2 := 0;
        end;
        VV := B.W[I]*B.Y[I];
        N0 := N0+S0*VV;
        N1 := N1+S1*VV;
        N2 := N2+S2*VV;
        VV := B.W[I];
        D0 := D0+S0*VV;
        D1 := D1+S1*VV;
        D2 := D2+S2*VV;
        Inc(I);
    end;
    F := B.SY*N0/D0;
    DF := B.SY*(N1*D0-N0*D1)/AP_Sqr(D0);
    D2F := B.SY*((N2*D0-N0*D2)*AP_Sqr(D0)-(N1*D0-N0*D1)*2*D0*D1)/AP_Sqr(AP_Sqr(D0));
end;


(*************************************************************************
This subroutine performs linear transformation of the argument.

INPUT PARAMETERS:
    B       -   rational interpolant in barycentric form
    CA, CB  -   transformation coefficients: x = CA*t + CB

OUTPUT PARAMETERS:
    B       -   transformed interpolant with X replaced by T

  -- ALGLIB PROJECT --
     Copyright 19.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricLinTransX(var B : BarycentricInterpolant;
     CA : Double;
     CB : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    
    //
    // special case, replace by constant F(CB)
    //
    if AP_FP_Eq(CA,0) then
    begin
        B.SY := BarycentricCalc(B, CB);
        V := 1;
        I:=0;
        while I<=B.N-1 do
        begin
            B.Y[I] := 1;
            B.W[I] := V;
            V := -V;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // general case: CA<>0
    //
    I:=0;
    while I<=B.N-1 do
    begin
        B.X[I] := (B.X[I]-CB)/CA;
        Inc(I);
    end;
    if AP_FP_Less(CA,0) then
    begin
        I:=0;
        while I<=B.N-1 do
        begin
            if I<B.N-1-I then
            begin
                J := B.N-1-I;
                V := B.X[I];
                B.X[I] := B.X[J];
                B.X[J] := V;
                V := B.Y[I];
                B.Y[I] := B.Y[J];
                B.Y[J] := V;
                V := B.W[I];
                B.W[I] := B.W[J];
                B.W[J] := V;
            end
            else
            begin
                Break;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
This  subroutine   performs   linear  transformation  of  the  barycentric
interpolant.

INPUT PARAMETERS:
    B       -   rational interpolant in barycentric form
    CA, CB  -   transformation coefficients: B2(x) = CA*B(x) + CB

OUTPUT PARAMETERS:
    B       -   transformed interpolant

  -- ALGLIB PROJECT --
     Copyright 19.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricLinTransY(var B : BarycentricInterpolant;
     CA : Double;
     CB : Double);
var
    I : AlglibInteger;
    V : Double;
begin
    I:=0;
    while I<=B.N-1 do
    begin
        B.Y[I] := CA*B.SY*B.Y[I]+CB;
        Inc(I);
    end;
    B.SY := 0;
    I:=0;
    while I<=B.N-1 do
    begin
        B.SY := Max(B.SY, AbsReal(B.Y[I]));
        Inc(I);
    end;
    if AP_FP_Greater(B.SY,0) then
    begin
        V := 1/B.SY;
        APVMul(@B.Y[0], 0, B.N-1, V);
    end;
end;


(*************************************************************************
Extracts X/Y/W arrays from rational interpolant

INPUT PARAMETERS:
    B   -   barycentric interpolant

OUTPUT PARAMETERS:
    N   -   nodes count, N>0
    X   -   interpolation nodes, array[0..N-1]
    F   -   function values, array[0..N-1]
    W   -   barycentric weights, array[0..N-1]

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricUnpack(const B : BarycentricInterpolant;
     var N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray;
     var W : TReal1DArray);
var
    V : Double;
begin
    N := B.N;
    SetLength(X, N);
    SetLength(Y, N);
    SetLength(W, N);
    V := B.SY;
    APVMove(@X[0], 0, N-1, @B.X[0], 0, N-1);
    APVMove(@Y[0], 0, N-1, @B.Y[0], 0, N-1, V);
    APVMove(@W[0], 0, N-1, @B.W[0], 0, N-1);
end;


(*************************************************************************
Serialization of the barycentric interpolant

INPUT PARAMETERS:
    B   -   barycentric interpolant

OUTPUT PARAMETERS:
    RA      -   array of real numbers which contains interpolant,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricSerialize(const B : BarycentricInterpolant;
     var RA : TReal1DArray;
     var RALen : AlglibInteger);
begin
    RALen := 2+2+3*B.N;
    SetLength(RA, RALen);
    RA[0] := RALen;
    RA[1] := BRCVNum;
    RA[2] := B.N;
    RA[3] := B.SY;
    APVMove(@RA[0], 4, 4+B.N-1, @B.X[0], 0, B.N-1);
    APVMove(@RA[0], 4+B.N, 4+2*B.N-1, @B.Y[0], 0, B.N-1);
    APVMove(@RA[0], 4+2*B.N, 4+3*B.N-1, @B.W[0], 0, B.N-1);
end;


(*************************************************************************
Unserialization of the barycentric interpolant

INPUT PARAMETERS:
    RA  -   array of real numbers which contains interpolant,

OUTPUT PARAMETERS:
    B   -   barycentric interpolant

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricUnserialize(const RA : TReal1DArray;
     var B : BarycentricInterpolant);
begin
    Assert(Round(RA[1])=BRCVNum, 'BarycentricUnserialize: corrupted array!');
    B.N := Round(RA[2]);
    B.SY := RA[3];
    SetLength(B.X, B.N);
    SetLength(B.Y, B.N);
    SetLength(B.W, B.N);
    APVMove(@B.X[0], 0, B.N-1, @RA[0], 4, 4+B.N-1);
    APVMove(@B.Y[0], 0, B.N-1, @RA[0], 4+B.N, 4+2*B.N-1);
    APVMove(@B.W[0], 0, B.N-1, @RA[0], 4+2*B.N, 4+3*B.N-1);
end;


(*************************************************************************
Copying of the barycentric interpolant

INPUT PARAMETERS:
    B   -   barycentric interpolant

OUTPUT PARAMETERS:
    B2  -   copy(B1)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricCopy(const B : BarycentricInterpolant;
     var B2 : BarycentricInterpolant);
begin
    B2.N := B.N;
    B2.SY := B.SY;
    SetLength(B2.X, B2.N);
    SetLength(B2.Y, B2.N);
    SetLength(B2.W, B2.N);
    APVMove(@B2.X[0], 0, B2.N-1, @B.X[0], 0, B2.N-1);
    APVMove(@B2.Y[0], 0, B2.N-1, @B.Y[0], 0, B2.N-1);
    APVMove(@B2.W[0], 0, B2.N-1, @B.W[0], 0, B2.N-1);
end;


(*************************************************************************
Rational interpolant from X/Y/W arrays

F(t) = SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))

INPUT PARAMETERS:
    X   -   interpolation nodes, array[0..N-1]
    F   -   function values, array[0..N-1]
    W   -   barycentric weights, array[0..N-1]
    N   -   nodes count, N>0

OUTPUT PARAMETERS:
    B   -   barycentric interpolant built from (X, Y, W)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricBuildXYW(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     var B : BarycentricInterpolant);
begin
    Assert(N>0, 'BarycentricBuildXYW: incorrect N!');
    
    //
    // fill X/Y/W
    //
    SetLength(B.X, N);
    SetLength(B.Y, N);
    SetLength(B.W, N);
    APVMove(@B.X[0], 0, N-1, @X[0], 0, N-1);
    APVMove(@B.Y[0], 0, N-1, @Y[0], 0, N-1);
    APVMove(@B.W[0], 0, N-1, @W[0], 0, N-1);
    B.N := N;
    
    //
    // Normalize
    //
    BarycentricNormalize(B);
end;


(*************************************************************************
Rational interpolant without poles

The subroutine constructs the rational interpolating function without real
poles  (see  'Barycentric rational interpolation with no  poles  and  high
rates of approximation', Michael S. Floater. and  Kai  Hormann,  for  more
information on this subject).

Input parameters:
    X   -   interpolation nodes, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of nodes, N>0.
    D   -   order of the interpolation scheme, 0 <= D <= N-1.
            D<0 will cause an error.
            D>=N it will be replaced with D=N-1.
            if you don't know what D to choose, use small value about 3-5.

Output parameters:
    B   -   barycentric interpolant.

Note:
    this algorithm always succeeds and calculates the weights  with  close
    to machine precision.

  -- ALGLIB PROJECT --
     Copyright 17.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricBuildFloaterHormann(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     D : AlglibInteger;
     var B : BarycentricInterpolant);
var
    S0 : Double;
    S : Double;
    V : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    Perm : TInteger1DArray;
    WTemp : TReal1DArray;
begin
    Assert(N>0, 'BarycentricFloaterHormann: N<=0!');
    Assert(D>=0, 'BarycentricFloaterHormann: incorrect D!');
    
    //
    // Prepare
    //
    if D>N-1 then
    begin
        D := N-1;
    end;
    B.N := N;
    
    //
    // special case: N=1
    //
    if N=1 then
    begin
        SetLength(B.X, N);
        SetLength(B.Y, N);
        SetLength(B.W, N);
        B.X[0] := X[0];
        B.Y[0] := Y[0];
        B.W[0] := 1;
        BarycentricNormalize(B);
        Exit;
    end;
    
    //
    // Fill X/Y
    //
    SetLength(B.X, N);
    SetLength(B.Y, N);
    APVMove(@B.X[0], 0, N-1, @X[0], 0, N-1);
    APVMove(@B.Y[0], 0, N-1, @Y[0], 0, N-1);
    TagSortFastR(B.X, B.Y, N);
    
    //
    // Calculate Wk
    //
    SetLength(B.W, N);
    S0 := 1;
    K:=1;
    while K<=D do
    begin
        S0 := -S0;
        Inc(K);
    end;
    K:=0;
    while K<=N-1 do
    begin
        
        //
        // Wk
        //
        S := 0;
        I:=Max(K-D, 0);
        while I<=Min(K, N-1-D) do
        begin
            V := 1;
            J:=I;
            while J<=I+D do
            begin
                if J<>K then
                begin
                    V := V/AbsReal(B.X[K]-B.X[J]);
                end;
                Inc(J);
            end;
            S := S+V;
            Inc(I);
        end;
        B.W[K] := S0*S;
        
        //
        // Next S0
        //
        S0 := -S0;
        Inc(K);
    end;
    
    //
    // Normalize
    //
    BarycentricNormalize(B);
end;


(*************************************************************************
Weghted rational least  squares  fitting  using  Floater-Hormann  rational
functions  with  optimal  D  chosen  from  [0,9],  with  constraints   and
individual weights.

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal D (least WEIGHTED root
mean square error) is chosen.  Task  is  linear,  so  linear least squares
solver  is  used.  Complexity  of  this  computational  scheme is O(N*M^2)
(mostly dominated by the least squares solver).

SEE ALSO
* BarycentricFitFloaterHormann(), "lightweight" fitting without invididual
  weights and constraints.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where function values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained barycentric interpolants:
* excessive  constraints  can  be  inconsistent.   Floater-Hormann   basis
  functions aren't as flexible as splines (although they are very smooth).
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints on the function  VALUES at the interval
  boundaries. Note that consustency of the  constraints  on  the  function
  DERIVATIVES is NOT guaranteed (you can use in such cases  cubic  splines
  which are more flexible).
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricFitFloaterHormannWC(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var B : BarycentricInterpolant;
     var Rep : BarycentricFitReport);
var
    D : AlglibInteger;
    I : AlglibInteger;
    WRMSCur : Double;
    WRMSBest : Double;
    LocB : BarycentricInterpolant;
    LocRep : BarycentricFitReport;
    LocInfo : AlglibInteger;
begin
    if (N<1) or (M<2) or (K<0) or (K>=M) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // Find optimal D
    //
    // Info is -3 by default (degenerate constraints).
    // If LocInfo will always be equal to -3, Info will remain equal to -3.
    // If at least once LocInfo will be -4, Info will be -4.
    //
    WRMSBest := MaxRealNumber;
    Rep.DBest := -1;
    Info := -3;
    D:=0;
    while D<=Min(9, N-1) do
    begin
        BarycentricFitWCFixedD(X, Y, W, N, XC, YC, DC, K, M, D, LocInfo, LocB, LocRep);
        Assert((LocInfo=-4) or (LocInfo=-3) or (LocInfo>0), 'BarycentricFitFloaterHormannWC: unexpected result from BarycentricFitWCFixedD!');
        if LocInfo>0 then
        begin
            
            //
            // Calculate weghted RMS
            //
            WRMSCur := 0;
            I:=0;
            while I<=N-1 do
            begin
                WRMSCur := WRMSCur+AP_Sqr(W[I]*(Y[I]-BarycentricCalc(LocB, X[I])));
                Inc(I);
            end;
            WRMSCur := Sqrt(WRMSCur/N);
            if AP_FP_Less(WRMSCur,WRMSBest) or (Rep.DBest<0) then
            begin
                BarycentricCopy(LocB, B);
                Rep.DBest := D;
                Info := 1;
                Rep.RMSError := LocRep.RMSError;
                Rep.AvgError := LocRep.AvgError;
                Rep.AvgRelError := LocRep.AvgRelError;
                Rep.MaxError := LocRep.MaxError;
                Rep.TaskRCond := LocRep.TaskRCond;
                WRMSBest := WRMSCur;
            end;
        end
        else
        begin
            if (LocInfo<>-3) and (Info<0) then
            begin
                Info := LocInfo;
            end;
        end;
        Inc(D);
    end;
end;


(*************************************************************************
Rational least squares fitting, without weights and constraints.

See BarycentricFitFloaterHormannWC() for more information.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricFitFloaterHormann(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var B : BarycentricInterpolant;
     var Rep : BarycentricFitReport);
var
    W : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    I : AlglibInteger;
begin
    if N<1 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(W, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := 1;
        Inc(I);
    end;
    BarycentricFitFloaterHormannWC(X, Y, W, N, XC, YC, DC, 0, M, Info, B, Rep);
end;


(*************************************************************************
Normalization of barycentric interpolant:
* B.N, B.X, B.Y and B.W are initialized
* B.SY is NOT initialized
* Y[] is normalized, scaling coefficient is stored in B.SY
* W[] is normalized, no scaling coefficient is stored
* X[] is sorted

Internal subroutine.
*************************************************************************)
procedure BarycentricNormalize(var B : BarycentricInterpolant);
var
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    J2 : AlglibInteger;
    V : Double;
begin
    
    //
    // Normalize task: |Y|<=1, |W|<=1, sort X[]
    //
    B.SY := 0;
    I:=0;
    while I<=B.N-1 do
    begin
        B.SY := Max(B.SY, AbsReal(B.Y[I]));
        Inc(I);
    end;
    if AP_FP_Greater(B.SY,0) and AP_FP_Greater(AbsReal(B.SY-1),10*MachineEpsilon) then
    begin
        V := 1/B.SY;
        APVMul(@B.Y[0], 0, B.N-1, V);
    end;
    V := 0;
    I:=0;
    while I<=B.N-1 do
    begin
        V := Max(V, AbsReal(B.W[I]));
        Inc(I);
    end;
    if AP_FP_Greater(V,0) and AP_FP_Greater(AbsReal(V-1),10*MachineEpsilon) then
    begin
        V := 1/V;
        APVMul(@B.W[0], 0, B.N-1, V);
    end;
    I:=0;
    while I<=B.N-2 do
    begin
        if AP_FP_Less(B.X[I+1],B.X[I]) then
        begin
            TagSort(B.X, B.N, P1, P2);
            J:=0;
            while J<=B.N-1 do
            begin
                J2 := P2[J];
                V := B.Y[J];
                B.Y[J] := B.Y[J2];
                B.Y[J2] := V;
                V := B.W[J];
                B.W[J] := B.W[J2];
                B.W[J2] := V;
                Inc(J);
            end;
            Break;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal subroutine, calculates barycentric basis functions.
Used for efficient simultaneous calculation of N basis functions.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure BarycentricCalcBasis(const B : BarycentricInterpolant;
     T : Double;
     var Y : TReal1DArray);
var
    S2 : Double;
    S : Double;
    V : Double;
    I : AlglibInteger;
    J : AlglibInteger;
begin
    
    //
    // special case: N=1
    //
    if B.N=1 then
    begin
        Y[0] := 1;
        Exit;
    end;
    
    //
    // Here we assume that task is normalized, i.e.:
    // 1. abs(Y[i])<=1
    // 2. abs(W[i])<=1
    // 3. X[] is ordered
    //
    // First, we decide: should we use "safe" formula (guarded
    // against overflow) or fast one?
    //
    S := AbsReal(T-B.X[0]);
    I:=0;
    while I<=B.N-1 do
    begin
        V := B.X[I];
        if AP_FP_Eq(V,T) then
        begin
            J:=0;
            while J<=B.N-1 do
            begin
                Y[J] := 0;
                Inc(J);
            end;
            Y[I] := 1;
            Exit;
        end;
        V := AbsReal(T-V);
        if AP_FP_Less(V,S) then
        begin
            S := V;
        end;
        Inc(I);
    end;
    S2 := 0;
    I:=0;
    while I<=B.N-1 do
    begin
        V := S/(T-B.X[I]);
        V := V*B.W[I];
        Y[I] := V;
        S2 := S2+V;
        Inc(I);
    end;
    V := 1/S2;
    APVMul(@Y[0], 0, B.N-1, V);
end;


(*************************************************************************
Internal Floater-Hormann fitting subroutine for fixed D
*************************************************************************)
procedure BarycentricFitWCFixedD(X : TReal1DArray;
     Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     XC : TReal1DArray;
     YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     D : AlglibInteger;
     var Info : AlglibInteger;
     var B : BarycentricInterpolant;
     var Rep : BarycentricFitReport);
var
    FMatrix : TReal2DArray;
    CMatrix : TReal2DArray;
    Y2 : TReal1DArray;
    W2 : TReal1DArray;
    SX : TReal1DArray;
    SY : TReal1DArray;
    SBF : TReal1DArray;
    XOriginal : TReal1DArray;
    YOriginal : TReal1DArray;
    Tmp : TReal1DArray;
    LRep : LSFitReport;
    V0 : Double;
    V1 : Double;
    MX : Double;
    B2 : BarycentricInterpolant;
    I : AlglibInteger;
    J : AlglibInteger;
    RelCnt : AlglibInteger;
    XA : Double;
    XB : Double;
    SA : Double;
    SB : Double;
    Decay : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    XC := DynamicArrayCopy(XC);
    YC := DynamicArrayCopy(YC);
    if (N<1) or (M<2) or (K<0) or (K>=M) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=K-1 do
    begin
        Info := 0;
        if DC[I]<0 then
        begin
            Info := -1;
        end;
        if DC[I]>1 then
        begin
            Info := -1;
        end;
        if Info<0 then
        begin
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // weight decay for correct handling of task which becomes
    // degenerate after constraints are applied
    //
    Decay := 10000*MachineEpsilon;
    
    //
    // Scale X, Y, XC, YC
    //
    LSFitScaleXY(X, Y, N, XC, YC, DC, K, XA, XB, SA, SB, XOriginal, YOriginal);
    
    //
    // allocate space, initialize:
    // * FMatrix-   values of basis functions at X[]
    // * CMatrix-   values (derivatives) of basis functions at XC[]
    //
    SetLength(Y2, N+M);
    SetLength(W2, N+M);
    SetLength(FMatrix, N+M, M);
    if K>0 then
    begin
        SetLength(CMatrix, K, M+1);
    end;
    SetLength(Y2, N+M);
    SetLength(W2, N+M);
    
    //
    // Prepare design and constraints matrices:
    // * fill constraints matrix
    // * fill first N rows of design matrix with values
    // * fill next M rows of design matrix with regularizing term
    // * append M zeros to Y
    // * append M elements, mean(abs(W)) each, to W
    //
    SetLength(SX, M);
    SetLength(SY, M);
    SetLength(SBF, M);
    J:=0;
    while J<=M-1 do
    begin
        SX[J] := AP_Double(2*J)/(M-1)-1;
        Inc(J);
    end;
    I:=0;
    while I<=M-1 do
    begin
        SY[I] := 1;
        Inc(I);
    end;
    BarycentricBuildFloaterHormann(SX, SY, M, D, B2);
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        BarycentricCalcBasis(B2, X[I], SBF);
        APVMove(@FMatrix[I][0], 0, M-1, @SBF[0], 0, M-1);
        Y2[I] := Y[I];
        W2[I] := W[I];
        MX := MX+AbsReal(W[I])/N;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=M-1 do
        begin
            if I=J then
            begin
                FMatrix[N+I,J] := Decay;
            end
            else
            begin
                FMatrix[N+I,J] := 0;
            end;
            Inc(J);
        end;
        Y2[N+I] := 0;
        W2[N+I] := MX;
        Inc(I);
    end;
    if K>0 then
    begin
        J:=0;
        while J<=M-1 do
        begin
            I:=0;
            while I<=M-1 do
            begin
                SY[I] := 0;
                Inc(I);
            end;
            SY[J] := 1;
            BarycentricBuildFloaterHormann(SX, SY, M, D, B2);
            I:=0;
            while I<=K-1 do
            begin
                Assert((DC[I]>=0) and (DC[I]<=1), 'BarycentricFit: internal error!');
                BarycentricDiff1(B2, XC[I], V0, V1);
                if DC[I]=0 then
                begin
                    CMatrix[I,J] := V0;
                end;
                if DC[I]=1 then
                begin
                    CMatrix[I,J] := V1;
                end;
                Inc(I);
            end;
            Inc(J);
        end;
        I:=0;
        while I<=K-1 do
        begin
            CMatrix[I,M] := YC[I];
            Inc(I);
        end;
    end;
    
    //
    // Solve constrained task
    //
    if K>0 then
    begin
        
        //
        // solve using regularization
        //
        LSFitLinearWC(Y2, W2, FMatrix, CMatrix, N+M, M, K, Info, Tmp, LRep);
    end
    else
    begin
        
        //
        // no constraints, no regularization needed
        //
        LSFitLinearWC(Y, W, FMatrix, CMatrix, N, M, K, Info, Tmp, LRep);
    end;
    if Info<0 then
    begin
        Exit;
    end;
    
    //
    // Generate interpolant and scale it
    //
    APVMove(@SY[0], 0, M-1, @Tmp[0], 0, M-1);
    BarycentricBuildFloaterHormann(SX, SY, M, D, B);
    BarycentricLinTransX(B, 2/(XB-XA), -(XA+XB)/(XB-XA));
    BarycentricLinTransY(B, SB-SA, SA);
    
    //
    // Scale absolute errors obtained from LSFitLinearW.
    // Relative error should be calculated separately
    // (because of shifting/scaling of the task)
    //
    Rep.TaskRCond := LRep.TaskRCond;
    Rep.RMSError := LRep.RMSError*(SB-SA);
    Rep.AvgError := LRep.AvgError*(SB-SA);
    Rep.MaxError := LRep.MaxError*(SB-SA);
    Rep.AvgRelError := 0;
    RelCnt := 0;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Neq(YOriginal[I],0) then
        begin
            Rep.AvgRelError := Rep.AvgRelError+AbsReal(BarycentricCalc(B, XOriginal[I])-YOriginal[I])/AbsReal(YOriginal[I]);
            RelCnt := RelCnt+1;
        end;
        Inc(I);
    end;
    if RelCnt<>0 then
    begin
        Rep.AvgRelError := Rep.AvgRelError/RelCnt;
    end;
end;


end.