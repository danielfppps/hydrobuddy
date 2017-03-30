{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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
unit polint;
interface
uses Math, Sysutils, Ap, tsort, ratinterpolation, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, ratint, apserv;

type
(*************************************************************************
Polynomial fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************)
PolynomialFitReport = record
    TaskRCond : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    MaxError : Double;
end;



procedure PolynomialBuild(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
procedure PolynomialBuildEqDist(A : Double;
     B : Double;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
procedure PolynomialBuildCheb1(A : Double;
     B : Double;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
procedure PolynomialBuildCheb2(A : Double;
     B : Double;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
function PolynomialCalcEqDist(A : Double;
     B : Double;
     const F : TReal1DArray;
     N : AlglibInteger;
     T : Double):Double;
function PolynomialCalcCheb1(A : Double;
     B : Double;
     const F : TReal1DArray;
     N : AlglibInteger;
     T : Double):Double;
function PolynomialCalcCheb2(A : Double;
     B : Double;
     const F : TReal1DArray;
     N : AlglibInteger;
     T : Double):Double;
procedure PolynomialFit(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var P : BarycentricInterpolant;
     var Rep : PolynomialFitReport);
procedure PolynomialFitWC(X : TReal1DArray;
     Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     XC : TReal1DArray;
     YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var P : BarycentricInterpolant;
     var Rep : PolynomialFitReport);

implementation

(*************************************************************************
Lagrange intepolant: generation of the model on the general grid.
This function has O(N^2) complexity.

INPUT PARAMETERS:
    X   -   abscissas, array[0..N-1]
    Y   -   function values, array[0..N-1]
    N   -   number of points, N>=1

OIYTPUT PARAMETERS
    P   -   barycentric model which represents Lagrange interpolant
            (see ratint unit info and BarycentricCalc() description for
            more information).

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure PolynomialBuild(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
var
    J : AlglibInteger;
    K : AlglibInteger;
    W : TReal1DArray;
    B : Double;
    A : Double;
    V : Double;
    MX : Double;
begin
    Assert(N>0, 'PolIntBuild: N<=0!');
    
    //
    // calculate W[j]
    // multi-pass algorithm is used to avoid overflow
    //
    SetLength(W, N);
    A := X[0];
    B := X[0];
    J:=0;
    while J<=N-1 do
    begin
        W[J] := 1;
        A := Min(A, X[J]);
        B := Max(B, X[J]);
        Inc(J);
    end;
    K:=0;
    while K<=N-1 do
    begin
        
        //
        // W[K] is used instead of 0.0 because
        // cycle on J does not touch K-th element
        // and we MUST get maximum from ALL elements
        //
        MX := AbsReal(W[K]);
        J:=0;
        while J<=N-1 do
        begin
            if J<>K then
            begin
                V := (B-A)/(X[J]-X[K]);
                W[J] := W[J]*V;
                MX := Max(MX, AbsReal(W[J]));
            end;
            Inc(J);
        end;
        if K mod 5=0 then
        begin
            
            //
            // every 5-th run we renormalize W[]
            //
            V := 1/MX;
            APVMul(@W[0], 0, N-1, V);
        end;
        Inc(K);
    end;
    BarycentricBuildXYW(X, Y, W, N, P);
end;


(*************************************************************************
Lagrange intepolant: generation of the model on equidistant grid.
This function has O(N) complexity.

INPUT PARAMETERS:
    A   -   left boundary of [A,B]
    B   -   right boundary of [A,B]
    Y   -   function values at the nodes, array[0..N-1]
    N   -   number of points, N>=1
            for N=1 a constant model is constructed.

OIYTPUT PARAMETERS
    P   -   barycentric model which represents Lagrange interpolant
            (see ratint unit info and BarycentricCalc() description for
            more information).

  -- ALGLIB --
     Copyright 03.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure PolynomialBuildEqDist(A : Double;
     B : Double;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
var
    I : AlglibInteger;
    W : TReal1DArray;
    X : TReal1DArray;
    V : Double;
begin
    Assert(N>0, 'PolIntBuildEqDist: N<=0!');
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        SetLength(X, 1);
        SetLength(W, 1);
        X[0] := Double(0.5)*(B+A);
        W[0] := 1;
        BarycentricBuildXYW(X, Y, W, 1, P);
        Exit;
    end;
    
    //
    // general case
    //
    SetLength(X, N);
    SetLength(W, N);
    V := 1;
    I:=0;
    while I<=N-1 do
    begin
        W[I] := V;
        X[I] := A+(B-A)*I/(N-1);
        V := -V*(N-1-I);
        V := V/(I+1);
        Inc(I);
    end;
    BarycentricBuildXYW(X, Y, W, N, P);
end;


(*************************************************************************
Lagrange intepolant on Chebyshev grid (first kind).
This function has O(N) complexity.

INPUT PARAMETERS:
    A   -   left boundary of [A,B]
    B   -   right boundary of [A,B]
    Y   -   function values at the nodes, array[0..N-1],
            Y[I] = Y(0.5*(B+A) + 0.5*(B-A)*Cos(PI*(2*i+1)/(2*n)))
    N   -   number of points, N>=1
            for N=1 a constant model is constructed.

OIYTPUT PARAMETERS
    P   -   barycentric model which represents Lagrange interpolant
            (see ratint unit info and BarycentricCalc() description for
            more information).

  -- ALGLIB --
     Copyright 03.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure PolynomialBuildCheb1(A : Double;
     B : Double;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
var
    I : AlglibInteger;
    W : TReal1DArray;
    X : TReal1DArray;
    V : Double;
    T : Double;
begin
    Assert(N>0, 'PolIntBuildCheb1: N<=0!');
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        SetLength(X, 1);
        SetLength(W, 1);
        X[0] := Double(0.5)*(B+A);
        W[0] := 1;
        BarycentricBuildXYW(X, Y, W, 1, P);
        Exit;
    end;
    
    //
    // general case
    //
    SetLength(X, N);
    SetLength(W, N);
    V := 1;
    I:=0;
    while I<=N-1 do
    begin
        T := Tan(Double(0.5)*Pi*(2*i+1)/(2*n));
        W[I] := 2*V*T/(1+AP_Sqr(T));
        X[I] := Double(0.5)*(B+A)+Double(0.5)*(B-A)*(1-AP_Sqr(T))/(1+AP_Sqr(T));
        V := -V;
        Inc(I);
    end;
    BarycentricBuildXYW(X, Y, W, N, P);
end;


(*************************************************************************
Lagrange intepolant on Chebyshev grid (second kind).
This function has O(N) complexity.

INPUT PARAMETERS:
    A   -   left boundary of [A,B]
    B   -   right boundary of [A,B]
    Y   -   function values at the nodes, array[0..N-1],
            Y[I] = Y(0.5*(B+A) + 0.5*(B-A)*Cos(PI*i/(n-1)))
    N   -   number of points, N>=1
            for N=1 a constant model is constructed.

OIYTPUT PARAMETERS
    P   -   barycentric model which represents Lagrange interpolant
            (see ratint unit info and BarycentricCalc() description for
            more information).

  -- ALGLIB --
     Copyright 03.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure PolynomialBuildCheb2(A : Double;
     B : Double;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var P : BarycentricInterpolant);
var
    I : AlglibInteger;
    W : TReal1DArray;
    X : TReal1DArray;
    V : Double;
begin
    Assert(N>0, 'PolIntBuildCheb2: N<=0!');
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        SetLength(X, 1);
        SetLength(W, 1);
        X[0] := Double(0.5)*(B+A);
        W[0] := 1;
        BarycentricBuildXYW(X, Y, W, 1, P);
        Exit;
    end;
    
    //
    // general case
    //
    SetLength(X, N);
    SetLength(W, N);
    V := 1;
    I:=0;
    while I<=N-1 do
    begin
        if (I=0) or (I=N-1) then
        begin
            W[I] := V*Double(0.5);
        end
        else
        begin
            W[I] := V;
        end;
        X[I] := Double(0.5)*(B+A)+Double(0.5)*(B-A)*Cos(Pi*i/(n-1));
        V := -V;
        Inc(I);
    end;
    BarycentricBuildXYW(X, Y, W, N, P);
end;


(*************************************************************************
Fast equidistant polynomial interpolation function with O(N) complexity

INPUT PARAMETERS:
    A   -   left boundary of [A,B]
    B   -   right boundary of [A,B]
    F   -   function values, array[0..N-1]
    N   -   number of points on equidistant grid, N>=1
            for N=1 a constant model is constructed.
    T   -   position where P(x) is calculated

RESULT
    value of the Lagrange interpolant at T
    
IMPORTANT
    this function provides fast interface which is not overflow-safe
    nor it is very precise.
    the best option is to use  PolynomialBuildEqDist()/BarycentricCalc()
    subroutines unless you are pretty sure that your data will not result
    in overflow.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function PolynomialCalcEqDist(A : Double;
     B : Double;
     const F : TReal1DArray;
     N : AlglibInteger;
     T : Double):Double;
var
    S1 : Double;
    S2 : Double;
    V : Double;
    Threshold : Double;
    S : Double;
    H : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    W : Double;
    X : Double;
begin
    Assert(N>0, 'PolIntEqDist: N<=0!');
    Threshold := Sqrt(MinRealNumber);
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        Result := F[0];
        Exit;
    end;
    
    //
    // First, decide: should we use "safe" formula (guarded
    // against overflow) or fast one?
    //
    J := 0;
    S := T-A;
    I:=1;
    while I<=N-1 do
    begin
        X := A+AP_Double(I)/(N-1)*(B-A);
        if AP_FP_Less(AbsReal(T-X),AbsReal(S)) then
        begin
            S := T-X;
            J := I;
        end;
        Inc(I);
    end;
    if AP_FP_Eq(S,0) then
    begin
        Result := F[J];
        Exit;
    end;
    if AP_FP_Greater(AbsReal(S),Threshold) then
    begin
        
        //
        // use fast formula
        //
        J := -1;
        S := Double(1.0);
    end;
    
    //
    // Calculate using safe or fast barycentric formula
    //
    S1 := 0;
    S2 := 0;
    W := Double(1.0);
    H := (B-A)/(N-1);
    I:=0;
    while I<=N-1 do
    begin
        if I<>J then
        begin
            V := S*W/(T-(A+I*H));
            S1 := S1+V*F[I];
            S2 := S2+V;
        end
        else
        begin
            V := W;
            S1 := S1+V*F[I];
            S2 := S2+V;
        end;
        W := -W*(N-1-I);
        W := W/(I+1);
        Inc(I);
    end;
    Result := S1/S2;
end;


(*************************************************************************
Fast polynomial interpolation function on Chebyshev points (first kind)
with O(N) complexity.

INPUT PARAMETERS:
    A   -   left boundary of [A,B]
    B   -   right boundary of [A,B]
    F   -   function values, array[0..N-1]
    N   -   number of points on Chebyshev grid (first kind),
            X[i] = 0.5*(B+A) + 0.5*(B-A)*Cos(PI*(2*i+1)/(2*n))
            for N=1 a constant model is constructed.
    T   -   position where P(x) is calculated

RESULT
    value of the Lagrange interpolant at T

IMPORTANT
    this function provides fast interface which is not overflow-safe
    nor it is very precise.
    the best option is to use  PolIntBuildCheb1()/BarycentricCalc()
    subroutines unless you are pretty sure that your data will not result
    in overflow.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function PolynomialCalcCheb1(A : Double;
     B : Double;
     const F : TReal1DArray;
     N : AlglibInteger;
     T : Double):Double;
var
    S1 : Double;
    S2 : Double;
    V : Double;
    Threshold : Double;
    S : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    A0 : Double;
    Delta : Double;
    Alpha : Double;
    Beta : Double;
    CA : Double;
    SA : Double;
    TempC : Double;
    TempS : Double;
    X : Double;
    W : Double;
    P1 : Double;
begin
    Assert(N>0, 'PolIntCheb1: N<=0!');
    Threshold := Sqrt(MinRealNumber);
    T := (T-Double(0.5)*(A+B))/(Double(0.5)*(B-A));
    
    //
    // Fast exit
    //
    if N=1 then
    begin
        Result := F[0];
        Exit;
    end;
    
    //
    // Prepare information for the recurrence formula
    // used to calculate sin(pi*(2j+1)/(2n+2)) and
    // cos(pi*(2j+1)/(2n+2)):
    //
    // A0    = pi/(2n+2)
    // Delta = pi/(n+1)
    // Alpha = 2 sin^2 (Delta/2)
    // Beta  = sin(Delta)
    //
    // so that sin(..) = sin(A0+j*delta) and cos(..) = cos(A0+j*delta).
    // Then we use
    //
    // sin(x+delta) = sin(x) - (alpha*sin(x) - beta*cos(x))
    // cos(x+delta) = cos(x) - (alpha*cos(x) - beta*sin(x))
    //
    // to repeatedly calculate sin(..) and cos(..).
    //
    A0 := Pi/(2*(N-1)+2);
    Delta := 2*Pi/(2*(N-1)+2);
    Alpha := 2*AP_Sqr(Sin(Delta/2));
    Beta := Sin(Delta);
    
    //
    // First, decide: should we use "safe" formula (guarded
    // against overflow) or fast one?
    //
    CA := Cos(A0);
    SA := Sin(A0);
    J := 0;
    X := CA;
    S := T-X;
    I:=1;
    while I<=N-1 do
    begin
        
        //
        // Next X[i]
        //
        TempS := SA-(Alpha*SA-Beta*CA);
        TempC := CA-(Alpha*CA+Beta*SA);
        SA := TempS;
        CA := TempC;
        X := CA;
        
        //
        // Use X[i]
        //
        if AP_FP_Less(AbsReal(T-X),AbsReal(S)) then
        begin
            S := T-X;
            J := I;
        end;
        Inc(I);
    end;
    if AP_FP_Eq(S,0) then
    begin
        Result := F[J];
        Exit;
    end;
    if AP_FP_Greater(AbsReal(S),Threshold) then
    begin
        
        //
        // use fast formula
        //
        J := -1;
        S := Double(1.0);
    end;
    
    //
    // Calculate using safe or fast barycentric formula
    //
    S1 := 0;
    S2 := 0;
    CA := Cos(A0);
    SA := Sin(A0);
    P1 := Double(1.0);
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // Calculate X[i], W[i]
        //
        X := CA;
        W := P1*SA;
        
        //
        // Proceed
        //
        if I<>J then
        begin
            V := S*W/(T-X);
            S1 := S1+V*F[I];
            S2 := S2+V;
        end
        else
        begin
            V := W;
            S1 := S1+V*F[I];
            S2 := S2+V;
        end;
        
        //
        // Next CA, SA, P1
        //
        TempS := SA-(Alpha*SA-Beta*CA);
        TempC := CA-(Alpha*CA+Beta*SA);
        SA := TempS;
        CA := TempC;
        P1 := -P1;
        Inc(I);
    end;
    Result := S1/S2;
end;


(*************************************************************************
Fast polynomial interpolation function on Chebyshev points (second kind)
with O(N) complexity.

INPUT PARAMETERS:
    A   -   left boundary of [A,B]
    B   -   right boundary of [A,B]
    F   -   function values, array[0..N-1]
    N   -   number of points on Chebyshev grid (second kind),
            X[i] = 0.5*(B+A) + 0.5*(B-A)*Cos(PI*i/(n-1))
            for N=1 a constant model is constructed.
    T   -   position where P(x) is calculated

RESULT
    value of the Lagrange interpolant at T

IMPORTANT
    this function provides fast interface which is not overflow-safe
    nor it is very precise.
    the best option is to use PolIntBuildCheb2()/BarycentricCalc()
    subroutines unless you are pretty sure that your data will not result
    in overflow.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function PolynomialCalcCheb2(A : Double;
     B : Double;
     const F : TReal1DArray;
     N : AlglibInteger;
     T : Double):Double;
var
    S1 : Double;
    S2 : Double;
    V : Double;
    Threshold : Double;
    S : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    A0 : Double;
    Delta : Double;
    Alpha : Double;
    Beta : Double;
    CA : Double;
    SA : Double;
    TempC : Double;
    TempS : Double;
    X : Double;
    W : Double;
    P1 : Double;
begin
    Assert(N>0, 'PolIntCheb2: N<=0!');
    Threshold := Sqrt(MinRealNumber);
    T := (T-Double(0.5)*(A+B))/(Double(0.5)*(B-A));
    
    //
    // Fast exit
    //
    if N=1 then
    begin
        Result := F[0];
        Exit;
    end;
    
    //
    // Prepare information for the recurrence formula
    // used to calculate sin(pi*i/n) and
    // cos(pi*i/n):
    //
    // A0    = 0
    // Delta = pi/n
    // Alpha = 2 sin^2 (Delta/2)
    // Beta  = sin(Delta)
    //
    // so that sin(..) = sin(A0+j*delta) and cos(..) = cos(A0+j*delta).
    // Then we use
    //
    // sin(x+delta) = sin(x) - (alpha*sin(x) - beta*cos(x))
    // cos(x+delta) = cos(x) - (alpha*cos(x) - beta*sin(x))
    //
    // to repeatedly calculate sin(..) and cos(..).
    //
    A0 := Double(0.0);
    Delta := Pi/(N-1);
    Alpha := 2*AP_Sqr(Sin(Delta/2));
    Beta := Sin(Delta);
    
    //
    // First, decide: should we use "safe" formula (guarded
    // against overflow) or fast one?
    //
    CA := Cos(A0);
    SA := Sin(A0);
    J := 0;
    X := CA;
    S := T-X;
    I:=1;
    while I<=N-1 do
    begin
        
        //
        // Next X[i]
        //
        TempS := SA-(Alpha*SA-Beta*CA);
        TempC := CA-(Alpha*CA+Beta*SA);
        SA := TempS;
        CA := TempC;
        X := CA;
        
        //
        // Use X[i]
        //
        if AP_FP_Less(AbsReal(T-X),AbsReal(S)) then
        begin
            S := T-X;
            J := I;
        end;
        Inc(I);
    end;
    if AP_FP_Eq(S,0) then
    begin
        Result := F[J];
        Exit;
    end;
    if AP_FP_Greater(AbsReal(S),Threshold) then
    begin
        
        //
        // use fast formula
        //
        J := -1;
        S := Double(1.0);
    end;
    
    //
    // Calculate using safe or fast barycentric formula
    //
    S1 := 0;
    S2 := 0;
    CA := Cos(A0);
    SA := Sin(A0);
    P1 := Double(1.0);
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // Calculate X[i], W[i]
        //
        X := CA;
        if (I=0) or (I=N-1) then
        begin
            W := Double(0.5)*P1;
        end
        else
        begin
            W := Double(1.0)*P1;
        end;
        
        //
        // Proceed
        //
        if I<>J then
        begin
            V := S*W/(T-X);
            S1 := S1+V*F[I];
            S2 := S2+V;
        end
        else
        begin
            V := W;
            S1 := S1+V*F[I];
            S2 := S2+V;
        end;
        
        //
        // Next CA, SA, P1
        //
        TempS := SA-(Alpha*SA-Beta*CA);
        TempC := CA-(Alpha*CA+Beta*SA);
        SA := TempS;
        CA := TempC;
        P1 := -P1;
        Inc(I);
    end;
    Result := S1/S2;
end;


(*************************************************************************
Least squares fitting by polynomial.

This subroutine is "lightweight" alternative for more complex and feature-
rich PolynomialFitWC().  See  PolynomialFitWC() for more information about
subroutine parameters (we don't duplicate it here because of length)

  -- ALGLIB PROJECT --
     Copyright 12.10.2009 by Bochkanov Sergey
*************************************************************************)
procedure PolynomialFit(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var P : BarycentricInterpolant;
     var Rep : PolynomialFitReport);
var
    I : AlglibInteger;
    W : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
begin
    if N>0 then
    begin
        SetLength(W, N);
        I:=0;
        while I<=N-1 do
        begin
            W[I] := 1;
            Inc(I);
        end;
    end;
    PolynomialFitWC(X, Y, W, N, XC, YC, DC, 0, M, Info, P, Rep);
end;


(*************************************************************************
Weighted  fitting  by  Chebyshev  polynomial  in  barycentric  form,  with
constraints on function values or first derivatives.

Small regularizing term is used when solving constrained tasks (to improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO:
    PolynomialFit()

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where polynomial values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that P(XC[i])=YC[i]
            * DC[i]=1   means that P'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions (= polynomial_degree + 1), M>=1

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    P   -   interpolant in barycentric form.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
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
constrained regression splines:
* even simple constraints can be inconsistent, see  Wikipedia  article  on
  this subject: http://en.wikipedia.org/wiki/Birkhoff_interpolation
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the one special cases, however, we can  guarantee  consistency.  This
  case  is:  M>1  and constraints on the function values (NOT DERIVATIVES)

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 10.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure PolynomialFitWC(X : TReal1DArray;
     Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     XC : TReal1DArray;
     YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var P : BarycentricInterpolant;
     var Rep : PolynomialFitReport);
var
    XA : Double;
    XB : Double;
    SA : Double;
    SB : Double;
    XOriginal : TReal1DArray;
    YOriginal : TReal1DArray;
    Y2 : TReal1DArray;
    W2 : TReal1DArray;
    Tmp : TReal1DArray;
    Tmp2 : TReal1DArray;
    TmpDiff : TReal1DArray;
    BX : TReal1DArray;
    BY : TReal1DArray;
    BW : TReal1DArray;
    FMatrix : TReal2DArray;
    CMatrix : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MX : Double;
    Decay : Double;
    U : Double;
    V : Double;
    S : Double;
    RelCnt : AlglibInteger;
    LRep : LSFitReport;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    XC := DynamicArrayCopy(XC);
    YC := DynamicArrayCopy(YC);
    if (M<1) or (N<1) or (K<0) or (K>=M) then
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
    // allocate space, initialize/fill:
    // * FMatrix-   values of basis functions at X[]
    // * CMatrix-   values (derivatives) of basis functions at XC[]
    // * fill constraints matrix
    // * fill first N rows of design matrix with values
    // * fill next M rows of design matrix with regularizing term
    // * append M zeros to Y
    // * append M elements, mean(abs(W)) each, to W
    //
    SetLength(Y2, N+M);
    SetLength(W2, N+M);
    SetLength(Tmp, M);
    SetLength(TmpDiff, M);
    SetLength(FMatrix, N+M, M);
    if K>0 then
    begin
        SetLength(CMatrix, K, M+1);
    end;
    
    //
    // Fill design matrix, Y2, W2:
    // * first N rows with basis functions for original points
    // * next M rows with decay terms
    //
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // prepare Ith row
        // use Tmp for calculations to avoid multidimensional arrays overhead
        //
        J:=0;
        while J<=M-1 do
        begin
            if J=0 then
            begin
                Tmp[J] := 1;
            end
            else
            begin
                if J=1 then
                begin
                    Tmp[J] := X[I];
                end
                else
                begin
                    Tmp[J] := 2*X[I]*Tmp[J-1]-Tmp[J-2];
                end;
            end;
            Inc(J);
        end;
        APVMove(@FMatrix[I][0], 0, M-1, @Tmp[0], 0, M-1);
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
        Inc(I);
    end;
    APVMove(@Y2[0], 0, N-1, @Y[0], 0, N-1);
    APVMove(@W2[0], 0, N-1, @W[0], 0, N-1);
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        MX := MX+AbsReal(W[I]);
        Inc(I);
    end;
    MX := MX/N;
    I:=0;
    while I<=M-1 do
    begin
        Y2[N+I] := 0;
        W2[N+I] := MX;
        Inc(I);
    end;
    
    //
    // fill constraints matrix
    //
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // prepare Ith row
        // use Tmp for basis function values,
        // TmpDiff for basos function derivatives
        //
        J:=0;
        while J<=M-1 do
        begin
            if J=0 then
            begin
                Tmp[J] := 1;
                TmpDiff[J] := 0;
            end
            else
            begin
                if J=1 then
                begin
                    Tmp[J] := XC[I];
                    TmpDiff[J] := 1;
                end
                else
                begin
                    Tmp[J] := 2*XC[I]*Tmp[J-1]-Tmp[J-2];
                    TmpDiff[J] := 2*(Tmp[J-1]+XC[I]*TmpDiff[J-1])-TmpDiff[J-2];
                end;
            end;
            Inc(J);
        end;
        if DC[I]=0 then
        begin
            APVMove(@CMatrix[I][0], 0, M-1, @Tmp[0], 0, M-1);
        end;
        if DC[I]=1 then
        begin
            APVMove(@CMatrix[I][0], 0, M-1, @TmpDiff[0], 0, M-1);
        end;
        CMatrix[I,M] := YC[I];
        Inc(I);
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
        LSFitLinearWC(Y, W, FMatrix, CMatrix, N, M, 0, Info, Tmp, LRep);
    end;
    if Info<0 then
    begin
        Exit;
    end;
    
    //
    // Generate barycentric model and scale it
    // * BX, BY store barycentric model nodes
    // * FMatrix is reused (remember - it is at least MxM, what we need)
    //
    // Model intialization is done in O(M^2). In principle, it can be
    // done in O(M*log(M)), but before it we solved task with O(N*M^2)
    // complexity, so it is only a small amount of total time spent.
    //
    SetLength(BX, M);
    SetLength(BY, M);
    SetLength(BW, M);
    SetLength(Tmp2, M);
    S := 1;
    I:=0;
    while I<=M-1 do
    begin
        if M<>1 then
        begin
            U := Cos(Pi*I/(M-1));
        end
        else
        begin
            U := 0;
        end;
        V := 0;
        J:=0;
        while J<=M-1 do
        begin
            if J=0 then
            begin
                Tmp2[J] := 1;
            end
            else
            begin
                if J=1 then
                begin
                    Tmp2[J] := U;
                end
                else
                begin
                    Tmp2[J] := 2*U*Tmp2[J-1]-Tmp2[J-2];
                end;
            end;
            V := V+Tmp[J]*Tmp2[J];
            Inc(J);
        end;
        BX[I] := U;
        BY[I] := V;
        BW[I] := S;
        if (I=0) or (I=M-1) then
        begin
            BW[I] := Double(0.5)*BW[I];
        end;
        S := -S;
        Inc(I);
    end;
    BarycentricBuildXYW(BX, BY, BW, M, P);
    BarycentricLinTransX(P, 2/(XB-XA), -(XA+XB)/(XB-XA));
    BarycentricLinTransY(P, SB-SA, SA);
    
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
            Rep.AvgRelError := Rep.AvgRelError+AbsReal(BarycentricCalc(P, XOriginal[I])-YOriginal[I])/AbsReal(YOriginal[I]);
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