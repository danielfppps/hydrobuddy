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
unit spline1d;
interface
uses Math, Sysutils, Ap, spline3, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, apserv;

type
(*************************************************************************
1-dimensional spline inteprolant
*************************************************************************)
Spline1DInterpolant = record
    Periodic : Boolean;
    N : AlglibInteger;
    K : AlglibInteger;
    X : TReal1DArray;
    C : TReal1DArray;
end;


(*************************************************************************
Spline fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************)
Spline1DFitReport = record
    TaskRCond : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    MaxError : Double;
end;



procedure Spline1DBuildLinear(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : Spline1DInterpolant);
procedure Spline1DBuildCubic(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     BoundLType : AlglibInteger;
     BoundL : Double;
     BoundRType : AlglibInteger;
     BoundR : Double;
     var C : Spline1DInterpolant);
procedure Spline1DBuildCatmullRom(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     BoundType : AlglibInteger;
     Tension : Double;
     var C : Spline1DInterpolant);
procedure Spline1DBuildHermite(X : TReal1DArray;
     Y : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var C : Spline1DInterpolant);
procedure Spline1DBuildAkima(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : Spline1DInterpolant);
procedure Spline1DFitCubicWC(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
procedure Spline1DFitHermiteWC(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
procedure Spline1DFitCubic(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
procedure Spline1DFitHermite(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
function Spline1DCalc(const C : Spline1DInterpolant; X : Double):Double;
procedure Spline1DDiff(const C : Spline1DInterpolant;
     X : Double;
     var S : Double;
     var DS : Double;
     var D2S : Double);
procedure Spline1DCopy(const C : Spline1DInterpolant;
     var CC : Spline1DInterpolant);
procedure Spline1DUnpack(const C : Spline1DInterpolant;
     var N : AlglibInteger;
     var Tbl : TReal2DArray);
procedure Spline1DLinTransX(var C : Spline1DInterpolant;
     A : Double;
     B : Double);
procedure Spline1DLinTransY(var C : Spline1DInterpolant;
     A : Double;
     B : Double);
function Spline1DIntegrate(const C : Spline1DInterpolant; X : Double):Double;

implementation

const
    Spline1DVNum = 11;

procedure Spline1DFitInternal(ST : AlglibInteger;
     X : TReal1DArray;
     Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     XC : TReal1DArray;
     YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);forward;
procedure HeapSortPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger);forward;
procedure HeapSortDPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     var D : TReal1DArray;
     N : AlglibInteger);forward;
procedure SolveTridiagonal(A : TReal1DArray;
     B : TReal1DArray;
     C : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var X : TReal1DArray);forward;
procedure SolveCyclicTridiagonal(const A : TReal1DArray;
     B : TReal1DArray;
     const C : TReal1DArray;
     const D : TReal1DArray;
     N : AlglibInteger;
     var X : TReal1DArray);forward;
function DiffThreePoint(T : Double;
     X0 : Double;
     F0 : Double;
     X1 : Double;
     F1 : Double;
     X2 : Double;
     F2 : Double):Double;forward;


(*************************************************************************
This subroutine builds linear spline interpolant

INPUT PARAMETERS:
    X   -   spline nodes, array[0..N-1]
    Y   -   function values, array[0..N-1]
    N   -   points count, N>=2
    
OUTPUT PARAMETERS:
    C   -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

  -- ALGLIB PROJECT --
     Copyright 24.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DBuildLinear(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : Spline1DInterpolant);
var
    I : AlglibInteger;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>1, 'Spline1DBuildLinear: N<2!');
    
    //
    // Sort points
    //
    HeapSortPoints(X, Y, N);
    
    //
    // Build
    //
    C.Periodic := False;
    C.N := N;
    C.K := 3;
    SetLength(C.X, N);
    SetLength(C.C, 4*(N-1));
    I:=0;
    while I<=N-1 do
    begin
        C.X[I] := X[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        C.C[4*I+0] := Y[I];
        C.C[4*I+1] := (Y[I+1]-Y[I])/(X[I+1]-X[I]);
        C.C[4*I+2] := 0;
        C.C[4*I+3] := 0;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine builds cubic spline interpolant.

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1].
    Y           -   function values, array[0..N-1].
    N           -   points count, N>=2
    BoundLType  -   boundary condition type for the left boundary
    BoundL      -   left boundary condition (first or second derivative,
                    depending on the BoundLType)
    BoundRType  -   boundary condition type for the right boundary
    BoundR      -   right boundary condition (first or second derivative,
                    depending on the BoundRType)

OUTPUT PARAMETERS:
    C           -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING BOUNDARY VALUES:

The BoundLType/BoundRType parameters can have the following values:
    * -1, which corresonds to the periodic (cyclic) boundary conditions.
          In this case:
          * both BoundLType and BoundRType must be equal to -1.
          * BoundL/BoundR are ignored
          * Y[last] is ignored (it is assumed to be equal to Y[first]).
    *  0, which  corresponds  to  the  parabolically   terminated  spline
          (BoundL and/or BoundR are ignored).
    *  1, which corresponds to the first derivative boundary condition
    *  2, which corresponds to the second derivative boundary condition

PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:

Problems with periodic boundary conditions have Y[first_point]=Y[last_point].
However, this subroutine doesn't require you to specify equal  values  for
the first and last points - it automatically forces them to be equal.

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DBuildCubic(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     BoundLType : AlglibInteger;
     BoundL : Double;
     BoundRType : AlglibInteger;
     BoundR : Double;
     var C : Spline1DInterpolant);
var
    A1 : TReal1DArray;
    A2 : TReal1DArray;
    A3 : TReal1DArray;
    B : TReal1DArray;
    D : TReal1DArray;
    DT : TReal1DArray;
    I : AlglibInteger;
    V : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>=2, 'Spline1DBuildCubic: N<2!');
    Assert((BoundLType=-1) or (BoundLType=0) or (BoundLType=1) or (BoundLType=2), 'Spline1DBuildCubic: incorrect BoundLType!');
    Assert((BoundRType=-1) or (BoundRType=0) or (BoundRType=1) or (BoundRType=2), 'Spline1DBuildCubic: incorrect BoundRType!');
    Assert((BoundRType=-1) and (BoundLType=-1) or (BoundRType<>-1) and (BoundLType<>-1), 'Spline1DBuildCubic: incorrect BoundLType/BoundRType!');
    
    //
    // Special cases:
    // * N=2, parabolic terminated boundary condition on both ends
    // * N=2, periodic boundary condition
    //
    if (N=2) and (BoundLType=0) and (BoundRType=0) then
    begin
        
        //
        // Change task type
        //
        BoundLType := 2;
        BoundL := 0;
        BoundRType := 2;
        BoundR := 0;
    end;
    if (N=2) and (BoundLType=-1) and (BoundRType=-1) then
    begin
        
        //
        // Change task type
        //
        BoundLType := 1;
        BoundL := 0;
        BoundRType := 1;
        BoundR := 0;
        Y[1] := Y[0];
    end;
    
    //
    // Periodic and non-periodic boundary conditions are
    // two separate classes
    //
    if (BoundRType=-1) and (BoundLType=-1) then
    begin
        
        //
        // Periodic boundary conditions
        //
        SetLength(A1, N-1);
        SetLength(A2, N-1);
        SetLength(A3, N-1);
        SetLength(B, N-1);
        
        //
        // Sort points.
        //
        HeapSortPoints(X, Y, N);
        Y[N-1] := Y[0];
        
        //
        // Boundary conditions at N-1 points
        // (one point less because last point is the same as first point).
        //
        A1[0] := X[1]-X[0];
        A2[0] := 2*(X[1]-X[0]+X[N-1]-X[N-2]);
        A3[0] := X[N-1]-X[N-2];
        B[0] := 3*(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2])*(X[1]-X[0])+3*(Y[1]-Y[0])/(X[1]-X[0])*(X[N-1]-X[N-2]);
        I:=1;
        while I<=N-2 do
        begin
            
            //
            // Altough last point is [N-2], we use X[N-1] and Y[N-1]
            // (because of periodicity)
            //
            A1[I] := X[I+1]-X[I];
            A2[I] := 2*(X[I+1]-X[I-1]);
            A3[I] := X[I]-X[I-1];
            B[I] := 3*(Y[I]-Y[I-1])/(X[I]-X[I-1])*(X[I+1]-X[I])+3*(Y[I+1]-Y[I])/(X[I+1]-X[I])*(X[I]-X[I-1]);
            Inc(I);
        end;
        
        //
        // Solve, add last point (with index N-1)
        //
        SolveCyclicTridiagonal(A1, A2, A3, B, N-1, DT);
        SetLength(D, N);
        APVMove(@D[0], 0, N-2, @DT[0], 0, N-2);
        D[N-1] := D[0];
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        Spline1DBuildHermite(X, Y, D, N, C);
        C.Periodic := True;
    end
    else
    begin
        
        //
        // Non-periodic boundary condition
        //
        SetLength(A1, N);
        SetLength(A2, N);
        SetLength(A3, N);
        SetLength(B, N);
        
        //
        // Sort points.
        //
        HeapSortPoints(X, Y, N);
        
        //
        // Left boundary conditions
        //
        if BoundLType=0 then
        begin
            A1[0] := 0;
            A2[0] := 1;
            A3[0] := 1;
            B[0] := 2*(Y[1]-Y[0])/(X[1]-X[0]);
        end;
        if BoundLType=1 then
        begin
            A1[0] := 0;
            A2[0] := 1;
            A3[0] := 0;
            B[0] := BoundL;
        end;
        if BoundLType=2 then
        begin
            A1[0] := 0;
            A2[0] := 2;
            A3[0] := 1;
            B[0] := 3*(Y[1]-Y[0])/(X[1]-X[0])-Double(0.5)*BoundL*(X[1]-X[0]);
        end;
        
        //
        // Central conditions
        //
        I:=1;
        while I<=N-2 do
        begin
            A1[I] := X[I+1]-X[I];
            A2[I] := 2*(X[I+1]-X[I-1]);
            A3[I] := X[I]-X[I-1];
            B[I] := 3*(Y[I]-Y[I-1])/(X[I]-X[I-1])*(X[I+1]-X[I])+3*(Y[I+1]-Y[I])/(X[I+1]-X[I])*(X[I]-X[I-1]);
            Inc(I);
        end;
        
        //
        // Right boundary conditions
        //
        if BoundRType=0 then
        begin
            A1[N-1] := 1;
            A2[N-1] := 1;
            A3[N-1] := 0;
            B[N-1] := 2*(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2]);
        end;
        if BoundRType=1 then
        begin
            A1[N-1] := 0;
            A2[N-1] := 1;
            A3[N-1] := 0;
            B[N-1] := BoundR;
        end;
        if BoundRType=2 then
        begin
            A1[N-1] := 1;
            A2[N-1] := 2;
            A3[N-1] := 0;
            B[N-1] := 3*(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2])+Double(0.5)*BoundR*(X[N-1]-X[N-2]);
        end;
        
        //
        // Solve
        //
        SolveTridiagonal(A1, A2, A3, B, N, D);
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        Spline1DBuildHermite(X, Y, D, N, C);
    end;
end;


(*************************************************************************
This subroutine builds Catmull-Rom spline interpolant.

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1].
    Y           -   function values, array[0..N-1].
    N           -   points count, N>=2
    BoundType   -   boundary condition type:
                    * -1 for periodic boundary condition
                    *  0 for parabolically terminated spline
    Tension     -   tension parameter:
                    * tension=0   corresponds to classic Catmull-Rom spline
                    * 0<tension<1 corresponds to more general form - cardinal spline

OUTPUT PARAMETERS:
    C           -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:

Problems with periodic boundary conditions have Y[first_point]=Y[last_point].
However, this subroutine doesn't require you to specify equal  values  for
the first and last points - it automatically forces them to be equal.

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DBuildCatmullRom(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     BoundType : AlglibInteger;
     Tension : Double;
     var C : Spline1DInterpolant);
var
    A1 : TReal1DArray;
    A2 : TReal1DArray;
    A3 : TReal1DArray;
    B : TReal1DArray;
    D : TReal1DArray;
    DT : TReal1DArray;
    I : AlglibInteger;
    V : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>=2, 'Spline1DBuildCatmullRom: N<2!');
    Assert((BoundType=-1) or (BoundType=0), 'Spline1DBuildCatmullRom: incorrect BoundType!');
    
    //
    // Special cases:
    // * N=2, parabolic terminated boundary condition on both ends
    // * N=2, periodic boundary condition
    //
    if (N=2) and (BoundType=0) then
    begin
        
        //
        // Just linear spline
        //
        Spline1DBuildLinear(X, Y, N, C);
        Exit;
    end;
    if (N=2) and (BoundType=-1) then
    begin
        
        //
        // Same as cubic spline with periodic conditions
        //
        Spline1DBuildCubic(X, Y, N, -1, Double(0.0), -1, Double(0.0), C);
        Exit;
    end;
    
    //
    // Periodic or non-periodic boundary conditions
    //
    if BoundType=-1 then
    begin
        
        //
        // Sort points.
        //
        HeapSortPoints(X, Y, N);
        Y[N-1] := Y[0];
        
        //
        // Periodic boundary conditions
        //
        SetLength(D, N);
        D[0] := (Y[1]-Y[N-2])/(2*(X[1]-X[0]+X[N-1]-X[N-2]));
        I:=1;
        while I<=N-2 do
        begin
            D[I] := (1-Tension)*(Y[I+1]-Y[I-1])/(X[I+1]-X[I-1]);
            Inc(I);
        end;
        D[N-1] := D[0];
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        Spline1DBuildHermite(X, Y, D, N, C);
        C.Periodic := True;
    end
    else
    begin
        
        //
        // Sort points.
        //
        HeapSortPoints(X, Y, N);
        
        //
        // Non-periodic boundary conditions
        //
        SetLength(D, N);
        I:=1;
        while I<=N-2 do
        begin
            D[I] := (1-Tension)*(Y[I+1]-Y[I-1])/(X[I+1]-X[I-1]);
            Inc(I);
        end;
        D[0] := 2*(Y[1]-Y[0])/(X[1]-X[0])-D[1];
        D[N-1] := 2*(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2])-D[N-2];
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        Spline1DBuildHermite(X, Y, D, N, C);
    end;
end;


(*************************************************************************
This subroutine builds Hermite spline interpolant.

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1]
    Y           -   function values, array[0..N-1]
    D           -   derivatives, array[0..N-1]
    N           -   points count, N>=2

OUTPUT PARAMETERS:
    C           -   spline interpolant.


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DBuildHermite(X : TReal1DArray;
     Y : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var C : Spline1DInterpolant);
var
    I : AlglibInteger;
    Delta : Double;
    Delta2 : Double;
    Delta3 : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    D := DynamicArrayCopy(D);
    Assert(N>=2, 'BuildHermiteSpline: N<2!');
    
    //
    // Sort points
    //
    HeapSortDPoints(X, Y, D, N);
    
    //
    // Build
    //
    SetLength(C.X, N);
    SetLength(C.C, 4*(N-1));
    C.Periodic := False;
    C.K := 3;
    C.N := N;
    I:=0;
    while I<=N-1 do
    begin
        C.X[I] := X[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        Delta := X[I+1]-X[I];
        Delta2 := AP_Sqr(Delta);
        Delta3 := Delta*Delta2;
        C.C[4*I+0] := Y[I];
        C.C[4*I+1] := D[I];
        C.C[4*I+2] := (3*(Y[I+1]-Y[I])-2*D[I]*Delta-D[I+1]*Delta)/Delta2;
        C.C[4*I+3] := (2*(Y[I]-Y[I+1])+D[I]*Delta+D[I+1]*Delta)/Delta3;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine builds Akima spline interpolant

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1]
    Y           -   function values, array[0..N-1]
    N           -   points count, N>=5

OUTPUT PARAMETERS:
    C           -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

  -- ALGLIB PROJECT --
     Copyright 24.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DBuildAkima(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : Spline1DInterpolant);
var
    I : AlglibInteger;
    D : TReal1DArray;
    W : TReal1DArray;
    Diff : TReal1DArray;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>=5, 'BuildAkimaSpline: N<5!');
    
    //
    // Sort points
    //
    HeapSortPoints(X, Y, N);
    
    //
    // Prepare W (weights), Diff (divided differences)
    //
    SetLength(W, N-1);
    SetLength(Diff, N-1);
    I:=0;
    while I<=N-2 do
    begin
        Diff[I] := (Y[I+1]-Y[I])/(X[I+1]-X[I]);
        Inc(I);
    end;
    I:=1;
    while I<=N-2 do
    begin
        W[I] := AbsReal(Diff[I]-Diff[I-1]);
        Inc(I);
    end;
    
    //
    // Prepare Hermite interpolation scheme
    //
    SetLength(D, N);
    I:=2;
    while I<=N-3 do
    begin
        if AP_FP_Neq(AbsReal(W[I-1])+AbsReal(W[I+1]),0) then
        begin
            D[I] := (W[I+1]*Diff[I-1]+W[I-1]*Diff[I])/(W[I+1]+W[I-1]);
        end
        else
        begin
            D[I] := ((X[I+1]-X[I])*Diff[I-1]+(X[I]-X[I-1])*Diff[I])/(X[I+1]-X[I-1]);
        end;
        Inc(I);
    end;
    D[0] := DiffThreePoint(X[0], X[0], Y[0], X[1], Y[1], X[2], Y[2]);
    D[1] := DiffThreePoint(X[1], X[0], Y[0], X[1], Y[1], X[2], Y[2]);
    D[N-2] := DiffThreePoint(X[N-2], X[N-3], Y[N-3], X[N-2], Y[N-2], X[N-1], Y[N-1]);
    D[N-1] := DiffThreePoint(X[N-1], X[N-3], Y[N-3], X[N-2], Y[N-2], X[N-1], Y[N-1]);
    
    //
    // Build Akima spline using Hermite interpolation scheme
    //
    Spline1DBuildHermite(X, Y, D, N, C);
end;


(*************************************************************************
Weighted fitting by cubic  spline,  with constraints on function values or
derivatives.

Equidistant grid with M-2 nodes on [min(x,xc),max(x,xc)] is  used to build
basis functions. Basis functions are cubic splines with continuous  second
derivatives  and  non-fixed first  derivatives  at  interval  ends.  Small
regularizing term is used  when  solving  constrained  tasks  (to  improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitHermiteWC()  -   fitting by Hermite splines (more flexible,
                                less smooth)
    Spline1DFitCubic()      -   "lightweight" fitting  by  cubic  splines,
                                without invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions ( = number_of_nodes+2), M>=4.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints  on  the  function  values  AND/OR  its
  derivatives at the interval boundaries.
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.


  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DFitCubicWC(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
begin
    Spline1DFitInternal(0, X, Y, W, N, XC, YC, DC, K, M, Info, S, Rep);
end;


(*************************************************************************
Weighted  fitting  by Hermite spline,  with constraints on function values
or first derivatives.

Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is  used  to  build
basis functions. Basis functions are Hermite splines.  Small  regularizing
term is used when solving constrained tasks (to improve stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitCubicWC()    -   fitting by Cubic splines (less flexible,
                                more smooth)
    Spline1DFitHermite()    -   "lightweight" Hermite fitting, without
                                invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions (= 2 * number of nodes),
            M>=4,
            M IS EVEN!

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -2 means odd M was passed (which is not supported)
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

IMPORTANT:
    this subroitine supports only even M's


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the several special cases, however, we can guarantee consistency.
* one of this cases is  M>=4  and   constraints  on   the  function  value
  (AND/OR its derivative) at the interval boundaries.
* another special case is M>=4  and  ONE  constraint on the function value
  (OR, BUT NOT AND, derivative) anywhere in [min(x),max(x)]

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DFitHermiteWC(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
begin
    Spline1DFitInternal(1, X, Y, W, N, XC, YC, DC, K, M, Info, S, Rep);
end;


(*************************************************************************
Least squares fitting by cubic spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitCubicWC().  See  Spline1DFitCubicWC() for more information
about subroutine parameters (we don't duplicate it here because of length)

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DFitCubic(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
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
    Spline1DFitCubicWC(X, Y, W, N, XC, YC, DC, 0, M, Info, S, Rep);
end;


(*************************************************************************
Least squares fitting by Hermite spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitHermiteWC().  See Spline1DFitHermiteWC()  description  for
more information about subroutine parameters (we don't duplicate  it  here
because of length).

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DFitHermite(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
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
    Spline1DFitHermiteWC(X, Y, W, N, XC, YC, DC, 0, M, Info, S, Rep);
end;


(*************************************************************************
This subroutine calculates the value of the spline at the given point X.

INPUT PARAMETERS:
    C   -   spline interpolant
    X   -   point

Result:
    S(x)

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************)
function Spline1DCalc(const C : Spline1DInterpolant; X : Double):Double;
var
    L : AlglibInteger;
    R : AlglibInteger;
    M : AlglibInteger;
    T : Double;
begin
    Assert(C.K=3, 'Spline1DCalc: internal error');
    
    //
    // correct if periodic
    //
    if C.Periodic then
    begin
        APPeriodicMap(X, C.X[0], C.X[C.N-1], T);
    end;
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    L := 0;
    R := C.N-2+1;
    while L<>R-1 do
    begin
        M := (L+R) div 2;
        if AP_FP_Greater_Eq(C.X[M],X) then
        begin
            R := M;
        end
        else
        begin
            L := M;
        end;
    end;
    
    //
    // Interpolation
    //
    X := X-C.X[L];
    M := 4*L;
    Result := C.C[M]+X*(C.C[M+1]+X*(C.C[M+2]+X*C.C[M+3]));
end;


(*************************************************************************
This subroutine differentiates the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.
    X   -   point

Result:
    S   -   S(x)
    DS  -   S'(x)
    D2S -   S''(x)

  -- ALGLIB PROJECT --
     Copyright 24.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DDiff(const C : Spline1DInterpolant;
     X : Double;
     var S : Double;
     var DS : Double;
     var D2S : Double);
var
    L : AlglibInteger;
    R : AlglibInteger;
    M : AlglibInteger;
    T : Double;
begin
    Assert(C.K=3, 'Spline1DCalc: internal error');
    
    //
    // correct if periodic
    //
    if C.Periodic then
    begin
        APPeriodicMap(X, C.X[0], C.X[C.N-1], T);
    end;
    
    //
    // Binary search
    //
    L := 0;
    R := C.N-2+1;
    while L<>R-1 do
    begin
        M := (L+R) div 2;
        if AP_FP_Greater_Eq(C.X[M],X) then
        begin
            R := M;
        end
        else
        begin
            L := M;
        end;
    end;
    
    //
    // Differentiation
    //
    X := X-C.X[L];
    M := 4*L;
    S := C.C[M]+X*(C.C[M+1]+X*(C.C[M+2]+X*C.C[M+3]));
    DS := C.C[M+1]+2*X*C.C[M+2]+3*AP_Sqr(X)*C.C[M+3];
    D2S := 2*C.C[M+2]+6*X*C.C[M+3];
end;


(*************************************************************************
This subroutine makes the copy of the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.

Result:
    CC  -   spline copy

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DCopy(const C : Spline1DInterpolant;
     var CC : Spline1DInterpolant);
begin
    CC.Periodic := C.Periodic;
    CC.N := C.N;
    CC.K := C.K;
    SetLength(CC.X, CC.N);
    APVMove(@CC.X[0], 0, CC.N-1, @C.X[0], 0, CC.N-1);
    SetLength(CC.C, (CC.K+1)*(CC.N-1));
    APVMove(@CC.C[0], 0, (CC.K+1)*(CC.N-1)-1, @C.C[0], 0, (CC.K+1)*(CC.N-1)-1);
end;


(*************************************************************************
This subroutine unpacks the spline into the coefficients table.

INPUT PARAMETERS:
    C   -   spline interpolant.
    X   -   point

Result:
    Tbl -   coefficients table, unpacked format, array[0..N-2, 0..5].
            For I = 0...N-2:
                Tbl[I,0] = X[i]
                Tbl[I,1] = X[i+1]
                Tbl[I,2] = C0
                Tbl[I,3] = C1
                Tbl[I,4] = C2
                Tbl[I,5] = C3
            On [x[i], x[i+1]] spline is equals to:
                S(x) = C0 + C1*t + C2*t^2 + C3*t^3
                t = x-x[i]

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DUnpack(const C : Spline1DInterpolant;
     var N : AlglibInteger;
     var Tbl : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(Tbl, C.N-2+1, 2+C.K+1);
    N := C.N;
    
    //
    // Fill
    //
    I:=0;
    while I<=N-2 do
    begin
        Tbl[I,0] := C.X[I];
        Tbl[I,1] := C.X[I+1];
        J:=0;
        while J<=C.K do
        begin
            Tbl[I,2+J] := C.C[(C.K+1)*I+J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine performs linear transformation of the spline argument.

INPUT PARAMETERS:
    C   -   spline interpolant.
    A, B-   transformation coefficients: x = A*t + B
Result:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DLinTransX(var C : Spline1DInterpolant;
     A : Double;
     B : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
    V : Double;
    DV : Double;
    D2V : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    D : TReal1DArray;
begin
    N := C.N;
    
    //
    // Special case: A=0
    //
    if AP_FP_Eq(A,0) then
    begin
        V := Spline1DCalc(C, B);
        I:=0;
        while I<=N-2 do
        begin
            C.C[(C.K+1)*I] := V;
            J:=1;
            while J<=C.K do
            begin
                C.C[(C.K+1)*I+J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // General case: A<>0.
    // Unpack, X, Y, dY/dX.
    // Scale and pack again.
    //
    Assert(C.K=3, 'Spline1DLinTransX: internal error');
    SetLength(X, N-1+1);
    SetLength(Y, N-1+1);
    SetLength(D, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := C.X[I];
        Spline1DDiff(C, X[I], V, DV, D2V);
        X[I] := (X[I]-B)/A;
        Y[I] := V;
        D[I] := A*DV;
        Inc(I);
    end;
    Spline1DBuildHermite(X, Y, D, N, C);
end;


(*************************************************************************
This subroutine performs linear transformation of the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.
    A, B-   transformation coefficients: S2(x) = A*S(x) + B
Result:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DLinTransY(var C : Spline1DInterpolant;
     A : Double;
     B : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
begin
    N := C.N;
    I:=0;
    while I<=N-2 do
    begin
        C.C[(C.K+1)*I] := A*C.C[(C.K+1)*I]+B;
        J:=1;
        while J<=C.K do
        begin
            C.C[(C.K+1)*I+J] := A*C.C[(C.K+1)*I+J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine integrates the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.
    X   -   right bound of the integration interval [a, x],
            here 'a' denotes min(x[])
Result:
    integral(S(t)dt,a,x)

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************)
function Spline1DIntegrate(const C : Spline1DInterpolant; X : Double):Double;
var
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    M : AlglibInteger;
    W : Double;
    V : Double;
    T : Double;
    IntAB : Double;
    AdditionalTerm : Double;
begin
    N := C.N;
    
    //
    // Periodic splines require special treatment. We make
    // following transformation:
    //
    //     integral(S(t)dt,A,X) = integral(S(t)dt,A,Z)+AdditionalTerm
    //
    // here X may lie outside of [A,B], Z lies strictly in [A,B],
    // AdditionalTerm is equals to integral(S(t)dt,A,B) times some
    // integer number (may be zero).
    //
    if C.Periodic and (AP_FP_Less(X,C.X[0]) or AP_FP_Greater(X,C.X[C.N-1])) then
    begin
        
        //
        // compute integral(S(x)dx,A,B)
        //
        IntAB := 0;
        I:=0;
        while I<=C.N-2 do
        begin
            W := C.X[I+1]-C.X[I];
            M := (C.K+1)*I;
            IntAB := IntAB+C.C[M]*W;
            V := W;
            J:=1;
            while J<=C.K do
            begin
                V := V*W;
                IntAB := IntAB+C.C[M+J]*V/(J+1);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // map X into [A,B]
        //
        APPeriodicMap(X, C.X[0], C.X[C.N-1], T);
        AdditionalTerm := T*IntAB;
    end
    else
    begin
        AdditionalTerm := 0;
    end;
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    L := 0;
    R := N-2+1;
    while L<>R-1 do
    begin
        M := (L+R) div 2;
        if AP_FP_Greater_Eq(C.X[M],X) then
        begin
            R := M;
        end
        else
        begin
            L := M;
        end;
    end;
    
    //
    // Integration
    //
    Result := 0;
    I:=0;
    while I<=L-1 do
    begin
        W := C.X[I+1]-C.X[I];
        M := (C.K+1)*I;
        Result := Result+C.C[M]*W;
        V := W;
        J:=1;
        while J<=C.K do
        begin
            V := V*W;
            Result := Result+C.C[M+J]*V/(J+1);
            Inc(J);
        end;
        Inc(I);
    end;
    W := X-C.X[L];
    M := (C.K+1)*L;
    V := W;
    Result := Result+C.C[M]*W;
    J:=1;
    while J<=C.K do
    begin
        V := V*W;
        Result := Result+C.C[M+J]*V/(J+1);
        Inc(J);
    end;
    Result := Result+AdditionalTerm;
end;


(*************************************************************************
Internal spline fitting subroutine

  -- ALGLIB PROJECT --
     Copyright 08.09.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline1DFitInternal(ST : AlglibInteger;
     X : TReal1DArray;
     Y : TReal1DArray;
     const W : TReal1DArray;
     N : AlglibInteger;
     XC : TReal1DArray;
     YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var S : Spline1DInterpolant;
     var Rep : Spline1DFitReport);
var
    FMatrix : TReal2DArray;
    CMatrix : TReal2DArray;
    Y2 : TReal1DArray;
    W2 : TReal1DArray;
    SX : TReal1DArray;
    SY : TReal1DArray;
    SD : TReal1DArray;
    Tmp : TReal1DArray;
    XOriginal : TReal1DArray;
    YOriginal : TReal1DArray;
    LRep : LSFitReport;
    V0 : Double;
    V1 : Double;
    V2 : Double;
    MX : Double;
    S2 : Spline1DInterpolant;
    I : AlglibInteger;
    J : AlglibInteger;
    RelCnt : AlglibInteger;
    XA : Double;
    XB : Double;
    SA : Double;
    SB : Double;
    BL : Double;
    BR : Double;
    Decay : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    XC := DynamicArrayCopy(XC);
    YC := DynamicArrayCopy(YC);
    Assert((ST=0) or (ST=1), 'Spline1DFit: internal error!');
    if (ST=0) and (M<4) then
    begin
        Info := -1;
        Exit;
    end;
    if (ST=1) and (M<4) then
    begin
        Info := -1;
        Exit;
    end;
    if (N<1) or (K<0) or (K>=M) then
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
    if (ST=1) and (M mod 2<>0) then
    begin
        
        //
        // Hermite fitter must have even number of basis functions
        //
        Info := -2;
        Exit;
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
    // * SX     -   grid for basis functions
    // * SY     -   values of basis functions at grid points
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
    if ST=0 then
    begin
        
        //
        // allocate space for cubic spline
        //
        SetLength(SX, M-2);
        SetLength(SY, M-2);
        J:=0;
        while J<=M-2-1 do
        begin
            SX[J] := AP_Double(2*J)/(M-2-1)-1;
            Inc(J);
        end;
    end;
    if ST=1 then
    begin
        
        //
        // allocate space for Hermite spline
        //
        SetLength(SX, M div 2);
        SetLength(SY, M div 2);
        SetLength(SD, M div 2);
        J:=0;
        while J<=M div 2-1 do
        begin
            SX[J] := AP_Double(2*J)/(M div 2-1)-1;
            Inc(J);
        end;
    end;
    
    //
    // Prepare design and constraints matrices:
    // * fill constraints matrix
    // * fill first N rows of design matrix with values
    // * fill next M rows of design matrix with regularizing term
    // * append M zeros to Y
    // * append M elements, mean(abs(W)) each, to W
    //
    J:=0;
    while J<=M-1 do
    begin
        
        //
        // prepare Jth basis function
        //
        if ST=0 then
        begin
            
            //
            // cubic spline basis
            //
            I:=0;
            while I<=M-2-1 do
            begin
                SY[I] := 0;
                Inc(I);
            end;
            BL := 0;
            BR := 0;
            if J<M-2 then
            begin
                SY[J] := 1;
            end;
            if J=M-2 then
            begin
                BL := 1;
            end;
            if J=M-1 then
            begin
                BR := 1;
            end;
            Spline1DBuildCubic(SX, SY, M-2, 1, BL, 1, BR, S2);
        end;
        if ST=1 then
        begin
            
            //
            // Hermite basis
            //
            I:=0;
            while I<=M div 2-1 do
            begin
                SY[I] := 0;
                SD[I] := 0;
                Inc(I);
            end;
            if J mod 2=0 then
            begin
                SY[J div 2] := 1;
            end
            else
            begin
                SD[J div 2] := 1;
            end;
            Spline1DBuildHermite(SX, SY, SD, M div 2, S2);
        end;
        
        //
        // values at X[], XC[]
        //
        I:=0;
        while I<=N-1 do
        begin
            FMatrix[I,J] := Spline1DCalc(S2, X[I]);
            Inc(I);
        end;
        I:=0;
        while I<=K-1 do
        begin
            Assert((DC[I]>=0) and (DC[I]<=2), 'Spline1DFit: internal error!');
            Spline1DDiff(S2, XC[I], V0, V1, V2);
            if DC[I]=0 then
            begin
                CMatrix[I,J] := V0;
            end;
            if DC[I]=1 then
            begin
                CMatrix[I,J] := V1;
            end;
            if DC[I]=2 then
            begin
                CMatrix[I,J] := V2;
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
    SetLength(Y2, N+M);
    SetLength(W2, N+M);
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
    // Generate spline and scale it
    //
    if ST=0 then
    begin
        
        //
        // cubic spline basis
        //
        APVMove(@SY[0], 0, M-2-1, @Tmp[0], 0, M-2-1);
        Spline1DBuildCubic(SX, SY, M-2, 1, Tmp[M-2], 1, Tmp[M-1], S);
    end;
    if ST=1 then
    begin
        
        //
        // Hermite basis
        //
        I:=0;
        while I<=M div 2-1 do
        begin
            SY[I] := Tmp[2*I];
            SD[I] := Tmp[2*I+1];
            Inc(I);
        end;
        Spline1DBuildHermite(SX, SY, SD, M div 2, S);
    end;
    Spline1DLinTransX(S, 2/(XB-XA), -(XA+XB)/(XB-XA));
    Spline1DLinTransY(S, SB-SA, SA);
    
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
            Rep.AvgRelError := Rep.AvgRelError+AbsReal(Spline1DCalc(S, XOriginal[I])-YOriginal[I])/AbsReal(YOriginal[I]);
            RelCnt := RelCnt+1;
        end;
        Inc(I);
    end;
    if RelCnt<>0 then
    begin
        Rep.AvgRelError := Rep.AvgRelError/RelCnt;
    end;
end;


(*************************************************************************
Internal subroutine. Heap sort.
*************************************************************************)
procedure HeapSortPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    IsAscending : Boolean;
    IsDescending : Boolean;
begin
    
    //
    // Test for already sorted set
    //
    IsAscending := True;
    IsDescending := True;
    I:=1;
    while I<=N-1 do
    begin
        IsAscending := IsAscending and AP_FP_Greater(X[I],X[I-1]);
        IsDescending := IsDescending and AP_FP_Less(X[I],X[I-1]);
        Inc(I);
    end;
    if IsAscending then
    begin
        Exit;
    end;
    if IsDescending then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J := N-1-I;
            if J<=I then
            begin
                Break;
            end;
            Tmp := X[I];
            X[I] := X[J];
            X[J] := Tmp;
            Tmp := Y[I];
            Y[I] := Y[J];
            Y[J] := Tmp;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        Exit;
    end;
    
    //
    // General case
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(X[k-1],X[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := X[k-1];
                X[k-1] := X[t-1];
                X[t-1] := Tmp;
                Tmp := Y[k-1];
                Y[k-1] := Y[t-1];
                Y[t-1] := Tmp;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := X[i];
        X[i] := X[0];
        X[0] := Tmp;
        Tmp := Y[i];
        Y[i] := Y[0];
        Y[0] := Tmp;
        t := 1;
        while t<>0 do
        begin
            k := 2*t;
            if k>i then
            begin
                t := 0;
            end
            else
            begin
                if k<i then
                begin
                    if AP_FP_Greater(X[k],X[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(X[t-1],X[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := X[k-1];
                    X[k-1] := X[t-1];
                    X[t-1] := Tmp;
                    Tmp := Y[k-1];
                    Y[k-1] := Y[t-1];
                    Y[t-1] := Tmp;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


(*************************************************************************
Internal subroutine. Heap sort.
*************************************************************************)
procedure HeapSortDPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     var D : TReal1DArray;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    IsAscending : Boolean;
    IsDescending : Boolean;
begin
    
    //
    // Test for already sorted set
    //
    IsAscending := True;
    IsDescending := True;
    I:=1;
    while I<=N-1 do
    begin
        IsAscending := IsAscending and AP_FP_Greater(X[I],X[I-1]);
        IsDescending := IsDescending and AP_FP_Less(X[I],X[I-1]);
        Inc(I);
    end;
    if IsAscending then
    begin
        Exit;
    end;
    if IsDescending then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J := N-1-I;
            if J<=I then
            begin
                Break;
            end;
            Tmp := X[I];
            X[I] := X[J];
            X[J] := Tmp;
            Tmp := Y[I];
            Y[I] := Y[J];
            Y[J] := Tmp;
            Tmp := D[I];
            D[I] := D[J];
            D[J] := Tmp;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        Exit;
    end;
    
    //
    // General case
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(X[k-1],X[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := X[k-1];
                X[k-1] := X[t-1];
                X[t-1] := Tmp;
                Tmp := Y[k-1];
                Y[k-1] := Y[t-1];
                Y[t-1] := Tmp;
                Tmp := D[k-1];
                D[k-1] := D[t-1];
                D[t-1] := Tmp;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := X[i];
        X[i] := X[0];
        X[0] := Tmp;
        Tmp := Y[i];
        Y[i] := Y[0];
        Y[0] := Tmp;
        Tmp := D[i];
        D[i] := D[0];
        D[0] := Tmp;
        t := 1;
        while t<>0 do
        begin
            k := 2*t;
            if k>i then
            begin
                t := 0;
            end
            else
            begin
                if k<i then
                begin
                    if AP_FP_Greater(X[k],X[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(X[t-1],X[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := X[k-1];
                    X[k-1] := X[t-1];
                    X[t-1] := Tmp;
                    Tmp := Y[k-1];
                    Y[k-1] := Y[t-1];
                    Y[t-1] := Tmp;
                    Tmp := D[k-1];
                    D[k-1] := D[t-1];
                    D[t-1] := Tmp;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


(*************************************************************************
Internal subroutine. Tridiagonal solver. Solves

( B[0] C[0]                      )
( A[1] B[1] C[1]                 )
(      A[2] B[2] C[2]            )
(            ..........          ) * X = D
(            ..........          )
(           A[N-2] B[N-2] C[N-2] )
(                  A[N-1] B[N-1] )

*************************************************************************)
procedure SolveTridiagonal(A : TReal1DArray;
     B : TReal1DArray;
     C : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var X : TReal1DArray);
var
    K : AlglibInteger;
    T : Double;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    C := DynamicArrayCopy(C);
    D := DynamicArrayCopy(D);
    SetLength(X, N-1+1);
    A[0] := 0;
    C[N-1] := 0;
    K:=1;
    while K<=N-1 do
    begin
        T := A[K]/B[K-1];
        B[K] := B[K]-T*C[K-1];
        D[K] := D[K]-T*D[K-1];
        Inc(K);
    end;
    X[N-1] := D[N-1]/B[N-1];
    K:=N-2;
    while K>=0 do
    begin
        X[K] := (D[K]-C[K]*X[K+1])/B[K];
        Dec(K);
    end;
end;


(*************************************************************************
Internal subroutine. Cyclic tridiagonal solver. Solves

( B[0] C[0]                 A[0] )
( A[1] B[1] C[1]                 )
(      A[2] B[2] C[2]            )
(            ..........          ) * X = D
(            ..........          )
(           A[N-2] B[N-2] C[N-2] )
( C[N-1]           A[N-1] B[N-1] )
*************************************************************************)
procedure SolveCyclicTridiagonal(const A : TReal1DArray;
     B : TReal1DArray;
     const C : TReal1DArray;
     const D : TReal1DArray;
     N : AlglibInteger;
     var X : TReal1DArray);
var
    K : AlglibInteger;
    T : Double;
    Alpha : Double;
    Beta : Double;
    Gamma : Double;
    Y : TReal1DArray;
    Z : TReal1DArray;
    U : TReal1DArray;
begin
    B := DynamicArrayCopy(B);
    Beta := A[0];
    Alpha := C[N-1];
    Gamma := -B[0];
    B[0] := 2*B[0];
    B[N-1] := B[N-1]-Alpha*Beta/Gamma;
    SetLength(U, N);
    K:=0;
    while K<=N-1 do
    begin
        U[K] := 0;
        Inc(K);
    end;
    U[0] := Gamma;
    U[N-1] := Alpha;
    SolveTridiagonal(A, B, C, D, N, Y);
    SolveTridiagonal(A, B, C, U, N, Z);
    SetLength(X, N);
    K:=0;
    while K<=N-1 do
    begin
        X[K] := Y[K]-(Y[0]+Beta/Gamma*Y[N-1])/(1+Z[0]+Beta/Gamma*Z[N-1])*Z[K];
        Inc(K);
    end;
end;


(*************************************************************************
Internal subroutine. Three-point differentiation
*************************************************************************)
function DiffThreePoint(T : Double;
     X0 : Double;
     F0 : Double;
     X1 : Double;
     F1 : Double;
     X2 : Double;
     F2 : Double):Double;
var
    A : Double;
    B : Double;
begin
    T := T-X0;
    X1 := X1-X0;
    X2 := X2-X0;
    A := (F2-F0-X2/X1*(F1-F0))/(AP_Sqr(X2)-X1*X2);
    B := (F1-F0-A*AP_Sqr(X1))/X1;
    Result := 2*A*T+B;
end;


end.