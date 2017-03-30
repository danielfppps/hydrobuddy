{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2006-2010, Sergey Bochkanov (ALGLIB project).

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
unit pspline;
interface
uses Math, Sysutils, Ap, spline3, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, apserv, spline1d, tsort, hsschur, evd, gammafunc, gq, gkq, autogk;

type
(*************************************************************************
Parametric spline inteprolant: 2-dimensional curve.

You should not try to access its members directly - use PSpline2XXXXXXXX()
functions instead.
*************************************************************************)
PSpline2Interpolant = record
    N : AlglibInteger;
    Periodic : Boolean;
    P : TReal1DArray;
    X : Spline1DInterpolant;
    Y : Spline1DInterpolant;
end;


(*************************************************************************
Parametric spline inteprolant: 3-dimensional curve.

You should not try to access its members directly - use PSpline3XXXXXXXX()
functions instead.
*************************************************************************)
PSpline3Interpolant = record
    N : AlglibInteger;
    Periodic : Boolean;
    P : TReal1DArray;
    X : Spline1DInterpolant;
    Y : Spline1DInterpolant;
    Z : Spline1DInterpolant;
end;



procedure PSpline2Build(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline2Interpolant);
procedure PSpline3Build(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline3Interpolant);
procedure PSpline2BuildPeriodic(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline2Interpolant);
procedure PSpline3BuildPeriodic(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline3Interpolant);
procedure PSpline2ParameterValues(const P : PSpline2Interpolant;
     var N : AlglibInteger;
     var T : TReal1DArray);
procedure PSpline3ParameterValues(const P : PSpline3Interpolant;
     var N : AlglibInteger;
     var T : TReal1DArray);
procedure PSpline2Calc(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var Y : Double);
procedure PSpline3Calc(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var Y : Double;
     var Z : Double);
procedure PSpline2Tangent(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var Y : Double);
procedure PSpline3Tangent(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var Y : Double;
     var Z : Double);
procedure PSpline2Diff(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var Y : Double;
     var DY : Double);
procedure PSpline3Diff(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var Y : Double;
     var DY : Double;
     var Z : Double;
     var DZ : Double);
procedure PSpline2Diff2(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var D2X : Double;
     var Y : Double;
     var DY : Double;
     var D2Y : Double);
procedure PSpline3Diff2(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var D2X : Double;
     var Y : Double;
     var DY : Double;
     var D2Y : Double;
     var Z : Double;
     var DZ : Double;
     var D2Z : Double);
function PSpline2ArcLength(const P : PSpline2Interpolant;
     A : Double;
     B : Double):Double;
function PSpline3ArcLength(const P : PSpline3Interpolant;
     A : Double;
     B : Double):Double;

implementation

procedure PSpline2Par(const XY : TReal2DArray;
     N : AlglibInteger;
     PT : AlglibInteger;
     var P : TReal1DArray);forward;
procedure PSpline3Par(const XY : TReal2DArray;
     N : AlglibInteger;
     PT : AlglibInteger;
     var P : TReal1DArray);forward;


(*************************************************************************
This function  builds  non-periodic 2-dimensional parametric spline  which
starts at (X[0],Y[0]) and ends at (X[N-1],Y[N-1]).

INPUT PARAMETERS:
    XY  -   points, array[0..N-1,0..1].
            XY[I,0:1] corresponds to the Ith point.
            Order of points is important!
    N   -   points count, N>=5 for Akima splines, N>=2 for other types  of
            splines.
    ST  -   spline type:
            * 0     Akima spline
            * 1     parabolically terminated Catmull-Rom spline (Tension=0)
            * 2     parabolically terminated cubic spline
    PT  -   parameterization type:
            * 0     uniform
            * 1     chord length
            * 2     centripetal

OUTPUT PARAMETERS:
    P   -   parametric spline interpolant


NOTES:
* this function  assumes  that  there all consequent points  are distinct.
  I.e. (x0,y0)<>(x1,y1),  (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so on.
  However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
  =(x2,y2).

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2Build(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline2Interpolant);
var
    Tmp : TReal1DArray;
    V : Double;
    I : AlglibInteger;
    i_ : AlglibInteger;
begin
    XY := DynamicArrayCopy(XY);
    Assert((ST>=0) and (ST<=2), 'PSpline2Build: incorrect spline type!');
    Assert((PT>=0) and (PT<=2), 'PSpline2Build: incorrect parameterization type!');
    if ST=0 then
    begin
        Assert(N>=5, 'PSpline2Build: N<5 (minimum value for Akima splines)!');
    end
    else
    begin
        Assert(N>=2, 'PSpline2Build: N<2!');
    end;
    
    //
    // Prepare
    //
    P.N := N;
    P.Periodic := False;
    SetLength(Tmp, N);
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    PSpline2Par(XY, N, PT, P.P);
    Assert(APSERVAreDistinct(P.P, N), 'PSpline2Build: consequent points are too close!');
    
    //
    // Build splines
    //
    if ST=0 then
    begin
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,0];
        end;
        Spline1DBuildAkima(P.P, Tmp, N, P.X);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,1];
        end;
        Spline1DBuildAkima(P.P, Tmp, N, P.Y);
    end;
    if ST=1 then
    begin
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,0];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N, 0, Double(0.0), P.X);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,1];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N, 0, Double(0.0), P.Y);
    end;
    if ST=2 then
    begin
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,0];
        end;
        Spline1DBuildCubic(P.P, Tmp, N, 0, Double(0.0), 0, Double(0.0), P.X);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,1];
        end;
        Spline1DBuildCubic(P.P, Tmp, N, 0, Double(0.0), 0, Double(0.0), P.Y);
    end;
end;


(*************************************************************************
This function  builds  non-periodic 3-dimensional parametric spline  which
starts at (X[0],Y[0],Z[0]) and ends at (X[N-1],Y[N-1],Z[N-1]).

Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
description here.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3Build(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline3Interpolant);
var
    Tmp : TReal1DArray;
    V : Double;
    I : AlglibInteger;
    i_ : AlglibInteger;
begin
    XY := DynamicArrayCopy(XY);
    Assert((ST>=0) and (ST<=2), 'PSpline3Build: incorrect spline type!');
    Assert((PT>=0) and (PT<=2), 'PSpline3Build: incorrect parameterization type!');
    if ST=0 then
    begin
        Assert(N>=5, 'PSpline3Build: N<5 (minimum value for Akima splines)!');
    end
    else
    begin
        Assert(N>=2, 'PSpline3Build: N<2!');
    end;
    
    //
    // Prepare
    //
    P.N := N;
    P.Periodic := False;
    SetLength(Tmp, N);
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    PSpline3Par(XY, N, PT, P.P);
    Assert(APSERVAreDistinct(P.P, N), 'PSpline3Build: consequent points are too close!');
    
    //
    // Build splines
    //
    if ST=0 then
    begin
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,0];
        end;
        Spline1DBuildAkima(P.P, Tmp, N, P.X);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,1];
        end;
        Spline1DBuildAkima(P.P, Tmp, N, P.Y);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,2];
        end;
        Spline1DBuildAkima(P.P, Tmp, N, P.Z);
    end;
    if ST=1 then
    begin
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,0];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N, 0, Double(0.0), P.X);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,1];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N, 0, Double(0.0), P.Y);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,2];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N, 0, Double(0.0), P.Z);
    end;
    if ST=2 then
    begin
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,0];
        end;
        Spline1DBuildCubic(P.P, Tmp, N, 0, Double(0.0), 0, Double(0.0), P.X);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,1];
        end;
        Spline1DBuildCubic(P.P, Tmp, N, 0, Double(0.0), 0, Double(0.0), P.Y);
        for i_ := 0 to N-1 do
        begin
            Tmp[i_] := XY[i_,2];
        end;
        Spline1DBuildCubic(P.P, Tmp, N, 0, Double(0.0), 0, Double(0.0), P.Z);
    end;
end;


(*************************************************************************
This  function  builds  periodic  2-dimensional  parametric  spline  which
starts at (X[0],Y[0]), goes through all points to (X[N-1],Y[N-1]) and then
back to (X[0],Y[0]).

INPUT PARAMETERS:
    XY  -   points, array[0..N-1,0..1].
            XY[I,0:1] corresponds to the Ith point.
            XY[N-1,0:1] must be different from XY[0,0:1].
            Order of points is important!
    N   -   points count, N>=3 for other types of splines.
    ST  -   spline type:
            * 1     Catmull-Rom spline (Tension=0) with cyclic boundary conditions
            * 2     cubic spline with cyclic boundary conditions
    PT  -   parameterization type:
            * 0     uniform
            * 1     chord length
            * 2     centripetal

OUTPUT PARAMETERS:
    P   -   parametric spline interpolant


NOTES:
* this function  assumes  that there all consequent points  are  distinct.
  I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so  on.
  However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
  =(x2,y2).
* last point of sequence is NOT equal to the first  point.  You  shouldn't
  make curve "explicitly periodic" by making them equal.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2BuildPeriodic(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline2Interpolant);
var
    XYP : TReal2DArray;
    Tmp : TReal1DArray;
    V : Double;
    I : AlglibInteger;
    i_ : AlglibInteger;
begin
    XY := DynamicArrayCopy(XY);
    Assert((ST>=1) and (ST<=2), 'PSpline2BuildPeriodic: incorrect spline type!');
    Assert((PT>=0) and (PT<=2), 'PSpline2BuildPeriodic: incorrect parameterization type!');
    Assert(N>=3, 'PSpline2BuildPeriodic: N<3!');
    
    //
    // Prepare
    //
    P.N := N;
    P.Periodic := True;
    SetLength(Tmp, N+1);
    SetLength(XYP, N+1, 2);
    for i_ := 0 to N-1 do
    begin
        XYP[i_,0] := XY[i_,0];
    end;
    for i_ := 0 to N-1 do
    begin
        XYP[i_,1] := XY[i_,1];
    end;
    APVMove(@XYP[N][0], 0, 1, @XY[0][0], 0, 1);
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    PSpline2Par(XYP, N+1, PT, P.P);
    Assert(APSERVAreDistinct(P.P, N+1), 'PSpline2BuildPeriodic: consequent (or first and last) points are too close!');
    
    //
    // Build splines
    //
    if ST=1 then
    begin
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,0];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N+1, -1, Double(0.0), P.X);
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,1];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N+1, -1, Double(0.0), P.Y);
    end;
    if ST=2 then
    begin
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,0];
        end;
        Spline1DBuildCubic(P.P, Tmp, N+1, -1, Double(0.0), -1, Double(0.0), P.X);
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,1];
        end;
        Spline1DBuildCubic(P.P, Tmp, N+1, -1, Double(0.0), -1, Double(0.0), P.Y);
    end;
end;


(*************************************************************************
This  function  builds  periodic  3-dimensional  parametric  spline  which
starts at (X[0],Y[0],Z[0]), goes through all points to (X[N-1],Y[N-1],Z[N-1])
and then back to (X[0],Y[0],Z[0]).

Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
description here.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3BuildPeriodic(XY : TReal2DArray;
     N : AlglibInteger;
     ST : AlglibInteger;
     PT : AlglibInteger;
     var P : PSpline3Interpolant);
var
    XYP : TReal2DArray;
    Tmp : TReal1DArray;
    V : Double;
    I : AlglibInteger;
    i_ : AlglibInteger;
begin
    XY := DynamicArrayCopy(XY);
    Assert((ST>=1) and (ST<=2), 'PSpline3BuildPeriodic: incorrect spline type!');
    Assert((PT>=0) and (PT<=2), 'PSpline3BuildPeriodic: incorrect parameterization type!');
    Assert(N>=3, 'PSpline3BuildPeriodic: N<3!');
    
    //
    // Prepare
    //
    P.N := N;
    P.Periodic := True;
    SetLength(Tmp, N+1);
    SetLength(XYP, N+1, 3);
    for i_ := 0 to N-1 do
    begin
        XYP[i_,0] := XY[i_,0];
    end;
    for i_ := 0 to N-1 do
    begin
        XYP[i_,1] := XY[i_,1];
    end;
    for i_ := 0 to N-1 do
    begin
        XYP[i_,2] := XY[i_,2];
    end;
    APVMove(@XYP[N][0], 0, 2, @XY[0][0], 0, 2);
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    PSpline3Par(XYP, N+1, PT, P.P);
    Assert(APSERVAreDistinct(P.P, N+1), 'PSplineBuild2Periodic: consequent (or first and last) points are too close!');
    
    //
    // Build splines
    //
    if ST=1 then
    begin
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,0];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N+1, -1, Double(0.0), P.X);
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,1];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N+1, -1, Double(0.0), P.Y);
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,2];
        end;
        Spline1DBuildCatmullRom(P.P, Tmp, N+1, -1, Double(0.0), P.Z);
    end;
    if ST=2 then
    begin
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,0];
        end;
        Spline1DBuildCubic(P.P, Tmp, N+1, -1, Double(0.0), -1, Double(0.0), P.X);
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,1];
        end;
        Spline1DBuildCubic(P.P, Tmp, N+1, -1, Double(0.0), -1, Double(0.0), P.Y);
        for i_ := 0 to N do
        begin
            Tmp[i_] := XYP[i_,2];
        end;
        Spline1DBuildCubic(P.P, Tmp, N+1, -1, Double(0.0), -1, Double(0.0), P.Z);
    end;
end;


(*************************************************************************
This function returns vector of parameter values correspoding to points.

I.e. for P created from (X[0],Y[0])...(X[N-1],Y[N-1]) and U=TValues(P)  we
have
    (X[0],Y[0]) = PSpline2Calc(P,U[0]),
    (X[1],Y[1]) = PSpline2Calc(P,U[1]),
    (X[2],Y[2]) = PSpline2Calc(P,U[2]),
    ...

INPUT PARAMETERS:
    P   -   parametric spline interpolant

OUTPUT PARAMETERS:
    N   -   array size
    T   -   array[0..N-1]


NOTES:
* for non-periodic splines U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]=1
* for periodic splines     U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]<1

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2ParameterValues(const P : PSpline2Interpolant;
     var N : AlglibInteger;
     var T : TReal1DArray);
begin
    Assert(P.N>=2, 'PSpline2ParameterValues: internal error!');
    N := P.N;
    SetLength(T, N);
    APVMove(@T[0], 0, N-1, @P.P[0], 0, N-1);
    T[0] := 0;
    if  not P.Periodic then
    begin
        T[N-1] := 1;
    end;
end;


(*************************************************************************
This function returns vector of parameter values correspoding to points.

Same as PSpline2ParameterValues(), but for 3D.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3ParameterValues(const P : PSpline3Interpolant;
     var N : AlglibInteger;
     var T : TReal1DArray);
begin
    Assert(P.N>=2, 'PSpline3ParameterValues: internal error!');
    N := P.N;
    SetLength(T, N);
    APVMove(@T[0], 0, N-1, @P.P[0], 0, N-1);
    T[0] := 0;
    if  not P.Periodic then
    begin
        T[N-1] := 1;
    end;
end;


(*************************************************************************
This function  calculates  the value of the parametric spline for a  given
value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-position
    Y   -   Y-position


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2Calc(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var Y : Double);
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    X := Spline1DCalc(P.X, T);
    Y := Spline1DCalc(P.Y, T);
end;


(*************************************************************************
This function  calculates  the value of the parametric spline for a  given
value of parameter T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-position
    Y   -   Y-position
    Z   -   Z-position


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3Calc(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var Y : Double;
     var Z : Double);
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    X := Spline1DCalc(P.X, T);
    Y := Spline1DCalc(P.Y, T);
    Z := Spline1DCalc(P.Z, T);
end;


(*************************************************************************
This function  calculates  tangent vector for a given value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X    -   X-component of tangent vector (normalized)
    Y    -   Y-component of tangent vector (normalized)
    
NOTE:
    X^2+Y^2 is either 1 (for non-zero tangent vector) or 0.


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2Tangent(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var Y : Double);
var
    V : Double;
    V0 : Double;
    V1 : Double;
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    PSpline2Diff(P, T, V0, X, V1, Y);
    if AP_FP_Neq(X,0) or AP_FP_Neq(Y,0) then
    begin
        
        //
        // this code is a bit more complex than X^2+Y^2 to avoid
        // overflow for large values of X and Y.
        //
        V := SafePythag2(X, Y);
        X := X/V;
        Y := Y/V;
    end;
end;


(*************************************************************************
This function  calculates  tangent vector for a given value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X    -   X-component of tangent vector (normalized)
    Y    -   Y-component of tangent vector (normalized)
    Z    -   Z-component of tangent vector (normalized)

NOTE:
    X^2+Y^2+Z^2 is either 1 (for non-zero tangent vector) or 0.


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3Tangent(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var Y : Double;
     var Z : Double);
var
    V : Double;
    V0 : Double;
    V1 : Double;
    V2 : Double;
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    PSpline3Diff(P, T, V0, X, V1, Y, V2, Z);
    if AP_FP_Neq(X,0) or AP_FP_Neq(Y,0) or AP_FP_Neq(Z,0) then
    begin
        V := SafePythag3(X, Y, Z);
        X := X/V;
        Y := Y/V;
        Z := Z/V;
    end;
end;


(*************************************************************************
This function calculates derivative, i.e. it returns (dX/dT,dY/dT).

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   X-derivative
    Y   -   Y-value
    DY  -   Y-derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2Diff(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var Y : Double;
     var DY : Double);
var
    D2S : Double;
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    Spline1DDiff(P.X, T, X, DX, D2S);
    Spline1DDiff(P.Y, T, Y, DY, D2S);
end;


(*************************************************************************
This function calculates derivative, i.e. it returns (dX/dT,dY/dT,dZ/dT).

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   X-derivative
    Y   -   Y-value
    DY  -   Y-derivative
    Z   -   Z-value
    DZ  -   Z-derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3Diff(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var Y : Double;
     var DY : Double;
     var Z : Double;
     var DZ : Double);
var
    D2S : Double;
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    Spline1DDiff(P.X, T, X, DX, D2S);
    Spline1DDiff(P.Y, T, Y, DY, D2S);
    Spline1DDiff(P.Z, T, Z, DZ, D2S);
end;


(*************************************************************************
This function calculates first and second derivative with respect to T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   derivative
    D2X -   second derivative
    Y   -   Y-value
    DY  -   derivative
    D2Y -   second derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline2Diff2(const P : PSpline2Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var D2X : Double;
     var Y : Double;
     var DY : Double;
     var D2Y : Double);
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    Spline1DDiff(P.X, T, X, DX, D2X);
    Spline1DDiff(P.Y, T, Y, DY, D2Y);
end;


(*************************************************************************
This function calculates first and second derivative with respect to T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   derivative
    D2X -   second derivative
    Y   -   Y-value
    DY  -   derivative
    D2Y -   second derivative
    Z   -   Z-value
    DZ  -   derivative
    D2Z -   second derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************)
procedure PSpline3Diff2(const P : PSpline3Interpolant;
     T : Double;
     var X : Double;
     var DX : Double;
     var D2X : Double;
     var Y : Double;
     var DY : Double;
     var D2Y : Double;
     var Z : Double;
     var DZ : Double;
     var D2Z : Double);
begin
    if P.Periodic then
    begin
        T := T-Floor(T);
    end;
    Spline1DDiff(P.X, T, X, DX, D2X);
    Spline1DDiff(P.Y, T, Y, DY, D2Y);
    Spline1DDiff(P.Z, T, Z, DZ, D2Z);
end;


(*************************************************************************
This function  calculates  arc length, i.e. length of  curve  between  t=a
and t=b.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    A,B -   parameter values corresponding to arc ends:
            * B>A will result in positive length returned
            * B<A will result in negative length returned

RESULT:
    length of arc starting at T=A and ending at T=B.


  -- ALGLIB PROJECT --
     Copyright 30.05.2010 by Bochkanov Sergey
*************************************************************************)
function PSpline2ArcLength(const P : PSpline2Interpolant;
     A : Double;
     B : Double):Double;
var
    State : AutoGKState;
    Rep : AutoGKReport;
    SX : Double;
    DSX : Double;
    D2SX : Double;
    SY : Double;
    DSY : Double;
    D2SY : Double;
begin
    AutoGKSmooth(A, B, State);
    while AutoGKIteration(State) do
    begin
        Spline1DDiff(P.X, State.X, SX, DSX, D2SX);
        Spline1DDiff(P.Y, State.X, SY, DSY, D2SY);
        State.F := SafePythag2(DSX, DSY);
    end;
    AutoGKResults(State, Result, Rep);
    Assert(Rep.TerminationType>0, 'PSpline2ArcLength: internal error!');
end;


(*************************************************************************
This function  calculates  arc length, i.e. length of  curve  between  t=a
and t=b.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    A,B -   parameter values corresponding to arc ends:
            * B>A will result in positive length returned
            * B<A will result in negative length returned

RESULT:
    length of arc starting at T=A and ending at T=B.


  -- ALGLIB PROJECT --
     Copyright 30.05.2010 by Bochkanov Sergey
*************************************************************************)
function PSpline3ArcLength(const P : PSpline3Interpolant;
     A : Double;
     B : Double):Double;
var
    State : AutoGKState;
    Rep : AutoGKReport;
    SX : Double;
    DSX : Double;
    D2SX : Double;
    SY : Double;
    DSY : Double;
    D2SY : Double;
    SZ : Double;
    DSZ : Double;
    D2SZ : Double;
begin
    AutoGKSmooth(A, B, State);
    while AutoGKIteration(State) do
    begin
        Spline1DDiff(P.X, State.X, SX, DSX, D2SX);
        Spline1DDiff(P.Y, State.X, SY, DSY, D2SY);
        Spline1DDiff(P.Z, State.X, SZ, DSZ, D2SZ);
        State.F := SafePythag3(DSX, DSY, DSZ);
    end;
    AutoGKResults(State, Result, Rep);
    Assert(Rep.TerminationType>0, 'PSpline3ArcLength: internal error!');
end;


(*************************************************************************
Builds non-periodic parameterization for 2-dimensional spline
*************************************************************************)
procedure PSpline2Par(const XY : TReal2DArray;
     N : AlglibInteger;
     PT : AlglibInteger;
     var P : TReal1DArray);
var
    V : Double;
    I : AlglibInteger;
begin
    Assert((PT>=0) and (PT<=2), 'PSpline2Par: internal error!');
    
    //
    // Build parameterization:
    // * fill by non-normalized values
    // * normalize them so we have P[0]=0, P[N-1]=1.
    //
    SetLength(P, N);
    if PT=0 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            P[I] := I;
            Inc(I);
        end;
    end;
    if PT=1 then
    begin
        P[0] := 0;
        I:=1;
        while I<=N-1 do
        begin
            P[I] := P[I-1]+SafePythag2(XY[I,0]-XY[I-1,0], XY[I,1]-XY[I-1,1]);
            Inc(I);
        end;
    end;
    if PT=2 then
    begin
        P[0] := 0;
        I:=1;
        while I<=N-1 do
        begin
            P[I] := P[I-1]+Sqrt(SafePythag2(XY[I,0]-XY[I-1,0], XY[I,1]-XY[I-1,1]));
            Inc(I);
        end;
    end;
    V := 1/P[N-1];
    APVMul(@P[0], 0, N-1, V);
end;


(*************************************************************************
Builds non-periodic parameterization for 3-dimensional spline
*************************************************************************)
procedure PSpline3Par(const XY : TReal2DArray;
     N : AlglibInteger;
     PT : AlglibInteger;
     var P : TReal1DArray);
var
    V : Double;
    I : AlglibInteger;
begin
    Assert((PT>=0) and (PT<=2), 'PSpline3Par: internal error!');
    
    //
    // Build parameterization:
    // * fill by non-normalized values
    // * normalize them so we have P[0]=0, P[N-1]=1.
    //
    SetLength(P, N);
    if PT=0 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            P[I] := I;
            Inc(I);
        end;
    end;
    if PT=1 then
    begin
        P[0] := 0;
        I:=1;
        while I<=N-1 do
        begin
            P[I] := P[I-1]+SafePythag3(XY[I,0]-XY[I-1,0], XY[I,1]-XY[I-1,1], XY[I,2]-XY[I-1,2]);
            Inc(I);
        end;
    end;
    if PT=2 then
    begin
        P[0] := 0;
        I:=1;
        while I<=N-1 do
        begin
            P[I] := P[I-1]+Sqrt(SafePythag3(XY[I,0]-XY[I-1,0], XY[I,1]-XY[I-1,1], XY[I,2]-XY[I-1,2]));
            Inc(I);
        end;
    end;
    V := 1/P[N-1];
    APVMul(@P[0], 0, N-1, V);
end;


end.