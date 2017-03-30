{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
unit apserv;
interface
uses Math, Sysutils, Ap;

procedure TaskGenInt1D(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
procedure TaskGenInt1DEquidist(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
procedure TaskGenInt1DCheb1(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
procedure TaskGenInt1DCheb2(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
function APSERVAreDistinct(X : TReal1DArray; N : AlglibInteger):Boolean;
function SafePythag2(X : Double; Y : Double):Double;
function SafePythag3(X : Double; Y : Double; Z : Double):Double;
procedure APPeriodicMap(var X : Double;
     A : Double;
     B : Double;
     var K : Double);

implementation

(*************************************************************************
This  function  generates  1-dimensional  general  interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1D(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    H : Double;
begin
    Assert(N>=1, 'TaskGenInterpolationEqdist1D: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        X[0] := A;
        Y[0] := 2*RandomReal-1;
        H := (B-A)/(N-1);
        I:=1;
        while I<=N-1 do
        begin
            if I<>N-1 then
            begin
                X[I] := A+(I+Double(0.2)*(2*RandomReal-1))*H;
            end
            else
            begin
                X[I] := B;
            end;
            Y[I] := Y[I-1]+(2*RandomReal-1)*(X[I]-X[I-1]);
            Inc(I);
        end;
    end
    else
    begin
        X[0] := Double(0.5)*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function generates  1-dimensional equidistant interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1DEquidist(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    H : Double;
begin
    Assert(N>=1, 'TaskGenInterpolationEqdist1D: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        X[0] := A;
        Y[0] := 2*RandomReal-1;
        H := (B-A)/(N-1);
        I:=1;
        while I<=N-1 do
        begin
            X[I] := A+I*H;
            Y[I] := Y[I-1]+(2*RandomReal-1)*H;
            Inc(I);
        end;
    end
    else
    begin
        X[0] := Double(0.5)*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function generates  1-dimensional Chebyshev-1 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1DCheb1(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
begin
    Assert(N>=1, 'TaskGenInterpolation1DCheb1: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            X[I] := Double(0.5)*(B+A)+Double(0.5)*(B-A)*Cos(Pi*(2*i+1)/(2*N));
            if I=0 then
            begin
                Y[I] := 2*RandomReal-1;
            end
            else
            begin
                Y[I] := Y[I-1]+(2*RandomReal-1)*(X[I]-X[I-1]);
            end;
            Inc(I);
        end;
    end
    else
    begin
        X[0] := Double(0.5)*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function generates  1-dimensional Chebyshev-2 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1DCheb2(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
begin
    Assert(N>=1, 'TaskGenInterpolation1DCheb2: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            X[I] := Double(0.5)*(B+A)+Double(0.5)*(B-A)*Cos(Pi*i/(n-1));
            if I=0 then
            begin
                Y[I] := 2*RandomReal-1;
            end
            else
            begin
                Y[I] := Y[I-1]+(2*RandomReal-1)*(X[I]-X[I-1]);
            end;
            Inc(I);
        end;
    end
    else
    begin
        X[0] := Double(0.5)*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function checks that all values from X[] are distinct. It does more
than just usual floating point comparison:
* first, it calculates max(X) and min(X)
* second, it maps X[] from [min,max] to [1,2]
* only at this stage actual comparison is done

The meaning of such check is to ensure that all values are "distinct enough"
and will not cause interpolation subroutine to fail.

NOTE:
    X[] must be sorted by ascending (subroutine ASSERT's it)

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function APSERVAreDistinct(X : TReal1DArray; N : AlglibInteger):Boolean;
var
    IsSorted : Boolean;
    A : Double;
    B : Double;
    I : AlglibInteger;
begin
    X := DynamicArrayCopy(X);
    Assert(N>=1, 'APSERVIsDistinct: internal error!');
    if N=1 then
    begin
        
        //
        // everything is alright, it is up to caller to decide whether it
        // can interpolate something with just one point
        //
        Result := True;
        Exit;
    end;
    A := X[0];
    B := X[0];
    I:=1;
    while I<=N-1 do
    begin
        A := Min(A, X[I]);
        B := Max(B, X[I]);
        Assert(AP_FP_Greater_Eq(X[I],X[I-1]), 'APSERVIsDistinct: internal error!');
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        X[I] := (X[I]-A)/(B-A)+1;
        Inc(I);
    end;
    I:=1;
    while I<=N-1 do
    begin
        if AP_FP_Eq(X[I],X[I-1]) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    Result := True;
end;


(*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
function SafePythag2(X : Double; Y : Double):Double;
var
    W : Double;
    XABS : Double;
    YABS : Double;
    Z : Double;
begin
    XABS := AbsReal(X);
    YABS := AbsReal(Y);
    W := Max(XABS, YABS);
    Z := Min(XABS, YABS);
    if AP_FP_Eq(Z,0) then
    begin
        Result := W;
    end
    else
    begin
        Result := W*SQRT(1+AP_Sqr(Z/W));
    end;
end;


(*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
function SafePythag3(X : Double; Y : Double; Z : Double):Double;
var
    W : Double;
begin
    W := Max(AbsReal(X), Max(AbsReal(Y), AbsReal(Z)));
    if AP_FP_Eq(W,0) then
    begin
        Result := 0;
        Exit;
    end;
    X := X/W;
    Y := Y/W;
    Z := Z/W;
    Result := W*Sqrt(AP_Sqr(X)+AP_Sqr(Y)+AP_Sqr(Z));
end;


(*************************************************************************
This function makes periodic mapping of X to [A,B].

It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
such that X = T + K*(B-A).

NOTES:
* K is represented as real value, although actually it is integer
* T is guaranteed to be in [A,B]
* T replaces X

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure APPeriodicMap(var X : Double;
     A : Double;
     B : Double;
     var K : Double);
begin
    Assert(AP_FP_Less(A,B), 'APPeriodicMap: internal error!');
    K := Floor((X-A)/(B-A));
    X := X-K*(B-A);
    while AP_FP_Less(X,A) do
    begin
        X := X+(B-A);
        K := K-1;
    end;
    while AP_FP_Greater(X,B) do
    begin
        X := X-(B-A);
        K := K+1;
    end;
    X := Max(X, A);
    X := Min(X, B);
end;


end.