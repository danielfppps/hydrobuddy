{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
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
unit legendre;
interface
uses Math, Sysutils, Ap;

function LegendreCalculate(const n : AlglibInteger; const x : Double):Double;
function LegendreSum(const C : TReal1DArray;
     const n : AlglibInteger;
     const x : Double):Double;
procedure LegendreCoefficients(const N : AlglibInteger; var C : TReal1DArray);

implementation

(*************************************************************************
Calculation of the value of the Legendre polynomial Pn.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Legendre polynomial Pn at x
*************************************************************************)
function LegendreCalculate(const n : AlglibInteger; const x : Double):Double;
var
    a : Double;
    b : Double;
    i : AlglibInteger;
begin
    Result := 1;
    a := 1;
    b := x;
    if n=0 then
    begin
        Result := a;
        Exit;
    end;
    if n=1 then
    begin
        Result := b;
        Exit;
    end;
    i:=2;
    while i<=n do
    begin
        Result := ((2*i-1)*x*b-(i-1)*a)/i;
        a := b;
        b := Result;
        Inc(i);
    end;
end;


(*************************************************************************
Summation of Legendre polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*P0(x) + c[1]*P1(x) + ... + c[N]*PN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Legendre polynomial at x
*************************************************************************)
function LegendreSum(const C : TReal1DArray;
     const n : AlglibInteger;
     const x : Double):Double;
var
    b1 : Double;
    b2 : Double;
    i : AlglibInteger;
begin
    b1 := 0;
    b2 := 0;
    i:=n;
    while i>=0 do
    begin
        Result := (2*i+1)*x*b1/(i+1)-(i+1)*b2/(i+2)+C[i];
        b2 := b1;
        b1 := Result;
        Dec(i);
    end;
end;


(*************************************************************************
Representation of Pn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************)
procedure LegendreCoefficients(const N : AlglibInteger; var C : TReal1DArray);
var
    I : AlglibInteger;
begin
    SetLength(C, N+1);
    I:=0;
    while I<=N do
    begin
        C[I] := 0;
        Inc(I);
    end;
    C[N] := 1;
    i:=1;
    while i<=N do
    begin
        C[N] := C[N]*(n+i)/2/i;
        Inc(i);
    end;
    i:=0;
    while i<=n div 2-1 do
    begin
        C[N-2*(i+1)] := -C[N-2*i]*(n-2*i)*(n-2*i-1)/2/(i+1)/(2*(n-i)-1);
        Inc(i);
    end;
end;


end.