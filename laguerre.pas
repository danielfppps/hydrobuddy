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
unit laguerre;
interface
uses Math, Sysutils, Ap;

function LaguerreCalculate(const n : AlglibInteger; const x : Double):Double;
function LaguerreSum(const C : TReal1DArray;
     const n : AlglibInteger;
     const x : Double):Double;
procedure LaguerreCoefficients(const N : AlglibInteger; var C : TReal1DArray);

implementation

(*************************************************************************
Calculation of the value of the Laguerre polynomial.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Laguerre polynomial Ln at x
*************************************************************************)
function LaguerreCalculate(const n : AlglibInteger; const x : Double):Double;
var
    a : Double;
    b : Double;
    i : Double;
begin
    Result := 1;
    a := 1;
    b := 1-x;
    if n=1 then
    begin
        Result := b;
    end;
    i := 2;
    while AP_FP_Less_Eq(i,n) do
    begin
        Result := ((2*i-1-x)*b-(i-1)*a)/i;
        a := b;
        b := Result;
        i := i+1;
    end;
end;


(*************************************************************************
Summation of Laguerre polynomials using Clenshaw’s recurrence formula.

This routine calculates c[0]*L0(x) + c[1]*L1(x) + ... + c[N]*LN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Laguerre polynomial at x
*************************************************************************)
function LaguerreSum(const C : TReal1DArray;
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
        Result := (2*i+1-x)*b1/(i+1)-(i+1)*b2/(i+2)+C[i];
        b2 := b1;
        b1 := Result;
        Dec(i);
    end;
end;


(*************************************************************************
Representation of Ln as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************)
procedure LaguerreCoefficients(const N : AlglibInteger; var C : TReal1DArray);
var
    i : AlglibInteger;
begin
    SetLength(C, N+1);
    C[0] := 1;
    i:=0;
    while i<=n-1 do
    begin
        C[i+1] := -C[i]*(n-i)/(i+1)/(i+1);
        Inc(i);
    end;
end;


end.