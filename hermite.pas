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
unit hermite;
interface
uses Math, Sysutils, Ap;

function HermiteCalculate(const N : AlglibInteger; const X : Double):Double;
function HermiteSum(const C : TReal1DArray;
     const n : AlglibInteger;
     const x : Double):Double;
procedure HermiteCoefficients(const N : AlglibInteger; var C : TReal1DArray);

implementation

(*************************************************************************
Calculation of the value of the Hermite polynomial.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Hermite polynomial Hn at x
*************************************************************************)
function HermiteCalculate(const N : AlglibInteger; const X : Double):Double;
var
    I : AlglibInteger;
    a : Double;
    b : Double;
begin
    
    //
    // Prepare A and B
    //
    a := 1;
    b := 2*x;
    
    //
    // Special cases: N=0 or N=1
    //
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
    
    //
    // General case: N>=2
    //
    I:=2;
    while I<=N do
    begin
        Result := 2*x*b-2*(I-1)*a;
        a := b;
        b := Result;
        Inc(I);
    end;
end;


(*************************************************************************
Summation of Hermite polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*H0(x) + c[1]*H1(x) + ... + c[N]*HN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Hermite polynomial at x
*************************************************************************)
function HermiteSum(const C : TReal1DArray;
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
        Result := 2*(x*b1-(i+1)*b2)+C[i];
        b2 := b1;
        b1 := Result;
        Dec(i);
    end;
end;


(*************************************************************************
Representation of Hn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************)
procedure HermiteCoefficients(const N : AlglibInteger; var C : TReal1DArray);
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
    C[N] := Exp(N*Ln(2));
    i:=0;
    while i<=n div 2-1 do
    begin
        C[N-2*(i+1)] := -C[N-2*i]*(n-2*i)*(n-2*i-1)/4/(i+1);
        Inc(i);
    end;
end;


end.