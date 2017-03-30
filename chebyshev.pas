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
unit chebyshev;
interface
uses Math, Sysutils, Ap;

function ChebyshevCalculate(const r : AlglibInteger;
     const N : AlglibInteger;
     const X : Double):Double;
function ChebyshevSum(const C : TReal1DArray;
     const r : AlglibInteger;
     const n : AlglibInteger;
     const x : Double):Double;
procedure ChebyshevCoefficients(const N : AlglibInteger; var C : TReal1DArray);
procedure FromChebyshev(const A : TReal1DArray;
     const N : AlglibInteger;
     var B : TReal1DArray);

implementation

(*************************************************************************
Calculation of the value of the Chebyshev polynomials of the
first and second kinds.

Parameters:
    r   -   polynomial kind, either 1 or 2.
    n   -   degree, n>=0
    x   -   argument, -1 <= x <= 1

Result:
    the value of the Chebyshev polynomial at x
*************************************************************************)
function ChebyshevCalculate(const r : AlglibInteger;
     const N : AlglibInteger;
     const X : Double):Double;
var
    I : AlglibInteger;
    A : Double;
    B : Double;
begin
    
    //
    // Prepare A and B
    //
    if r=1 then
    begin
        a := 1;
        b := x;
    end
    else
    begin
        a := 1;
        b := 2*x;
    end;
    
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
        Result := 2*x*b-a;
        a := b;
        b := Result;
        Inc(I);
    end;
end;


(*************************************************************************
Summation of Chebyshev polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*T0(x) + c[1]*T1(x) + ... + c[N]*TN(x)
or
    c[0]*U0(x) + c[1]*U1(x) + ... + c[N]*UN(x)
depending on the R.

Parameters:
    r   -   polynomial kind, either 1 or 2.
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Chebyshev polynomial at x
*************************************************************************)
function ChebyshevSum(const C : TReal1DArray;
     const r : AlglibInteger;
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
    while i>=1 do
    begin
        Result := 2*x*b1-b2+C[i];
        b2 := b1;
        b1 := Result;
        Dec(i);
    end;
    if r=1 then
    begin
        Result := -b2+x*b1+C[0];
    end
    else
    begin
        Result := -b2+2*x*b1+C[0];
    end;
end;


(*************************************************************************
Representation of Tn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************)
procedure ChebyshevCoefficients(const N : AlglibInteger; var C : TReal1DArray);
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
    if (N=0) or (N=1) then
    begin
        C[N] := 1;
    end
    else
    begin
        C[N] := Exp((n-1)*Ln(2));
        I:=0;
        while I<=n div 2-1 do
        begin
            C[N-2*(i+1)] := -C[N-2*i]*(n-2*i)*(n-2*i-1)/4/(i+1)/(n-i-1);
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Conversion of a series of Chebyshev polynomials to a power series.

Represents A[0]*T0(x) + A[1]*T1(x) + ... + A[N]*Tn(x) as
B[0] + B[1]*X + ... + B[N]*X^N.

Input parameters:
    A   -   Chebyshev series coefficients
    N   -   degree, N>=0
    
Output parameters
    B   -   power series coefficients
*************************************************************************)
procedure FromChebyshev(const A : TReal1DArray;
     const N : AlglibInteger;
     var B : TReal1DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    E : Double;
    D : Double;
begin
    SetLength(b, N+1);
    I:=0;
    while I<=N do
    begin
        b[i] := 0;
        Inc(I);
    end;
    d := 0;
    i := 0;
    repeat
        k := i;
        repeat
            e := b[k];
            b[k] := 0;
            if (i<=1) and (k=i) then
            begin
                b[k] := 1;
            end
            else
            begin
                if i<>0 then
                begin
                    b[k] := 2*d;
                end;
                if k>i+1 then
                begin
                    b[k] := b[k]-b[k-2];
                end;
            end;
            d := e;
            k := k+1;
        until  not (k<=n);
        d := b[i];
        e := 0;
        k := i;
        while k<=n do
        begin
            e := e+b[k]*a[k];
            k := k+2;
        end;
        b[i] := e;
        i := i+1;
    until  not (i<=n);
end;


end.