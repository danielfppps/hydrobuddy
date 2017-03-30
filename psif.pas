{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
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
unit psif;
interface
uses Math, Sysutils, Ap;

function Psi(X : Double):Double;

implementation

(*************************************************************************
Psi (digamma) function

             d      -
  psi(x)  =  -- ln | (x)
             dx

is the logarithmic derivative of the gamma function.
For integer x,
                  n-1
                   -
psi(n) = -EUL  +   >  1/k.
                   -
                  k=1

This formula is used for 0 < n <= 10.  If x is negative, it
is transformed to a positive argument by the reflection
formula  psi(1-x) = psi(x) + pi cot(pi x).
For general positive x, the argument is made greater than 10
using the recurrence  psi(x+1) = psi(x) + 1/x.
Then the following asymptotic expansion is applied:

                          inf.   B
                           -      2k
psi(x) = log(x) - 1/2x -   >   -------
                           -        2k
                          k=1   2k x

where the B2k are Bernoulli numbers.

ACCURACY:
   Relative error (except absolute when |psi| < 1):
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        30000       1.3e-15     1.4e-16
   IEEE      -30,0       40000       1.5e-15     2.2e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
*************************************************************************)
function Psi(X : Double):Double;
var
    p : Double;
    q : Double;
    nz : Double;
    s : Double;
    w : Double;
    y : Double;
    z : Double;
    polv : Double;
    i : AlglibInteger;
    n : AlglibInteger;
    negative : AlglibInteger;
begin
    negative := 0;
    nz := Double(0.0);
    if AP_FP_Less_Eq(x,0) then
    begin
        negative := 1;
        q := x;
        p := Floor(q);
        if AP_FP_Eq(p,q) then
        begin
            Assert(False, 'Singularity in Psi(x)');
            Result := MaxRealNumber;
            Exit;
        end;
        nz := q-p;
        if AP_FP_Neq(nz,Double(0.5)) then
        begin
            if AP_FP_Greater(nz,Double(0.5)) then
            begin
                p := p+Double(1.0);
                nz := q-p;
            end;
            nz := PI/tan(PI*nz);
        end
        else
        begin
            nz := Double(0.0);
        end;
        x := Double(1.0)-x;
    end;
    if AP_FP_Less_Eq(x,Double(10.0)) and AP_FP_Eq(x,floor(x)) then
    begin
        y := Double(0.0);
        n := Floor(x);
        i:=1;
        while i<=n-1 do
        begin
            w := i;
            y := y+Double(1.0)/w;
            Inc(i);
        end;
        y := y-Double(0.57721566490153286061);
    end
    else
    begin
        s := x;
        w := Double(0.0);
        while AP_FP_Less(s,Double(10.0)) do
        begin
            w := w+Double(1.0)/s;
            s := s+Double(1.0);
        end;
        if AP_FP_Less(s,Double(1.0E17)) then
        begin
            z := Double(1.0)/(s*s);
            polv := Double(8.33333333333333333333E-2);
            polv := polv*z-Double(2.10927960927960927961E-2);
            polv := polv*z+Double(7.57575757575757575758E-3);
            polv := polv*z-Double(4.16666666666666666667E-3);
            polv := polv*z+Double(3.96825396825396825397E-3);
            polv := polv*z-Double(8.33333333333333333333E-3);
            polv := polv*z+Double(8.33333333333333333333E-2);
            y := z*polv;
        end
        else
        begin
            y := Double(0.0);
        end;
        y := ln(s)-Double(0.5)/s-y-w;
    end;
    if negative<>0 then
    begin
        y := y-nz;
    end;
    Result := y;
end;


end.