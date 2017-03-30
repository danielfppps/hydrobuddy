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
unit betaf;
interface
uses Math, Sysutils, Ap, gammafunc;

function Beta(a : Double; b : Double):Double;

implementation

(*************************************************************************
Beta function


                  -     -
                 | (a) | (b)
beta( a, b )  =  -----------.
                    -
                   | (a+b)

For large arguments the logarithm of the function is
evaluated using lgam(), then exponentiated.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,30       30000       8.1e-14     1.1e-14

Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
*************************************************************************)
function Beta(a : Double; b : Double):Double;
var
    y : Double;
    sg : Double;
    S : Double;
begin
    sg := 1;
    Assert(AP_FP_Greater(a,0) or AP_FP_Neq(a,Floor(a)), 'Overflow in Beta');
    Assert(AP_FP_Greater(b,0) or AP_FP_Neq(b,Floor(b)), 'Overflow in Beta');
    y := a+b;
    if AP_FP_Greater(absReal(y),Double(171.624376956302725)) then
    begin
        y := LnGamma(y, S);
        sg := sg*S;
        y := LnGamma(b, S)-y;
        sg := sg*S;
        y := LnGamma(a, S)+y;
        sg := sg*S;
        Assert(AP_FP_Less_Eq(y,Ln(MaxRealNumber)), 'Overflow in Beta');
        Result := sg*Exp(y);
        Exit;
    end;
    y := Gamma(y);
    Assert(AP_FP_Neq(y,0), 'Overflow in Beta');
    if AP_FP_Greater(a,b) then
    begin
        y := Gamma(a)/y;
        y := y*Gamma(b);
    end
    else
    begin
        y := Gamma(b)/y;
        y := y*Gamma(a);
    end;
    Result := y;
end;


end.