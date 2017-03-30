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
unit gammafunc;
interface
uses Math, Sysutils, Ap;

function Gamma(x : Double):Double;
function LnGamma(x : Double; var SgnGam : Double):Double;

implementation

function GammaStirF(X : Double):Double;forward;


(*************************************************************************
Gamma function

Input parameters:
    X   -   argument

Domain:
    0 < X < 171.6
    -170 < X < 0, X is not an integer.

Relative error:
 arithmetic   domain     # trials      peak         rms
    IEEE    -170,-33      20000       2.3e-15     3.3e-16
    IEEE     -33,  33     20000       9.4e-16     2.2e-16
    IEEE      33, 171.6   20000       2.3e-15     3.2e-16

Cephes Math Library Release 2.8:  June, 2000
Original copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************)
function Gamma(x : Double):Double;
var
    p : Double;
    PP : Double;
    q : Double;
    QQ : Double;
    z : Double;
    i : AlglibInteger;
    SgnGam : Double;
begin
    SgnGam := 1;
    q := AbsReal(x);
    if AP_FP_Greater(q,Double(33.0)) then
    begin
        if AP_FP_Less(x,Double(0.0)) then
        begin
            p := floor(q);
            i := Round(p);
            if i mod 2=0 then
            begin
                SgnGam := -1;
            end;
            z := q-p;
            if AP_FP_Greater(z,Double(0.5)) then
            begin
                p := p+1;
                z := q-p;
            end;
            z := q*Sin(Pi*z);
            z := AbsReal(z);
            z := Pi/(z*GammaStirF(q));
        end
        else
        begin
            z := GammaStirF(x);
        end;
        Result := SgnGam*z;
        Exit;
    end;
    z := 1;
    while AP_FP_Greater_Eq(x,3) do
    begin
        x := x-1;
        z := z*x;
    end;
    while AP_FP_Less(x,0) do
    begin
        if AP_FP_Greater(x,-Double(0.000000001)) then
        begin
            Result := z/((1+Double(0.5772156649015329)*x)*x);
            Exit;
        end;
        z := z/x;
        x := x+1;
    end;
    while AP_FP_Less(x,2) do
    begin
        if AP_FP_Less(x,Double(0.000000001)) then
        begin
            Result := z/((1+Double(0.5772156649015329)*x)*x);
            Exit;
        end;
        z := z/x;
        x := x+Double(1.0);
    end;
    if AP_FP_Eq(x,2) then
    begin
        Result := z;
        Exit;
    end;
    x := x-Double(2.0);
    PP := Double(1.60119522476751861407E-4);
    PP := Double(1.19135147006586384913E-3)+X*PP;
    PP := Double(1.04213797561761569935E-2)+X*PP;
    PP := Double(4.76367800457137231464E-2)+X*PP;
    PP := Double(2.07448227648435975150E-1)+X*PP;
    PP := Double(4.94214826801497100753E-1)+X*PP;
    PP := Double(9.99999999999999996796E-1)+X*PP;
    QQ := -Double(2.31581873324120129819E-5);
    QQ := Double(5.39605580493303397842E-4)+X*QQ;
    QQ := -Double(4.45641913851797240494E-3)+X*QQ;
    QQ := Double(1.18139785222060435552E-2)+X*QQ;
    QQ := Double(3.58236398605498653373E-2)+X*QQ;
    QQ := -Double(2.34591795718243348568E-1)+X*QQ;
    QQ := Double(7.14304917030273074085E-2)+X*QQ;
    QQ := Double(1.00000000000000000320)+X*QQ;
    Result := z*PP/QQ;
    Exit;
end;


(*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************)
function LnGamma(x : Double; var SgnGam : Double):Double;
var
    A : Double;
    B : Double;
    C : Double;
    p : Double;
    q : Double;
    u : Double;
    w : Double;
    z : Double;
    i : AlglibInteger;
    LogPi : Double;
    LS2PI : Double;
    Tmp : Double;
begin
    SgnGam := 1;
    LogPi := Double(1.14472988584940017414);
    LS2PI := Double(0.91893853320467274178);
    if AP_FP_Less(x,-Double(34.0)) then
    begin
        q := -x;
        w := LnGamma(q, Tmp);
        p := floor(q);
        i := Round(p);
        if i mod 2=0 then
        begin
            SgnGam := -1;
        end
        else
        begin
            SgnGam := 1;
        end;
        z := q-p;
        if AP_FP_Greater(z,Double(0.5)) then
        begin
            p := p+1;
            z := p-q;
        end;
        z := q*Sin(Pi*z);
        Result := LogPi-Ln(z)-w;
        Exit;
    end;
    if AP_FP_Less(x,13) then
    begin
        z := 1;
        p := 0;
        u := x;
        while AP_FP_Greater_Eq(u,3) do
        begin
            p := p-1;
            u := x+p;
            z := z*u;
        end;
        while AP_FP_Less(u,2) do
        begin
            z := z/u;
            p := p+1;
            u := x+p;
        end;
        if AP_FP_Less(z,0) then
        begin
            sgngam := -1;
            z := -z;
        end
        else
        begin
            sgngam := 1;
        end;
        if AP_FP_Eq(u,2) then
        begin
            Result := Ln(z);
            Exit;
        end;
        p := p-2;
        x := x+p;
        B := -Double(1378.25152569120859100);
        B := -Double(38801.6315134637840924)+X*B;
        B := -Double(331612.992738871184744)+X*B;
        B := -Double(1162370.97492762307383)+X*B;
        B := -Double(1721737.00820839662146)+X*B;
        B := -Double(853555.664245765465627)+X*B;
        C := 1;
        C := -Double(351.815701436523470549)+X*C;
        C := -Double(17064.2106651881159223)+X*C;
        C := -Double(220528.590553854454839)+X*C;
        C := -Double(1139334.44367982507207)+X*C;
        C := -Double(2532523.07177582951285)+X*C;
        C := -Double(2018891.41433532773231)+X*C;
        p := x*B/C;
        Result := Ln(z)+p;
        Exit;
    end;
    q := (x-Double(0.5))*Ln(x)-x+LS2PI;
    if AP_FP_Greater(x,100000000) then
    begin
        Result := q;
        Exit;
    end;
    p := 1/(x*x);
    if AP_FP_Greater_Eq(x,Double(1000.0)) then
    begin
        q := q+((Double(7.9365079365079365079365)*Double(0.0001)*p-Double(2.7777777777777777777778)*Double(0.001))*p+Double(0.0833333333333333333333))/x;
    end
    else
    begin
        A := Double(8.11614167470508450300)*Double(0.0001);
        A := -Double(5.95061904284301438324)*Double(0.0001)+p*A;
        A := Double(7.93650340457716943945)*Double(0.0001)+p*A;
        A := -Double(2.77777777730099687205)*Double(0.001)+p*A;
        A := Double(8.33333333333331927722)*Double(0.01)+p*A;
        q := q+A/x;
    end;
    Result := q;
end;


function GammaStirF(X : Double):Double;
var
    y : Double;
    w : Double;
    v : Double;
    Stir : Double;
begin
    w := 1/x;
    Stir := Double(7.87311395793093628397E-4);
    Stir := -Double(2.29549961613378126380E-4)+w*Stir;
    Stir := -Double(2.68132617805781232825E-3)+w*Stir;
    Stir := Double(3.47222221605458667310E-3)+w*Stir;
    Stir := Double(8.33333333333482257126E-2)+w*Stir;
    w := 1+w*Stir;
    y := Exp(x);
    if AP_FP_Greater(x,Double(143.01608)) then
    begin
        v := Power(x, Double(0.5)*x-Double(0.25));
        y := v*(v/y);
    end
    else
    begin
        y := Power(x, x-Double(0.5))/y;
    end;
    Result := Double(2.50662827463100050242)*y*w;
end;


end.