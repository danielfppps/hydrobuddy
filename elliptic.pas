{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier

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
unit elliptic;
interface
uses Math, Sysutils, Ap;

function EllipticIntegralK(m : Double):Double;
function EllipticIntegralKHighPrecision(m1 : Double):Double;
function IncompleteEllipticIntegralK(phi : Double; m : Double):Double;
function EllipticIntegralE(m : Double):Double;
function IncompleteEllipticIntegralE(phi : Double; m : Double):Double;

implementation

(*************************************************************************
Complete elliptic integral of the first kind

Approximates the integral



           pi/2
            -
           | |
           |           dt
K(m)  =    |    ------------------
           |                   2
         | |    sqrt( 1 - m sin t )
          -
           0

using the approximation

    P(x)  -  log x Q(x).

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,1        30000       2.5e-16     6.8e-17

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function EllipticIntegralK(m : Double):Double;
begin
    Result := EllipticIntegralKHighPrecision(Double(1.0)-m);
end;


(*************************************************************************
Complete elliptic integral of the first kind

Approximates the integral



           pi/2
            -
           | |
           |           dt
K(m)  =    |    ------------------
           |                   2
         | |    sqrt( 1 - m sin t )
          -
           0

where m = 1 - m1, using the approximation

    P(x)  -  log x Q(x).

The argument m1 is used rather than m so that the logarithmic
singularity at m = 1 will be shifted to the origin; this
preserves maximum accuracy.

K(0) = pi/2.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,1        30000       2.5e-16     6.8e-17

Алгоритм взят из библиотеки Cephes
*************************************************************************)
function EllipticIntegralKHighPrecision(m1 : Double):Double;
var
    P : Double;
    Q : Double;
begin
    if AP_FP_Less_Eq(m1,MachineEpsilon) then
    begin
        Result := Double(1.3862943611198906188E0)-Double(0.5)*Ln(m1);
    end
    else
    begin
        P := Double(1.37982864606273237150E-4);
        P := P*m1+Double(2.28025724005875567385E-3);
        P := P*m1+Double(7.97404013220415179367E-3);
        P := P*m1+Double(9.85821379021226008714E-3);
        P := P*m1+Double(6.87489687449949877925E-3);
        P := P*m1+Double(6.18901033637687613229E-3);
        P := P*m1+Double(8.79078273952743772254E-3);
        P := P*m1+Double(1.49380448916805252718E-2);
        P := P*m1+Double(3.08851465246711995998E-2);
        P := P*m1+Double(9.65735902811690126535E-2);
        P := P*m1+Double(1.38629436111989062502E0);
        Q := Double(2.94078955048598507511E-5);
        Q := Q*m1+Double(9.14184723865917226571E-4);
        Q := Q*m1+Double(5.94058303753167793257E-3);
        Q := Q*m1+Double(1.54850516649762399335E-2);
        Q := Q*m1+Double(2.39089602715924892727E-2);
        Q := Q*m1+Double(3.01204715227604046988E-2);
        Q := Q*m1+Double(3.73774314173823228969E-2);
        Q := Q*m1+Double(4.88280347570998239232E-2);
        Q := Q*m1+Double(7.03124996963957469739E-2);
        Q := Q*m1+Double(1.24999999999870820058E-1);
        Q := Q*m1+Double(4.99999999999999999821E-1);
        Result := P-Q*Ln(m1);
    end;
end;


(*************************************************************************
Incomplete elliptic integral of the first kind F(phi|m)

Approximates the integral



               phi
                -
               | |
               |           dt
F(phi_\m)  =    |    ------------------
               |                   2
             | |    sqrt( 1 - m sin t )
              -
               0

of amplitude phi and modulus m, using the arithmetic -
geometric mean algorithm.




ACCURACY:

Tested at random points with m in [0, 1] and phi as indicated.

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -10,10       200000      7.4e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteEllipticIntegralK(phi : Double; m : Double):Double;
var
    a : Double;
    b : Double;
    c : Double;
    e : Double;
    temp : Double;
    PIO2 : Double;
    t : Double;
    K : Double;
    d : AlglibInteger;
    md : AlglibInteger;
    s : AlglibInteger;
    npio2 : AlglibInteger;
begin
    PIO2 := Double(1.57079632679489661923);
    if AP_FP_Eq(m,0) then
    begin
        Result := phi;
        Exit;
    end;
    a := 1-m;
    if AP_FP_Eq(a,0) then
    begin
        Result := Ln(Tan(Double(0.5)*(PIO2+phi)));
        Exit;
    end;
    npio2 := Floor(phi/PIO2);
    if npio2 mod 2<>0 then
    begin
        npio2 := npio2+1;
    end;
    if npio2<>0 then
    begin
        K := EllipticIntegralK(1-a);
        phi := phi-npio2*PIO2;
    end
    else
    begin
        K := 0;
    end;
    if AP_FP_Less(phi,0) then
    begin
        phi := -phi;
        s := -1;
    end
    else
    begin
        s := 0;
    end;
    b := sqrt(a);
    t := tan(phi);
    if AP_FP_Greater(AbsReal(t),10) then
    begin
        e := Double(1.0)/(b*t);
        if AP_FP_Less(AbsReal(e),10) then
        begin
            e := arctan(e);
            if npio2=0 then
            begin
                K := EllipticIntegralK(1-a);
            end;
            temp := K-IncompleteEllipticIntegralK(e, m);
            if s<0 then
            begin
                temp := -temp;
            end;
            Result := temp+npio2*K;
            Exit;
        end;
    end;
    a := Double(1.0);
    c := sqrt(m);
    d := 1;
    md := 0;
    while AP_FP_Greater(AbsReal(c/a),MachineEpsilon) do
    begin
        temp := b/a;
        phi := phi+arctan(t*temp)+md*Pi;
        md := Trunc((phi+PIO2)/Pi);
        t := t*(Double(1.0)+temp)/(Double(1.0)-temp*t*t);
        c := Double(0.5)*(a-b);
        temp := sqrt(a*b);
        a := Double(0.5)*(a+b);
        b := temp;
        d := d+d;
    end;
    temp := (arctan(t)+md*Pi)/(d*a);
    if s<0 then
    begin
        temp := -temp;
    end;
    Result := temp+npio2*K;
end;


(*************************************************************************
Complete elliptic integral of the second kind

Approximates the integral


           pi/2
            -
           | |                 2
E(m)  =    |    sqrt( 1 - m sin t ) dt
         | |
          -
           0

using the approximation

     P(x)  -  x log x Q(x).

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0, 1       10000       2.1e-16     7.3e-17

Cephes Math Library, Release 2.8: June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
function EllipticIntegralE(m : Double):Double;
var
    P : Double;
    Q : Double;
begin
    Assert(AP_FP_Greater_Eq(m,0) and AP_FP_Less_Eq(m,1), 'Domain error in EllipticIntegralE: m<0 or m>1');
    m := 1-m;
    if AP_FP_Eq(m,0) then
    begin
        Result := 1;
        Exit;
    end;
    P := Double(1.53552577301013293365E-4);
    P := P*m+Double(2.50888492163602060990E-3);
    P := P*m+Double(8.68786816565889628429E-3);
    P := P*m+Double(1.07350949056076193403E-2);
    P := P*m+Double(7.77395492516787092951E-3);
    P := P*m+Double(7.58395289413514708519E-3);
    P := P*m+Double(1.15688436810574127319E-2);
    P := P*m+Double(2.18317996015557253103E-2);
    P := P*m+Double(5.68051945617860553470E-2);
    P := P*m+Double(4.43147180560990850618E-1);
    P := P*m+Double(1.00000000000000000299E0);
    Q := Double(3.27954898576485872656E-5);
    Q := Q*m+Double(1.00962792679356715133E-3);
    Q := Q*m+Double(6.50609489976927491433E-3);
    Q := Q*m+Double(1.68862163993311317300E-2);
    Q := Q*m+Double(2.61769742454493659583E-2);
    Q := Q*m+Double(3.34833904888224918614E-2);
    Q := Q*m+Double(4.27180926518931511717E-2);
    Q := Q*m+Double(5.85936634471101055642E-2);
    Q := Q*m+Double(9.37499997197644278445E-2);
    Q := Q*m+Double(2.49999999999888314361E-1);
    Result := P-Q*m*Ln(m);
end;


(*************************************************************************
Incomplete elliptic integral of the second kind

Approximates the integral


               phi
                -
               | |
               |                   2
E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
               |
             | |
              -
               0

of amplitude phi and modulus m, using the arithmetic -
geometric mean algorithm.

ACCURACY:

Tested at random arguments with phi in [-10, 10] and m in
[0, 1].
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -10,10      150000       3.3e-15     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1993, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteEllipticIntegralE(phi : Double; m : Double):Double;
var
    PIO2 : Double;
    a : Double;
    b : Double;
    c : Double;
    e : Double;
    temp : Double;
    lphi : Double;
    t : Double;
    EBig : Double;
    d : AlglibInteger;
    md : AlglibInteger;
    npio2 : AlglibInteger;
    s : AlglibInteger;
begin
    PIO2 := Double(1.57079632679489661923);
    if AP_FP_Eq(m,0) then
    begin
        Result := phi;
        Exit;
    end;
    lphi := phi;
    npio2 := Floor(lphi/PIO2);
    if npio2 mod 2<>0 then
    begin
        npio2 := npio2+1;
    end;
    lphi := lphi-npio2*PIO2;
    if AP_FP_Less(lphi,0) then
    begin
        lphi := -lphi;
        s := -1;
    end
    else
    begin
        s := 1;
    end;
    a := Double(1.0)-m;
    EBig := EllipticIntegralE(m);
    if AP_FP_Eq(a,0) then
    begin
        temp := sin(lphi);
        if s<0 then
        begin
            temp := -temp;
        end;
        Result := temp+npio2*Ebig;
        Exit;
    end;
    t := tan(lphi);
    b := sqrt(a);
    
    //
    // Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
    // for pointing out an instability near odd multiples of pi/2
    //
    if AP_FP_Greater(AbsReal(t),10) then
    begin
        
        //
        // Transform the amplitude
        //
        e := Double(1.0)/(b*t);
        
        //
        // ... but avoid multiple recursions.
        //
        if AP_FP_Less(AbsReal(e),10) then
        begin
            e := arctan(e);
            temp := EBig+m*sin(lphi)*sin(e)-IncompleteEllipticIntegralE(e, m);
            if s<0 then
            begin
                temp := -temp;
            end;
            Result := temp+npio2*Ebig;
            Exit;
        end;
    end;
    c := sqrt(m);
    a := Double(1.0);
    d := 1;
    e := Double(0.0);
    md := 0;
    while AP_FP_Greater(AbsReal(c/a),MachineEpsilon) do
    begin
        temp := b/a;
        lphi := lphi+arctan(t*temp)+md*PI;
        md := Trunc((lphi+PIO2)/PI);
        t := t*(Double(1.0)+temp)/(Double(1.0)-temp*t*t);
        c := Double(0.5)*(a-b);
        temp := sqrt(a*b);
        a := Double(0.5)*(a+b);
        b := temp;
        d := d+d;
        e := e+c*sin(lphi);
    end;
    temp := EBig/EllipticIntegralK(m);
    temp := temp*((arctan(t)+md*PI)/(d*a));
    temp := temp+e;
    if s<0 then
    begin
        temp := -temp;
    end;
    Result := temp+npio2*EBig;
    Exit;
end;


end.