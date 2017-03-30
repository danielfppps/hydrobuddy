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
unit bessel;
interface
uses Math, Sysutils, Ap;

function BesselJ0(X : Double):Double;
function BesselJ1(X : Double):Double;
function BesselJN(n : AlglibInteger; x : Double):Double;
function BesselY0(X : Double):Double;
function BesselY1(X : Double):Double;
function BesselYN(N : AlglibInteger; X : Double):Double;
function BesselI0(X : Double):Double;
function BesselI1(x : Double):Double;
function BesselK0(X : Double):Double;
function BesselK1(x : Double):Double;
function BesselKN(nn : AlglibInteger; x : Double):Double;

implementation

procedure BesselMFirstCheb(c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);forward;
procedure BesselMNextCheb(x : Double;
     c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);forward;
procedure BesselM1FirstCheb(c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);forward;
procedure BesselM1NextCheb(x : Double;
     c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);forward;
procedure BesselAsympt0(X : Double;
     var PZero : Double;
     var QZero : Double);forward;
procedure BesselAsympt1(X : Double;
     var PZero : Double;
     var QZero : Double);forward;


(*************************************************************************
Bessel function of order zero

Returns Bessel function of order zero of the argument.

The domain is divided into the intervals [0, 5] and
(5, infinity). In the first interval the following rational
approximation is used:


       2         2
(w - r  ) (w - r  ) P (w) / Q (w)
      1         2    3       8

           2
where w = x  and the two r's are zeros of the function.

In the second interval, the Hankel asymptotic expansion
is employed with two rational functions of degree 6/6
and 7/7.

ACCURACY:

                     Absolute error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       60000       4.2e-16     1.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselJ0(X : Double):Double;
var
    XSq : Double;
    NN : Double;
    PZero : Double;
    QZero : Double;
    P1 : Double;
    Q1 : Double;
begin
    if AP_FP_Less(X,0) then
    begin
        X := -X;
    end;
    if AP_FP_Greater(X,Double(8.0)) then
    begin
        BesselAsympt0(X, PZero, QZero);
        NN := X-Pi/4;
        Result := Sqrt(2/Pi/X)*(PZero*Cos(NN)-QZero*Sin(NN));
        Exit;
    end;
    XSq := AP_Sqr(X);
    P1 := Double(26857.86856980014981415848441);
    P1 := -Double(40504123.71833132706360663322)+XSq*P1;
    P1 := Double(25071582855.36881945555156435)+XSq*P1;
    P1 := -Double(8085222034853.793871199468171)+XSq*P1;
    P1 := Double(1434354939140344.111664316553)+XSq*P1;
    P1 := -Double(136762035308817138.6865416609)+XSq*P1;
    P1 := Double(6382059341072356562.289432465)+XSq*P1;
    P1 := -Double(117915762910761053603.8440800)+XSq*P1;
    P1 := Double(493378725179413356181.6813446)+XSq*P1;
    Q1 := Double(1.0);
    Q1 := Double(1363.063652328970604442810507)+XSq*Q1;
    Q1 := Double(1114636.098462985378182402543)+XSq*Q1;
    Q1 := Double(669998767.2982239671814028660)+XSq*Q1;
    Q1 := Double(312304311494.1213172572469442)+XSq*Q1;
    Q1 := Double(112775673967979.8507056031594)+XSq*Q1;
    Q1 := Double(30246356167094626.98627330784)+XSq*Q1;
    Q1 := Double(5428918384092285160.200195092)+XSq*Q1;
    Q1 := Double(493378725179413356211.3278438)+XSq*Q1;
    Result := P1/Q1;
end;


(*************************************************************************
Bessel function of order one

Returns Bessel function of order one of the argument.

The domain is divided into the intervals [0, 8] and
(8, infinity). In the first interval a 24 term Chebyshev
expansion is used. In the second, the asymptotic
trigonometric representation is employed using two
rational functions of degree 5/5.

ACCURACY:

                     Absolute error:
arithmetic   domain      # trials      peak         rms
   IEEE      0, 30       30000       2.6e-16     1.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselJ1(X : Double):Double;
var
    S : Double;
    XSq : Double;
    NN : Double;
    PZero : Double;
    QZero : Double;
    P1 : Double;
    Q1 : Double;
begin
    S := Sign(X);
    if AP_FP_Less(X,0) then
    begin
        X := -X;
    end;
    if AP_FP_Greater(X,Double(8.0)) then
    begin
        BesselAsympt1(X, PZero, QZero);
        NN := X-3*Pi/4;
        Result := Sqrt(2/Pi/X)*(PZero*Cos(NN)-QZero*Sin(NN));
        if AP_FP_Less(S,0) then
        begin
            Result := -Result;
        end;
        Exit;
    end;
    XSq := AP_Sqr(X);
    P1 := Double(2701.122710892323414856790990);
    P1 := -Double(4695753.530642995859767162166)+XSq*P1;
    P1 := Double(3413234182.301700539091292655)+XSq*P1;
    P1 := -Double(1322983480332.126453125473247)+XSq*P1;
    P1 := Double(290879526383477.5409737601689)+XSq*P1;
    P1 := -Double(35888175699101060.50743641413)+XSq*P1;
    P1 := Double(2316433580634002297.931815435)+XSq*P1;
    P1 := -Double(66721065689249162980.20941484)+XSq*P1;
    P1 := Double(581199354001606143928.050809)+XSq*P1;
    Q1 := Double(1.0);
    Q1 := Double(1606.931573481487801970916749)+XSq*Q1;
    Q1 := Double(1501793.594998585505921097578)+XSq*Q1;
    Q1 := Double(1013863514.358673989967045588)+XSq*Q1;
    Q1 := Double(524371026216.7649715406728642)+XSq*Q1;
    Q1 := Double(208166122130760.7351240184229)+XSq*Q1;
    Q1 := Double(60920613989175217.46105196863)+XSq*Q1;
    Q1 := Double(11857707121903209998.37113348)+XSq*Q1;
    Q1 := Double(1162398708003212287858.529400)+XSq*Q1;
    Result := S*X*P1/Q1;
end;


(*************************************************************************
Bessel function of integer order

Returns Bessel function of order n, where n is a
(possibly negative) integer.

The ratio of jn(x) to j0(x) is computed by backward
recurrence.  First the ratio jn/jn-1 is found by a
continued fraction expansion.  Then the recurrence
relating successive orders is applied until j0 or j1 is
reached.

If n = 0 or 1 the routine for j0 or j1 is called
directly.

ACCURACY:

                     Absolute error:
arithmetic   range      # trials      peak         rms
   IEEE      0, 30        5000       4.4e-16     7.9e-17


Not suitable for large n or x. Use jv() (fractional order) instead.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselJN(n : AlglibInteger; x : Double):Double;
var
    pkm2 : Double;
    pkm1 : Double;
    pk : Double;
    xk : Double;
    r : Double;
    ans : Double;
    k : AlglibInteger;
    sg : AlglibInteger;
begin
    if n<0 then
    begin
        n := -n;
        if n mod 2=0 then
        begin
            sg := 1;
        end
        else
        begin
            sg := -1;
        end;
    end
    else
    begin
        sg := 1;
    end;
    if AP_FP_Less(x,0) then
    begin
        if n mod 2<>0 then
        begin
            sg := -sg;
        end;
        x := -x;
    end;
    if n=0 then
    begin
        Result := sg*BesselJ0(x);
        Exit;
    end;
    if n=1 then
    begin
        Result := sg*BesselJ1(x);
        Exit;
    end;
    if n=2 then
    begin
        if AP_FP_Eq(x,0) then
        begin
            Result := 0;
        end
        else
        begin
            Result := sg*(Double(2.0)*BesselJ1(x)/x-BesselJ0(x));
        end;
        Exit;
    end;
    if AP_FP_Less(x,MachineEpsilon) then
    begin
        Result := 0;
        Exit;
    end;
    k := 53;
    pk := 2*(n+k);
    ans := pk;
    xk := x*x;
    repeat
        pk := pk-Double(2.0);
        ans := pk-xk/ans;
        k := k-1;
    until k=0;
    ans := x/ans;
    pk := Double(1.0);
    pkm1 := Double(1.0)/ans;
    k := n-1;
    r := 2*k;
    repeat
        pkm2 := (pkm1*r-pk*x)/x;
        pk := pkm1;
        pkm1 := pkm2;
        r := r-Double(2.0);
        k := k-1;
    until k=0;
    if AP_FP_Greater(AbsReal(pk),AbsReal(pkm1)) then
    begin
        ans := BesselJ1(x)/pk;
    end
    else
    begin
        ans := BesselJ0(x)/pkm1;
    end;
    Result := sg*ans;
end;


(*************************************************************************
Bessel function of the second kind, order zero

Returns Bessel function of the second kind, of order
zero, of the argument.

The domain is divided into the intervals [0, 5] and
(5, infinity). In the first interval a rational approximation
R(x) is employed to compute
  y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
Thus a call to j0() is required.

In the second interval, the Hankel asymptotic expansion
is employed with two rational functions of degree 6/6
and 7/7.



ACCURACY:

 Absolute error, when y0(x) < 1; else relative error:

arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.3e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselY0(X : Double):Double;
var
    NN : Double;
    XSq : Double;
    PZero : Double;
    QZero : Double;
    P4 : Double;
    Q4 : Double;
begin
    if AP_FP_Greater(X,Double(8.0)) then
    begin
        BesselAsympt0(X, PZero, QZero);
        NN := X-Pi/4;
        Result := Sqrt(2/Pi/X)*(PZero*Sin(NN)+QZero*Cos(NN));
        Exit;
    end;
    XSq := AP_Sqr(X);
    P4 := -Double(41370.35497933148554125235152);
    P4 := Double(59152134.65686889654273830069)+XSq*P4;
    P4 := -Double(34363712229.79040378171030138)+XSq*P4;
    P4 := Double(10255208596863.94284509167421)+XSq*P4;
    P4 := -Double(1648605817185729.473122082537)+XSq*P4;
    P4 := Double(137562431639934407.8571335453)+XSq*P4;
    P4 := -Double(5247065581112764941.297350814)+XSq*P4;
    P4 := Double(65874732757195549259.99402049)+XSq*P4;
    P4 := -Double(27502866786291095837.01933175)+XSq*P4;
    Q4 := Double(1.0);
    Q4 := Double(1282.452772478993804176329391)+XSq*Q4;
    Q4 := Double(1001702.641288906265666651753)+XSq*Q4;
    Q4 := Double(579512264.0700729537480087915)+XSq*Q4;
    Q4 := Double(261306575504.1081249568482092)+XSq*Q4;
    Q4 := Double(91620380340751.85262489147968)+XSq*Q4;
    Q4 := Double(23928830434997818.57439356652)+XSq*Q4;
    Q4 := Double(4192417043410839973.904769661)+XSq*Q4;
    Q4 := Double(372645883898616588198.9980)+XSq*Q4;
    Result := P4/Q4+2/Pi*BesselJ0(X)*Ln(X);
end;


(*************************************************************************
Bessel function of second kind of order one

Returns Bessel function of the second kind of order one
of the argument.

The domain is divided into the intervals [0, 8] and
(8, infinity). In the first interval a 25 term Chebyshev
expansion is used, and a call to j1() is required.
In the second, the asymptotic trigonometric representation
is employed using two rational functions of degree 5/5.

ACCURACY:

                     Absolute error:
arithmetic   domain      # trials      peak         rms
   IEEE      0, 30       30000       1.0e-15     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselY1(X : Double):Double;
var
    NN : Double;
    XSq : Double;
    PZero : Double;
    QZero : Double;
    P4 : Double;
    Q4 : Double;
begin
    if AP_FP_Greater(X,Double(8.0)) then
    begin
        BesselAsympt1(X, PZero, QZero);
        NN := X-3*Pi/4;
        Result := Sqrt(2/Pi/X)*(PZero*Sin(NN)+QZero*Cos(NN));
        Exit;
    end;
    XSq := AP_Sqr(X);
    P4 := -Double(2108847.540133123652824139923);
    P4 := Double(3639488548.124002058278999428)+XSq*P4;
    P4 := -Double(2580681702194.450950541426399)+XSq*P4;
    P4 := Double(956993023992168.3481121552788)+XSq*P4;
    P4 := -Double(196588746272214065.8820322248)+XSq*P4;
    P4 := Double(21931073399177975921.11427556)+XSq*P4;
    P4 := -Double(1212297555414509577913.561535)+XSq*P4;
    P4 := Double(26554738314348543268942.48968)+XSq*P4;
    P4 := -Double(99637534243069222259967.44354)+XSq*P4;
    Q4 := Double(1.0);
    Q4 := Double(1612.361029677000859332072312)+XSq*Q4;
    Q4 := Double(1563282.754899580604737366452)+XSq*Q4;
    Q4 := Double(1128686837.169442121732366891)+XSq*Q4;
    Q4 := Double(646534088126.5275571961681500)+XSq*Q4;
    Q4 := Double(297663212564727.6729292742282)+XSq*Q4;
    Q4 := Double(108225825940881955.2553850180)+XSq*Q4;
    Q4 := Double(29549879358971486742.90758119)+XSq*Q4;
    Q4 := Double(5435310377188854170800.653097)+XSq*Q4;
    Q4 := Double(508206736694124324531442.4152)+XSq*Q4;
    Result := X*P4/Q4+2/Pi*(BesselJ1(X)*Ln(X)-1/X);
end;


(*************************************************************************
Bessel function of second kind of integer order

Returns Bessel function of order n, where n is a
(possibly negative) integer.

The function is evaluated by forward recurrence on
n, starting with values computed by the routines
y0() and y1().

If n = 0 or 1 the routine for y0 or y1 is called
directly.

ACCURACY:
                     Absolute error, except relative
                     when y > 1:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       3.4e-15     4.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselYN(N : AlglibInteger; X : Double):Double;
var
    I : AlglibInteger;
    A : Double;
    B : Double;
    Tmp : Double;
    S : Double;
begin
    S := 1;
    if N<0 then
    begin
        N := -N;
        if N mod 2<>0 then
        begin
            S := -1;
        end;
    end;
    if N=0 then
    begin
        Result := BesselY0(X);
        Exit;
    end;
    if N=1 then
    begin
        Result := S*BesselY1(X);
        Exit;
    end;
    A := BesselY0(X);
    B := BesselY1(X);
    I:=1;
    while I<=N-1 do
    begin
        Tmp := B;
        B := 2*I/X*B-A;
        A := Tmp;
        Inc(I);
    end;
    Result := S*B;
end;


(*************************************************************************
Modified Bessel function of order zero

Returns modified Bessel function of order zero of the
argument.

The function is defined as i0(x) = j0( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        30000       5.8e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselI0(X : Double):Double;
var
    y : Double;
    v : Double;
    z : Double;
    b0 : Double;
    b1 : Double;
    b2 : Double;
begin
    if AP_FP_Less(x,0) then
    begin
        x := -x;
    end;
    if AP_FP_Less_Eq(x,Double(8.0)) then
    begin
        y := x/Double(2.0)-Double(2.0);
        BesselMFirstCheb(-Double(4.41534164647933937950E-18), b0, b1, b2);
        BesselMNextCheb(y, Double(3.33079451882223809783E-17), b0, b1, b2);
        BesselMNextCheb(y, -Double(2.43127984654795469359E-16), b0, b1, b2);
        BesselMNextCheb(y, Double(1.71539128555513303061E-15), b0, b1, b2);
        BesselMNextCheb(y, -Double(1.16853328779934516808E-14), b0, b1, b2);
        BesselMNextCheb(y, Double(7.67618549860493561688E-14), b0, b1, b2);
        BesselMNextCheb(y, -Double(4.85644678311192946090E-13), b0, b1, b2);
        BesselMNextCheb(y, Double(2.95505266312963983461E-12), b0, b1, b2);
        BesselMNextCheb(y, -Double(1.72682629144155570723E-11), b0, b1, b2);
        BesselMNextCheb(y, Double(9.67580903537323691224E-11), b0, b1, b2);
        BesselMNextCheb(y, -Double(5.18979560163526290666E-10), b0, b1, b2);
        BesselMNextCheb(y, Double(2.65982372468238665035E-9), b0, b1, b2);
        BesselMNextCheb(y, -Double(1.30002500998624804212E-8), b0, b1, b2);
        BesselMNextCheb(y, Double(6.04699502254191894932E-8), b0, b1, b2);
        BesselMNextCheb(y, -Double(2.67079385394061173391E-7), b0, b1, b2);
        BesselMNextCheb(y, Double(1.11738753912010371815E-6), b0, b1, b2);
        BesselMNextCheb(y, -Double(4.41673835845875056359E-6), b0, b1, b2);
        BesselMNextCheb(y, Double(1.64484480707288970893E-5), b0, b1, b2);
        BesselMNextCheb(y, -Double(5.75419501008210370398E-5), b0, b1, b2);
        BesselMNextCheb(y, Double(1.88502885095841655729E-4), b0, b1, b2);
        BesselMNextCheb(y, -Double(5.76375574538582365885E-4), b0, b1, b2);
        BesselMNextCheb(y, Double(1.63947561694133579842E-3), b0, b1, b2);
        BesselMNextCheb(y, -Double(4.32430999505057594430E-3), b0, b1, b2);
        BesselMNextCheb(y, Double(1.05464603945949983183E-2), b0, b1, b2);
        BesselMNextCheb(y, -Double(2.37374148058994688156E-2), b0, b1, b2);
        BesselMNextCheb(y, Double(4.93052842396707084878E-2), b0, b1, b2);
        BesselMNextCheb(y, -Double(9.49010970480476444210E-2), b0, b1, b2);
        BesselMNextCheb(y, Double(1.71620901522208775349E-1), b0, b1, b2);
        BesselMNextCheb(y, -Double(3.04682672343198398683E-1), b0, b1, b2);
        BesselMNextCheb(y, Double(6.76795274409476084995E-1), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        Result := exp(x)*v;
        Exit;
    end;
    z := Double(32.0)/x-Double(2.0);
    BesselMFirstCheb(-Double(7.23318048787475395456E-18), b0, b1, b2);
    BesselMNextCheb(z, -Double(4.83050448594418207126E-18), b0, b1, b2);
    BesselMNextCheb(z, Double(4.46562142029675999901E-17), b0, b1, b2);
    BesselMNextCheb(z, Double(3.46122286769746109310E-17), b0, b1, b2);
    BesselMNextCheb(z, -Double(2.82762398051658348494E-16), b0, b1, b2);
    BesselMNextCheb(z, -Double(3.42548561967721913462E-16), b0, b1, b2);
    BesselMNextCheb(z, Double(1.77256013305652638360E-15), b0, b1, b2);
    BesselMNextCheb(z, Double(3.81168066935262242075E-15), b0, b1, b2);
    BesselMNextCheb(z, -Double(9.55484669882830764870E-15), b0, b1, b2);
    BesselMNextCheb(z, -Double(4.15056934728722208663E-14), b0, b1, b2);
    BesselMNextCheb(z, Double(1.54008621752140982691E-14), b0, b1, b2);
    BesselMNextCheb(z, Double(3.85277838274214270114E-13), b0, b1, b2);
    BesselMNextCheb(z, Double(7.18012445138366623367E-13), b0, b1, b2);
    BesselMNextCheb(z, -Double(1.79417853150680611778E-12), b0, b1, b2);
    BesselMNextCheb(z, -Double(1.32158118404477131188E-11), b0, b1, b2);
    BesselMNextCheb(z, -Double(3.14991652796324136454E-11), b0, b1, b2);
    BesselMNextCheb(z, Double(1.18891471078464383424E-11), b0, b1, b2);
    BesselMNextCheb(z, Double(4.94060238822496958910E-10), b0, b1, b2);
    BesselMNextCheb(z, Double(3.39623202570838634515E-9), b0, b1, b2);
    BesselMNextCheb(z, Double(2.26666899049817806459E-8), b0, b1, b2);
    BesselMNextCheb(z, Double(2.04891858946906374183E-7), b0, b1, b2);
    BesselMNextCheb(z, Double(2.89137052083475648297E-6), b0, b1, b2);
    BesselMNextCheb(z, Double(6.88975834691682398426E-5), b0, b1, b2);
    BesselMNextCheb(z, Double(3.36911647825569408990E-3), b0, b1, b2);
    BesselMNextCheb(z, Double(8.04490411014108831608E-1), b0, b1, b2);
    v := Double(0.5)*(b0-b2);
    Result := exp(x)*v/sqrt(x);
end;


(*************************************************************************
Modified Bessel function of order one

Returns modified Bessel function of order one of the
argument.

The function is defined as i1(x) = -i j1( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.9e-15     2.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselI1(x : Double):Double;
var
    y : Double;
    z : Double;
    v : Double;
    b0 : Double;
    b1 : Double;
    b2 : Double;
begin
    z := AbsReal(x);
    if AP_FP_Less_Eq(z,Double(8.0)) then
    begin
        y := z/Double(2.0)-Double(2.0);
        BesselM1FirstCheb(Double(2.77791411276104639959E-18), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.11142121435816608115E-17), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.55363195773620046921E-16), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.10559694773538630805E-15), b0, b1, b2);
        BesselM1NextCheb(y, Double(7.60068429473540693410E-15), b0, b1, b2);
        BesselM1NextCheb(y, -Double(5.04218550472791168711E-14), b0, b1, b2);
        BesselM1NextCheb(y, Double(3.22379336594557470981E-13), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.98397439776494371520E-12), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.17361862988909016308E-11), b0, b1, b2);
        BesselM1NextCheb(y, -Double(6.66348972350202774223E-11), b0, b1, b2);
        BesselM1NextCheb(y, Double(3.62559028155211703701E-10), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.88724975172282928790E-9), b0, b1, b2);
        BesselM1NextCheb(y, Double(9.38153738649577178388E-9), b0, b1, b2);
        BesselM1NextCheb(y, -Double(4.44505912879632808065E-8), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.00329475355213526229E-7), b0, b1, b2);
        BesselM1NextCheb(y, -Double(8.56872026469545474066E-7), b0, b1, b2);
        BesselM1NextCheb(y, Double(3.47025130813767847674E-6), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.32731636560394358279E-5), b0, b1, b2);
        BesselM1NextCheb(y, Double(4.78156510755005422638E-5), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.61760815825896745588E-4), b0, b1, b2);
        BesselM1NextCheb(y, Double(5.12285956168575772895E-4), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.51357245063125314899E-3), b0, b1, b2);
        BesselM1NextCheb(y, Double(4.15642294431288815669E-3), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.05640848946261981558E-2), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.47264490306265168283E-2), b0, b1, b2);
        BesselM1NextCheb(y, -Double(5.29459812080949914269E-2), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.02643658689847095384E-1), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.76416518357834055153E-1), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.52587186443633654823E-1), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        z := v*z*exp(z);
    end
    else
    begin
        y := Double(32.0)/z-Double(2.0);
        BesselM1FirstCheb(Double(7.51729631084210481353E-18), b0, b1, b2);
        BesselM1NextCheb(y, Double(4.41434832307170791151E-18), b0, b1, b2);
        BesselM1NextCheb(y, -Double(4.65030536848935832153E-17), b0, b1, b2);
        BesselM1NextCheb(y, -Double(3.20952592199342395980E-17), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.96262899764595013876E-16), b0, b1, b2);
        BesselM1NextCheb(y, Double(3.30820231092092828324E-16), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.88035477551078244854E-15), b0, b1, b2);
        BesselM1NextCheb(y, -Double(3.81440307243700780478E-15), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.04202769841288027642E-14), b0, b1, b2);
        BesselM1NextCheb(y, Double(4.27244001671195135429E-14), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.10154184277266431302E-14), b0, b1, b2);
        BesselM1NextCheb(y, -Double(4.08355111109219731823E-13), b0, b1, b2);
        BesselM1NextCheb(y, -Double(7.19855177624590851209E-13), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.03562854414708950722E-12), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.41258074366137813316E-11), b0, b1, b2);
        BesselM1NextCheb(y, Double(3.25260358301548823856E-11), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.89749581235054123450E-11), b0, b1, b2);
        BesselM1NextCheb(y, -Double(5.58974346219658380687E-10), b0, b1, b2);
        BesselM1NextCheb(y, -Double(3.83538038596423702205E-9), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.63146884688951950684E-8), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.51223623787020892529E-7), b0, b1, b2);
        BesselM1NextCheb(y, -Double(3.88256480887769039346E-6), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.10588938762623716291E-4), b0, b1, b2);
        BesselM1NextCheb(y, -Double(9.76109749136146840777E-3), b0, b1, b2);
        BesselM1NextCheb(y, Double(7.78576235018280120474E-1), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        z := v*Exp(z)/Sqrt(z);
    end;
    if AP_FP_Less(x,0) then
    begin
        z := -z;
    end;
    Result := z;
end;


(*************************************************************************
Modified Bessel function, second kind, order zero

Returns modified Bessel function of the second kind
of order zero of the argument.

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

Tested at 2000 random points between 0 and 8.  Peak absolute
error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.2e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselK0(X : Double):Double;
var
    y : Double;
    z : Double;
    v : Double;
    b0 : Double;
    b1 : Double;
    b2 : Double;
begin
    Assert(AP_FP_Greater(x,0), 'Domain error in BesselK0: x<=0');
    if AP_FP_Less_Eq(x,2) then
    begin
        y := x*x-Double(2.0);
        BesselMFirstCheb(Double(1.37446543561352307156E-16), b0, b1, b2);
        BesselMNextCheb(y, Double(4.25981614279661018399E-14), b0, b1, b2);
        BesselMNextCheb(y, Double(1.03496952576338420167E-11), b0, b1, b2);
        BesselMNextCheb(y, Double(1.90451637722020886025E-9), b0, b1, b2);
        BesselMNextCheb(y, Double(2.53479107902614945675E-7), b0, b1, b2);
        BesselMNextCheb(y, Double(2.28621210311945178607E-5), b0, b1, b2);
        BesselMNextCheb(y, Double(1.26461541144692592338E-3), b0, b1, b2);
        BesselMNextCheb(y, Double(3.59799365153615016266E-2), b0, b1, b2);
        BesselMNextCheb(y, Double(3.44289899924628486886E-1), b0, b1, b2);
        BesselMNextCheb(y, -Double(5.35327393233902768720E-1), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        v := v-Ln(Double(0.5)*x)*BesselI0(x);
    end
    else
    begin
        z := Double(8.0)/x-Double(2.0);
        BesselMFirstCheb(Double(5.30043377268626276149E-18), b0, b1, b2);
        BesselMNextCheb(z, -Double(1.64758043015242134646E-17), b0, b1, b2);
        BesselMNextCheb(z, Double(5.21039150503902756861E-17), b0, b1, b2);
        BesselMNextCheb(z, -Double(1.67823109680541210385E-16), b0, b1, b2);
        BesselMNextCheb(z, Double(5.51205597852431940784E-16), b0, b1, b2);
        BesselMNextCheb(z, -Double(1.84859337734377901440E-15), b0, b1, b2);
        BesselMNextCheb(z, Double(6.34007647740507060557E-15), b0, b1, b2);
        BesselMNextCheb(z, -Double(2.22751332699166985548E-14), b0, b1, b2);
        BesselMNextCheb(z, Double(8.03289077536357521100E-14), b0, b1, b2);
        BesselMNextCheb(z, -Double(2.98009692317273043925E-13), b0, b1, b2);
        BesselMNextCheb(z, Double(1.14034058820847496303E-12), b0, b1, b2);
        BesselMNextCheb(z, -Double(4.51459788337394416547E-12), b0, b1, b2);
        BesselMNextCheb(z, Double(1.85594911495471785253E-11), b0, b1, b2);
        BesselMNextCheb(z, -Double(7.95748924447710747776E-11), b0, b1, b2);
        BesselMNextCheb(z, Double(3.57739728140030116597E-10), b0, b1, b2);
        BesselMNextCheb(z, -Double(1.69753450938905987466E-9), b0, b1, b2);
        BesselMNextCheb(z, Double(8.57403401741422608519E-9), b0, b1, b2);
        BesselMNextCheb(z, -Double(4.66048989768794782956E-8), b0, b1, b2);
        BesselMNextCheb(z, Double(2.76681363944501510342E-7), b0, b1, b2);
        BesselMNextCheb(z, -Double(1.83175552271911948767E-6), b0, b1, b2);
        BesselMNextCheb(z, Double(1.39498137188764993662E-5), b0, b1, b2);
        BesselMNextCheb(z, -Double(1.28495495816278026384E-4), b0, b1, b2);
        BesselMNextCheb(z, Double(1.56988388573005337491E-3), b0, b1, b2);
        BesselMNextCheb(z, -Double(3.14481013119645005427E-2), b0, b1, b2);
        BesselMNextCheb(z, Double(2.44030308206595545468E0), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        v := v*Exp(-x)/Sqrt(x);
    end;
    Result := v;
end;


(*************************************************************************
Modified Bessel function, second kind, order one

Computes the modified Bessel function of the second kind
of order one of the argument.

The range is partitioned into the two intervals [0,2] and
(2, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.2e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselK1(x : Double):Double;
var
    y : Double;
    z : Double;
    v : Double;
    b0 : Double;
    b1 : Double;
    b2 : Double;
begin
    z := Double(0.5)*x;
    Assert(AP_FP_Greater(z,0), 'Domain error in K1');
    if AP_FP_Less_Eq(x,2) then
    begin
        y := x*x-Double(2.0);
        BesselM1FirstCheb(-Double(7.02386347938628759343E-18), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.42744985051936593393E-15), b0, b1, b2);
        BesselM1NextCheb(y, -Double(6.66690169419932900609E-13), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.41148839263352776110E-10), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.21338763073472585583E-8), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.43340614156596823496E-6), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.73028895751305206302E-4), b0, b1, b2);
        BesselM1NextCheb(y, -Double(6.97572385963986435018E-3), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.22611180822657148235E-1), b0, b1, b2);
        BesselM1NextCheb(y, -Double(3.53155960776544875667E-1), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.52530022733894777053E0), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        Result := Ln(z)*BesselI1(x)+v/x;
    end
    else
    begin
        y := Double(8.0)/x-Double(2.0);
        BesselM1FirstCheb(-Double(5.75674448366501715755E-18), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.79405087314755922667E-17), b0, b1, b2);
        BesselM1NextCheb(y, -Double(5.68946255844285935196E-17), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.83809354436663880070E-16), b0, b1, b2);
        BesselM1NextCheb(y, -Double(6.05704724837331885336E-16), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.03870316562433424052E-15), b0, b1, b2);
        BesselM1NextCheb(y, -Double(7.01983709041831346144E-15), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.47715442448130437068E-14), b0, b1, b2);
        BesselM1NextCheb(y, -Double(8.97670518232499435011E-14), b0, b1, b2);
        BesselM1NextCheb(y, Double(3.34841966607842919884E-13), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.28917396095102890680E-12), b0, b1, b2);
        BesselM1NextCheb(y, Double(5.13963967348173025100E-12), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.12996783842756842877E-11), b0, b1, b2);
        BesselM1NextCheb(y, Double(9.21831518760500529508E-11), b0, b1, b2);
        BesselM1NextCheb(y, -Double(4.19035475934189648750E-10), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.01504975519703286596E-9), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.03457624656780970260E-8), b0, b1, b2);
        BesselM1NextCheb(y, Double(5.74108412545004946722E-8), b0, b1, b2);
        BesselM1NextCheb(y, -Double(3.50196060308781257119E-7), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.40648494783721712015E-6), b0, b1, b2);
        BesselM1NextCheb(y, -Double(1.93619797416608296024E-5), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.95215518471351631108E-4), b0, b1, b2);
        BesselM1NextCheb(y, -Double(2.85781685962277938680E-3), b0, b1, b2);
        BesselM1NextCheb(y, Double(1.03923736576817238437E-1), b0, b1, b2);
        BesselM1NextCheb(y, Double(2.72062619048444266945E0), b0, b1, b2);
        v := Double(0.5)*(b0-b2);
        Result := exp(-x)*v/sqrt(x);
    end;
end;


(*************************************************************************
Modified Bessel function, second kind, integer order

Returns modified Bessel function of the second kind
of order n of the argument.

The range is partitioned into the two intervals [0,9.55] and
(9.55, infinity).  An ascending power series is used in the
low range, and an asymptotic expansion in the high range.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        90000       1.8e-8      3.0e-10

Error is high only near the crossover point x = 9.55
between the two expansions used.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
*************************************************************************)
function BesselKN(nn : AlglibInteger; x : Double):Double;
var
    k : Double;
    kf : Double;
    nk1f : Double;
    nkf : Double;
    zn : Double;
    t : Double;
    s : Double;
    z0 : Double;
    z : Double;
    ans : Double;
    fn : Double;
    pn : Double;
    pk : Double;
    zmn : Double;
    tlg : Double;
    tox : Double;
    i : AlglibInteger;
    n : AlglibInteger;
    EUL : Double;
begin
    EUL := Double(5.772156649015328606065e-1);
    if nn<0 then
    begin
        n := -nn;
    end
    else
    begin
        n := nn;
    end;
    Assert(n<=31, 'Overflow in BesselKN');
    Assert(AP_FP_Greater(x,0), 'Domain error in BesselKN');
    if AP_FP_Less_Eq(x,Double(9.55)) then
    begin
        ans := Double(0.0);
        z0 := Double(0.25)*x*x;
        fn := Double(1.0);
        pn := Double(0.0);
        zmn := Double(1.0);
        tox := Double(2.0)/x;
        if n>0 then
        begin
            pn := -EUL;
            k := Double(1.0);
            i:=1;
            while i<=n-1 do
            begin
                pn := pn+Double(1.0)/k;
                k := k+Double(1.0);
                fn := fn*k;
                Inc(i);
            end;
            zmn := tox;
            if n=1 then
            begin
                ans := Double(1.0)/x;
            end
            else
            begin
                nk1f := fn/n;
                kf := Double(1.0);
                s := nk1f;
                z := -z0;
                zn := Double(1.0);
                i:=1;
                while i<=n-1 do
                begin
                    nk1f := nk1f/(n-i);
                    kf := kf*i;
                    zn := zn*z;
                    t := nk1f*zn/kf;
                    s := s+t;
                    Assert(AP_FP_Greater(MaxRealNumber-absReal(t),absReal(s)), 'Overflow in BesselKN');
                    Assert( not (AP_FP_Greater(tox,Double(1.0)) and AP_FP_Less(MaxRealNumber/tox,zmn)), 'Overflow in BesselKN');
                    zmn := zmn*tox;
                    Inc(i);
                end;
                s := s*Double(0.5);
                t := absReal(s);
                Assert( not (AP_FP_Greater(zmn,Double(1.0)) and AP_FP_Less(MaxRealNumber/zmn,t)), 'Overflow in BesselKN');
                Assert( not (AP_FP_Greater(t,Double(1.0)) and AP_FP_Less(MaxRealNumber/t,zmn)), 'Overflow in BesselKN');
                ans := s*zmn;
            end;
        end;
        tlg := Double(2.0)*ln(Double(0.5)*x);
        pk := -EUL;
        if n=0 then
        begin
            pn := pk;
            t := Double(1.0);
        end
        else
        begin
            pn := pn+Double(1.0)/n;
            t := Double(1.0)/fn;
        end;
        s := (pk+pn-tlg)*t;
        k := Double(1.0);
        repeat
            t := t*(z0/(k*(k+n)));
            pk := pk+Double(1.0)/k;
            pn := pn+Double(1.0)/(k+n);
            s := s+(pk+pn-tlg)*t;
            k := k+Double(1.0);
        until AP_FP_Less_Eq(absReal(t/s),MachineEpsilon);
        s := Double(0.5)*s/zmn;
        if n mod 2<>0 then
        begin
            s := -s;
        end;
        ans := ans+s;
        Result := ans;
        Exit;
    end;
    if AP_FP_Greater(x,Ln(MaxRealNumber)) then
    begin
        Result := 0;
        Exit;
    end;
    k := n;
    pn := Double(4.0)*k*k;
    pk := Double(1.0);
    z0 := Double(8.0)*x;
    fn := Double(1.0);
    t := Double(1.0);
    s := t;
    nkf := MaxRealNumber;
    i := 0;
    repeat
        z := pn-pk*pk;
        t := t*z/(fn*z0);
        nk1f := absReal(t);
        if (i>=n) and AP_FP_Greater(nk1f,nkf) then
        begin
            Break;
        end;
        nkf := nk1f;
        s := s+t;
        fn := fn+Double(1.0);
        pk := pk+Double(2.0);
        i := i+1;
    until AP_FP_Less_Eq(absReal(t/s),MachineEpsilon);
    Result := exp(-x)*sqrt(PI/(Double(2.0)*x))*s;
end;


(*************************************************************************
Internal subroutine

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
procedure BesselMFirstCheb(c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);
begin
    b0 := c;
    b1 := Double(0.0);
    b2 := Double(0.0);
end;


(*************************************************************************
Internal subroutine

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
procedure BesselMNextCheb(x : Double;
     c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);
begin
    b2 := b1;
    b1 := b0;
    b0 := x*b1-b2+c;
end;


(*************************************************************************
Internal subroutine

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
procedure BesselM1FirstCheb(c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);
begin
    b0 := c;
    b1 := Double(0.0);
    b2 := Double(0.0);
end;


(*************************************************************************
Internal subroutine

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
procedure BesselM1NextCheb(x : Double;
     c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);
begin
    b2 := b1;
    b1 := b0;
    b0 := x*b1-b2+c;
end;


procedure BesselAsympt0(X : Double; var PZero : Double; var QZero : Double);
var
    XSq : Double;
    P2 : Double;
    Q2 : Double;
    P3 : Double;
    Q3 : Double;
begin
    XSq := Double(64.0)/(X*X);
    P2 := Double(0.0);
    P2 := Double(2485.271928957404011288128951)+XSq*P2;
    P2 := Double(153982.6532623911470917825993)+XSq*P2;
    P2 := Double(2016135.283049983642487182349)+XSq*P2;
    P2 := Double(8413041.456550439208464315611)+XSq*P2;
    P2 := Double(12332384.76817638145232406055)+XSq*P2;
    P2 := Double(5393485.083869438325262122897)+XSq*P2;
    Q2 := Double(1.0);
    Q2 := Double(2615.700736920839685159081813)+XSq*Q2;
    Q2 := Double(156001.7276940030940592769933)+XSq*Q2;
    Q2 := Double(2025066.801570134013891035236)+XSq*Q2;
    Q2 := Double(8426449.050629797331554404810)+XSq*Q2;
    Q2 := Double(12338310.22786324960844856182)+XSq*Q2;
    Q2 := Double(5393485.083869438325560444960)+XSq*Q2;
    P3 := -Double(0.0);
    P3 := -Double(4.887199395841261531199129300)+XSq*P3;
    P3 := -Double(226.2630641933704113967255053)+XSq*P3;
    P3 := -Double(2365.956170779108192723612816)+XSq*P3;
    P3 := -Double(8239.066313485606568803548860)+XSq*P3;
    P3 := -Double(10381.41698748464093880530341)+XSq*P3;
    P3 := -Double(3984.617357595222463506790588)+XSq*P3;
    Q3 := Double(1.0);
    Q3 := Double(408.7714673983499223402830260)+XSq*Q3;
    Q3 := Double(15704.89191515395519392882766)+XSq*Q3;
    Q3 := Double(156021.3206679291652539287109)+XSq*Q3;
    Q3 := Double(533291.3634216897168722255057)+XSq*Q3;
    Q3 := Double(666745.4239319826986004038103)+XSq*Q3;
    Q3 := Double(255015.5108860942382983170882)+XSq*Q3;
    Pzero := P2/Q2;
    Qzero := 8*P3/Q3/X;
end;


procedure BesselAsympt1(X : Double; var PZero : Double; var QZero : Double);
var
    XSq : Double;
    P2 : Double;
    Q2 : Double;
    P3 : Double;
    Q3 : Double;
begin
    XSq := Double(64.0)/(X*X);
    P2 := -Double(1611.616644324610116477412898);
    P2 := -Double(109824.0554345934672737413139)+XSq*P2;
    P2 := -Double(1523529.351181137383255105722)+XSq*P2;
    P2 := -Double(6603373.248364939109255245434)+XSq*P2;
    P2 := -Double(9942246.505077641195658377899)+XSq*P2;
    P2 := -Double(4435757.816794127857114720794)+XSq*P2;
    Q2 := Double(1.0);
    Q2 := -Double(1455.009440190496182453565068)+XSq*Q2;
    Q2 := -Double(107263.8599110382011903063867)+XSq*Q2;
    Q2 := -Double(1511809.506634160881644546358)+XSq*Q2;
    Q2 := -Double(6585339.479723087072826915069)+XSq*Q2;
    Q2 := -Double(9934124.389934585658967556309)+XSq*Q2;
    Q2 := -Double(4435757.816794127856828016962)+XSq*Q2;
    P3 := Double(35.26513384663603218592175580);
    P3 := Double(1706.375429020768002061283546)+XSq*P3;
    P3 := Double(18494.26287322386679652009819)+XSq*P3;
    P3 := Double(66178.83658127083517939992166)+XSq*P3;
    P3 := Double(85145.16067533570196555001171)+XSq*P3;
    P3 := Double(33220.91340985722351859704442)+XSq*P3;
    Q3 := Double(1.0);
    Q3 := Double(863.8367769604990967475517183)+XSq*Q3;
    Q3 := Double(37890.22974577220264142952256)+XSq*Q3;
    Q3 := Double(400294.4358226697511708610813)+XSq*Q3;
    Q3 := Double(1419460.669603720892855755253)+XSq*Q3;
    Q3 := Double(1819458.042243997298924553839)+XSq*Q3;
    Q3 := Double(708712.8194102874357377502472)+XSq*Q3;
    Pzero := P2/Q2;
    Qzero := 8*P3/Q3/X;
end;


end.