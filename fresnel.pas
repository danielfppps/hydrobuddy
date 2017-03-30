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
unit fresnel;
interface
uses Math, Sysutils, Ap;

procedure FresnelIntegral(X : Double; var C : Double; var S : Double);

implementation

(*************************************************************************
Fresnel integral

Evaluates the Fresnel integrals

          x
          -
         | |
C(x) =   |   cos(pi/2 t**2) dt,
       | |
        -
         0

          x
          -
         | |
S(x) =   |   sin(pi/2 t**2) dt.
       | |
        -
         0


The integrals are evaluated by a power series for x < 1.
For x >= 1 auxiliary functions f(x) and g(x) are employed
such that

C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )



ACCURACY:

 Relative error.

Arithmetic  function   domain     # trials      peak         rms
  IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
  IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
procedure FresnelIntegral(X : Double; var C : Double; var S : Double);
var
    XXA : Double;
    F : Double;
    G : Double;
    CC : Double;
    SS : Double;
    T : Double;
    U : Double;
    X2 : Double;
    SN : Double;
    SD : Double;
    CN : Double;
    CD : Double;
    FN : Double;
    FD : Double;
    GN : Double;
    GD : Double;
    MPI : Double;
    MPIO2 : Double;
begin
    MPI := Double(3.14159265358979323846);
    MPIO2 := Double(1.57079632679489661923);
    XXA := X;
    X := AbsReal(XXA);
    X2 := X*X;
    if AP_FP_Less(X2,Double(2.5625)) then
    begin
        T := x2*x2;
        SN := -Double(2.99181919401019853726E3);
        SN := SN*T+Double(7.08840045257738576863E5);
        SN := SN*T-Double(6.29741486205862506537E7);
        SN := SN*T+Double(2.54890880573376359104E9);
        SN := SN*T-Double(4.42979518059697779103E10);
        SN := SN*T+Double(3.18016297876567817986E11);
        SD := Double(1.00000000000000000000E0);
        SD := SD*T+Double(2.81376268889994315696E2);
        SD := SD*T+Double(4.55847810806532581675E4);
        SD := SD*T+Double(5.17343888770096400730E6);
        SD := SD*T+Double(4.19320245898111231129E8);
        SD := SD*T+Double(2.24411795645340920940E10);
        SD := SD*T+Double(6.07366389490084639049E11);
        CN := -Double(4.98843114573573548651E-8);
        CN := CN*T+Double(9.50428062829859605134E-6);
        CN := CN*T-Double(6.45191435683965050962E-4);
        CN := CN*T+Double(1.88843319396703850064E-2);
        CN := CN*T-Double(2.05525900955013891793E-1);
        CN := CN*T+Double(9.99999999999999998822E-1);
        CD := Double(3.99982968972495980367E-12);
        CD := CD*T+Double(9.15439215774657478799E-10);
        CD := CD*T+Double(1.25001862479598821474E-7);
        CD := CD*T+Double(1.22262789024179030997E-5);
        CD := CD*T+Double(8.68029542941784300606E-4);
        CD := CD*T+Double(4.12142090722199792936E-2);
        CD := CD*T+Double(1.00000000000000000118E0);
        S := Sign(XXA)*x*x2*SN/SD;
        C := Sign(XXA)*x*CN/CD;
        Exit;
    end;
    if AP_FP_Greater(x,Double(36974.0)) then
    begin
        c := Sign(XXA)*Double(0.5);
        s := Sign(XXA)*Double(0.5);
        Exit;
    end;
    x2 := x*x;
    t := MPI*x2;
    u := 1/(t*t);
    t := 1/t;
    FN := Double(4.21543555043677546506E-1);
    FN := FN*U+Double(1.43407919780758885261E-1);
    FN := FN*U+Double(1.15220955073585758835E-2);
    FN := FN*U+Double(3.45017939782574027900E-4);
    FN := FN*U+Double(4.63613749287867322088E-6);
    FN := FN*U+Double(3.05568983790257605827E-8);
    FN := FN*U+Double(1.02304514164907233465E-10);
    FN := FN*U+Double(1.72010743268161828879E-13);
    FN := FN*U+Double(1.34283276233062758925E-16);
    FN := FN*U+Double(3.76329711269987889006E-20);
    FD := Double(1.00000000000000000000E0);
    FD := FD*U+Double(7.51586398353378947175E-1);
    FD := FD*U+Double(1.16888925859191382142E-1);
    FD := FD*U+Double(6.44051526508858611005E-3);
    FD := FD*U+Double(1.55934409164153020873E-4);
    FD := FD*U+Double(1.84627567348930545870E-6);
    FD := FD*U+Double(1.12699224763999035261E-8);
    FD := FD*U+Double(3.60140029589371370404E-11);
    FD := FD*U+Double(5.88754533621578410010E-14);
    FD := FD*U+Double(4.52001434074129701496E-17);
    FD := FD*U+Double(1.25443237090011264384E-20);
    GN := Double(5.04442073643383265887E-1);
    GN := GN*U+Double(1.97102833525523411709E-1);
    GN := GN*U+Double(1.87648584092575249293E-2);
    GN := GN*U+Double(6.84079380915393090172E-4);
    GN := GN*U+Double(1.15138826111884280931E-5);
    GN := GN*U+Double(9.82852443688422223854E-8);
    GN := GN*U+Double(4.45344415861750144738E-10);
    GN := GN*U+Double(1.08268041139020870318E-12);
    GN := GN*U+Double(1.37555460633261799868E-15);
    GN := GN*U+Double(8.36354435630677421531E-19);
    GN := GN*U+Double(1.86958710162783235106E-22);
    GD := Double(1.00000000000000000000E0);
    GD := GD*U+Double(1.47495759925128324529E0);
    GD := GD*U+Double(3.37748989120019970451E-1);
    GD := GD*U+Double(2.53603741420338795122E-2);
    GD := GD*U+Double(8.14679107184306179049E-4);
    GD := GD*U+Double(1.27545075667729118702E-5);
    GD := GD*U+Double(1.04314589657571990585E-7);
    GD := GD*U+Double(4.60680728146520428211E-10);
    GD := GD*U+Double(1.10273215066240270757E-12);
    GD := GD*U+Double(1.38796531259578871258E-15);
    GD := GD*U+Double(8.39158816283118707363E-19);
    GD := GD*U+Double(1.86958710162783236342E-22);
    f := 1-u*FN/FD;
    g := t*GN/GD;
    t := MPIO2*x2;
    cc := cos(t);
    ss := sin(t);
    t := MPI*x;
    c := Double(0.5)+(f*ss-g*cc)/t;
    s := Double(0.5)-(f*cc+g*ss)/t;
    C := C*Sign(XXA);
    S := S*sign(XXA);
end;


end.