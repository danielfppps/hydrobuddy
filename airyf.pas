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
unit airyf;
interface
uses Math, Sysutils, Ap;

procedure Airy(x : Double;
     var Ai : Double;
     var Aip : Double;
     var Bi : Double;
     var Bip : Double);

implementation

(*************************************************************************
Airy function

Solution of the differential equation

y"(x) = xy.

The function returns the two independent solutions Ai, Bi
and their first derivatives Ai'(x), Bi'(x).

Evaluation is by power series summation for small x,
by rational minimax approximations for large x.



ACCURACY:
Error criterion is absolute when function <= 1, relative
when function > 1, except * denotes relative error criterion.
For large negative x, the absolute error increases as x^1.5.
For large positive x, the relative error increases as x^1.5.

Arithmetic  domain   function  # trials      peak         rms
IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
procedure Airy(x : Double;
     var Ai : Double;
     var Aip : Double;
     var Bi : Double;
     var Bip : Double);
var
    z : Double;
    zz : Double;
    t : Double;
    f : Double;
    g : Double;
    uf : Double;
    ug : Double;
    k : Double;
    zeta : Double;
    theta : Double;
    domflg : AlglibInteger;
    c1 : Double;
    c2 : Double;
    sqrt3 : Double;
    sqpii : Double;
    AFN : Double;
    AFD : Double;
    AGN : Double;
    AGD : Double;
    APFN : Double;
    APFD : Double;
    APGN : Double;
    APGD : Double;
    AN : Double;
    AD : Double;
    APN : Double;
    APD : Double;
    BN16 : Double;
    BD16 : Double;
    BPPN : Double;
    BPPD : Double;
begin
    sqpii := Double(5.64189583547756286948E-1);
    c1 := Double(0.35502805388781723926);
    c2 := Double(0.258819403792806798405);
    sqrt3 := Double(1.732050807568877293527);
    domflg := 0;
    if AP_FP_Greater(x,Double(25.77)) then
    begin
        ai := 0;
        aip := 0;
        bi := MaxRealNumber;
        bip := MaxRealNumber;
        Exit;
    end;
    if AP_FP_Less(x,-Double(2.09)) then
    begin
        domflg := 15;
        t := sqrt(-x);
        zeta := -Double(2.0)*x*t/Double(3.0);
        t := sqrt(t);
        k := sqpii/t;
        z := Double(1.0)/zeta;
        zz := z*z;
        AFN := -Double(1.31696323418331795333E-1);
        AFN := AFN*zz-Double(6.26456544431912369773E-1);
        AFN := AFN*zz-Double(6.93158036036933542233E-1);
        AFN := AFN*zz-Double(2.79779981545119124951E-1);
        AFN := AFN*zz-Double(4.91900132609500318020E-2);
        AFN := AFN*zz-Double(4.06265923594885404393E-3);
        AFN := AFN*zz-Double(1.59276496239262096340E-4);
        AFN := AFN*zz-Double(2.77649108155232920844E-6);
        AFN := AFN*zz-Double(1.67787698489114633780E-8);
        AFD := Double(1.00000000000000000000E0);
        AFD := AFD*zz+Double(1.33560420706553243746E1);
        AFD := AFD*zz+Double(3.26825032795224613948E1);
        AFD := AFD*zz+Double(2.67367040941499554804E1);
        AFD := AFD*zz+Double(9.18707402907259625840E0);
        AFD := AFD*zz+Double(1.47529146771666414581E0);
        AFD := AFD*zz+Double(1.15687173795188044134E-1);
        AFD := AFD*zz+Double(4.40291641615211203805E-3);
        AFD := AFD*zz+Double(7.54720348287414296618E-5);
        AFD := AFD*zz+Double(4.51850092970580378464E-7);
        uf := Double(1.0)+zz*AFN/AFD;
        AGN := Double(1.97339932091685679179E-2);
        AGN := AGN*zz+Double(3.91103029615688277255E-1);
        AGN := AGN*zz+Double(1.06579897599595591108E0);
        AGN := AGN*zz+Double(9.39169229816650230044E-1);
        AGN := AGN*zz+Double(3.51465656105547619242E-1);
        AGN := AGN*zz+Double(6.33888919628925490927E-2);
        AGN := AGN*zz+Double(5.85804113048388458567E-3);
        AGN := AGN*zz+Double(2.82851600836737019778E-4);
        AGN := AGN*zz+Double(6.98793669997260967291E-6);
        AGN := AGN*zz+Double(8.11789239554389293311E-8);
        AGN := AGN*zz+Double(3.41551784765923618484E-10);
        AGD := Double(1.00000000000000000000E0);
        AGD := AGD*zz+Double(9.30892908077441974853E0);
        AGD := AGD*zz+Double(1.98352928718312140417E1);
        AGD := AGD*zz+Double(1.55646628932864612953E1);
        AGD := AGD*zz+Double(5.47686069422975497931E0);
        AGD := AGD*zz+Double(9.54293611618961883998E-1);
        AGD := AGD*zz+Double(8.64580826352392193095E-2);
        AGD := AGD*zz+Double(4.12656523824222607191E-3);
        AGD := AGD*zz+Double(1.01259085116509135510E-4);
        AGD := AGD*zz+Double(1.17166733214413521882E-6);
        AGD := AGD*zz+Double(4.91834570062930015649E-9);
        ug := z*AGN/AGD;
        theta := zeta+Double(0.25)*PI;
        f := sin(theta);
        g := cos(theta);
        ai := k*(f*uf-g*ug);
        bi := k*(g*uf+f*ug);
        APFN := Double(1.85365624022535566142E-1);
        APFN := APFN*zz+Double(8.86712188052584095637E-1);
        APFN := APFN*zz+Double(9.87391981747398547272E-1);
        APFN := APFN*zz+Double(4.01241082318003734092E-1);
        APFN := APFN*zz+Double(7.10304926289631174579E-2);
        APFN := APFN*zz+Double(5.90618657995661810071E-3);
        APFN := APFN*zz+Double(2.33051409401776799569E-4);
        APFN := APFN*zz+Double(4.08718778289035454598E-6);
        APFN := APFN*zz+Double(2.48379932900442457853E-8);
        APFD := Double(1.00000000000000000000E0);
        APFD := APFD*zz+Double(1.47345854687502542552E1);
        APFD := APFD*zz+Double(3.75423933435489594466E1);
        APFD := APFD*zz+Double(3.14657751203046424330E1);
        APFD := APFD*zz+Double(1.09969125207298778536E1);
        APFD := APFD*zz+Double(1.78885054766999417817E0);
        APFD := APFD*zz+Double(1.41733275753662636873E-1);
        APFD := APFD*zz+Double(5.44066067017226003627E-3);
        APFD := APFD*zz+Double(9.39421290654511171663E-5);
        APFD := APFD*zz+Double(5.65978713036027009243E-7);
        uf := Double(1.0)+zz*APFN/APFD;
        APGN := -Double(3.55615429033082288335E-2);
        APGN := APGN*zz-Double(6.37311518129435504426E-1);
        APGN := APGN*zz-Double(1.70856738884312371053E0);
        APGN := APGN*zz-Double(1.50221872117316635393E0);
        APGN := APGN*zz-Double(5.63606665822102676611E-1);
        APGN := APGN*zz-Double(1.02101031120216891789E-1);
        APGN := APGN*zz-Double(9.48396695961445269093E-3);
        APGN := APGN*zz-Double(4.60325307486780994357E-4);
        APGN := APGN*zz-Double(1.14300836484517375919E-5);
        APGN := APGN*zz-Double(1.33415518685547420648E-7);
        APGN := APGN*zz-Double(5.63803833958893494476E-10);
        APGD := Double(1.00000000000000000000E0);
        APGD := APGD*zz+Double(9.85865801696130355144E0);
        APGD := APGD*zz+Double(2.16401867356585941885E1);
        APGD := APGD*zz+Double(1.73130776389749389525E1);
        APGD := APGD*zz+Double(6.17872175280828766327E0);
        APGD := APGD*zz+Double(1.08848694396321495475E0);
        APGD := APGD*zz+Double(9.95005543440888479402E-2);
        APGD := APGD*zz+Double(4.78468199683886610842E-3);
        APGD := APGD*zz+Double(1.18159633322838625562E-4);
        APGD := APGD*zz+Double(1.37480673554219441465E-6);
        APGD := APGD*zz+Double(5.79912514929147598821E-9);
        ug := z*APGN/APGD;
        k := sqpii*t;
        aip := -k*(g*uf+f*ug);
        bip := k*(f*uf-g*ug);
        Exit;
    end;
    if AP_FP_Greater_Eq(x,Double(2.09)) then
    begin
        domflg := 5;
        t := sqrt(x);
        zeta := Double(2.0)*x*t/Double(3.0);
        g := exp(zeta);
        t := sqrt(t);
        k := Double(2.0)*t*g;
        z := Double(1.0)/zeta;
        AN := Double(3.46538101525629032477E-1);
        AN := AN*z+Double(1.20075952739645805542E1);
        AN := AN*z+Double(7.62796053615234516538E1);
        AN := AN*z+Double(1.68089224934630576269E2);
        AN := AN*z+Double(1.59756391350164413639E2);
        AN := AN*z+Double(7.05360906840444183113E1);
        AN := AN*z+Double(1.40264691163389668864E1);
        AN := AN*z+Double(9.99999999999999995305E-1);
        AD := Double(5.67594532638770212846E-1);
        AD := AD*z+Double(1.47562562584847203173E1);
        AD := AD*z+Double(8.45138970141474626562E1);
        AD := AD*z+Double(1.77318088145400459522E2);
        AD := AD*z+Double(1.64234692871529701831E2);
        AD := AD*z+Double(7.14778400825575695274E1);
        AD := AD*z+Double(1.40959135607834029598E1);
        AD := AD*z+Double(1.00000000000000000470E0);
        f := AN/AD;
        ai := sqpii*f/k;
        k := -Double(0.5)*sqpii*t/g;
        APN := Double(6.13759184814035759225E-1);
        APN := APN*z+Double(1.47454670787755323881E1);
        APN := APN*z+Double(8.20584123476060982430E1);
        APN := APN*z+Double(1.71184781360976385540E2);
        APN := APN*z+Double(1.59317847137141783523E2);
        APN := APN*z+Double(6.99778599330103016170E1);
        APN := APN*z+Double(1.39470856980481566958E1);
        APN := APN*z+Double(1.00000000000000000550E0);
        APD := Double(3.34203677749736953049E-1);
        APD := APD*z+Double(1.11810297306158156705E1);
        APD := APD*z+Double(7.11727352147859965283E1);
        APD := APD*z+Double(1.58778084372838313640E2);
        APD := APD*z+Double(1.53206427475809220834E2);
        APD := APD*z+Double(6.86752304592780337944E1);
        APD := APD*z+Double(1.38498634758259442477E1);
        APD := APD*z+Double(9.99999999999999994502E-1);
        f := APN/APD;
        aip := f*k;
        if AP_FP_Greater(x,Double(8.3203353)) then
        begin
            BN16 := -Double(2.53240795869364152689E-1);
            BN16 := BN16*z+Double(5.75285167332467384228E-1);
            BN16 := BN16*z-Double(3.29907036873225371650E-1);
            BN16 := BN16*z+Double(6.44404068948199951727E-2);
            BN16 := BN16*z-Double(3.82519546641336734394E-3);
            BD16 := Double(1.00000000000000000000E0);
            BD16 := BD16*z-Double(7.15685095054035237902E0);
            BD16 := BD16*z+Double(1.06039580715664694291E1);
            BD16 := BD16*z-Double(5.23246636471251500874E0);
            BD16 := BD16*z+Double(9.57395864378383833152E-1);
            BD16 := BD16*z-Double(5.50828147163549611107E-2);
            f := z*BN16/BD16;
            k := sqpii*g;
            bi := k*(Double(1.0)+f)/t;
            BPPN := Double(4.65461162774651610328E-1);
            BPPN := BPPN*z-Double(1.08992173800493920734E0);
            BPPN := BPPN*z+Double(6.38800117371827987759E-1);
            BPPN := BPPN*z-Double(1.26844349553102907034E-1);
            BPPN := BPPN*z+Double(7.62487844342109852105E-3);
            BPPD := Double(1.00000000000000000000E0);
            BPPD := BPPD*z-Double(8.70622787633159124240E0);
            BPPD := BPPD*z+Double(1.38993162704553213172E1);
            BPPD := BPPD*z-Double(7.14116144616431159572E0);
            BPPD := BPPD*z+Double(1.34008595960680518666E0);
            BPPD := BPPD*z-Double(7.84273211323341930448E-2);
            f := z*BPPN/BPPD;
            bip := k*t*(Double(1.0)+f);
            Exit;
        end;
    end;
    f := Double(1.0);
    g := x;
    t := Double(1.0);
    uf := Double(1.0);
    ug := x;
    k := Double(1.0);
    z := x*x*x;
    while AP_FP_Greater(t,MachineEpsilon) do
    begin
        uf := uf*z;
        k := k+Double(1.0);
        uf := uf/k;
        ug := ug*z;
        k := k+Double(1.0);
        ug := ug/k;
        uf := uf/k;
        f := f+uf;
        k := k+Double(1.0);
        ug := ug/k;
        g := g+ug;
        t := AbsReal(uf/f);
    end;
    uf := c1*f;
    ug := c2*g;
    if domflg mod 2=0 then
    begin
        ai := uf-ug;
    end;
    if domflg div 2 mod 2=0 then
    begin
        bi := sqrt3*(uf+ug);
    end;
    k := Double(4.0);
    uf := x*x/Double(2.0);
    ug := z/Double(3.0);
    f := uf;
    g := Double(1.0)+ug;
    uf := uf/Double(3.0);
    t := Double(1.0);
    while AP_FP_Greater(t,MachineEpsilon) do
    begin
        uf := uf*z;
        ug := ug/k;
        k := k+Double(1.0);
        ug := ug*z;
        uf := uf/k;
        f := f+uf;
        k := k+Double(1.0);
        ug := ug/k;
        uf := uf/k;
        g := g+ug;
        k := k+Double(1.0);
        t := AbsReal(ug/g);
    end;
    uf := c1*f;
    ug := c2*g;
    if domflg div 4 mod 2=0 then
    begin
        aip := uf-ug;
    end;
    if domflg div 8 mod 2=0 then
    begin
        bip := sqrt3*(uf+ug);
    end;
end;


end.