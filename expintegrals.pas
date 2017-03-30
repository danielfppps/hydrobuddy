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
unit expintegrals;
interface
uses Math, Sysutils, Ap;

function ExponentialIntegralEI(X : Double):Double;
function ExponentialIntegralEN(X : Double; N : AlglibInteger):Double;

implementation

(*************************************************************************
Exponential integral Ei(x)

              x
               -     t
              | |   e
   Ei(x) =   -|-   ---  dt .
            | |     t
             -
            -inf

Not defined for x <= 0.
See also expn.c.



ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,100       50000      8.6e-16     1.3e-16

Cephes Math Library Release 2.8:  May, 1999
Copyright 1999 by Stephen L. Moshier
*************************************************************************)
function ExponentialIntegralEI(X : Double):Double;
var
    EUL : Double;
    f : Double;
    f1 : Double;
    f2 : Double;
    w : Double;
begin
    EUL := Double(0.5772156649015328606065);
    if AP_FP_Less_Eq(X,0) then
    begin
        Result := 0;
        Exit;
    end;
    if AP_FP_Less(x,2) then
    begin
        f1 := -Double(5.350447357812542947283);
        f1 := f1*x+Double(218.5049168816613393830);
        f1 := f1*x-Double(4176.572384826693777058);
        f1 := f1*x+Double(55411.76756393557601232);
        f1 := f1*x-Double(331338.1331178144034309);
        f1 := f1*x+Double(1592627.163384945414220);
        f2 := Double(1.000000000000000000000);
        f2 := f2*x-Double(52.50547959112862969197);
        f2 := f2*x+Double(1259.616186786790571525);
        f2 := f2*x-Double(17565.49581973534652631);
        f2 := f2*x+Double(149306.2117002725991967);
        f2 := f2*x-Double(729494.9239640527645655);
        f2 := f2*x+Double(1592627.163384945429726);
        f := f1/f2;
        Result := EUL+Ln(x)+x*f;
        Exit;
    end;
    if AP_FP_Less(x,4) then
    begin
        w := 1/x;
        f1 := Double(1.981808503259689673238E-2);
        f1 := f1*w-Double(1.271645625984917501326);
        f1 := f1*w-Double(2.088160335681228318920);
        f1 := f1*w+Double(2.755544509187936721172);
        f1 := f1*w-Double(4.409507048701600257171E-1);
        f1 := f1*w+Double(4.665623805935891391017E-2);
        f1 := f1*w-Double(1.545042679673485262580E-3);
        f1 := f1*w+Double(7.059980605299617478514E-5);
        f2 := Double(1.000000000000000000000);
        f2 := f2*w+Double(1.476498670914921440652);
        f2 := f2*w+Double(5.629177174822436244827E-1);
        f2 := f2*w+Double(1.699017897879307263248E-1);
        f2 := f2*w+Double(2.291647179034212017463E-2);
        f2 := f2*w+Double(4.450150439728752875043E-3);
        f2 := f2*w+Double(1.727439612206521482874E-4);
        f2 := f2*w+Double(3.953167195549672482304E-5);
        f := f1/f2;
        Result := Exp(x)*w*(1+w*f);
        Exit;
    end;
    if AP_FP_Less(x,8) then
    begin
        w := 1/x;
        f1 := -Double(1.373215375871208729803);
        f1 := f1*w-Double(7.084559133740838761406E-1);
        f1 := f1*w+Double(1.580806855547941010501);
        f1 := f1*w-Double(2.601500427425622944234E-1);
        f1 := f1*w+Double(2.994674694113713763365E-2);
        f1 := f1*w-Double(1.038086040188744005513E-3);
        f1 := f1*w+Double(4.371064420753005429514E-5);
        f1 := f1*w+Double(2.141783679522602903795E-6);
        f2 := Double(1.000000000000000000000);
        f2 := f2*w+Double(8.585231423622028380768E-1);
        f2 := f2*w+Double(4.483285822873995129957E-1);
        f2 := f2*w+Double(7.687932158124475434091E-2);
        f2 := f2*w+Double(2.449868241021887685904E-2);
        f2 := f2*w+Double(8.832165941927796567926E-4);
        f2 := f2*w+Double(4.590952299511353531215E-4);
        f2 := f2*w+-Double(4.729848351866523044863E-6);
        f2 := f2*w+Double(2.665195537390710170105E-6);
        f := f1/f2;
        Result := exp(x)*w*(1+w*f);
        Exit;
    end;
    if AP_FP_Less(x,16) then
    begin
        w := 1/x;
        f1 := -Double(2.106934601691916512584);
        f1 := f1*w+Double(1.732733869664688041885);
        f1 := f1*w-Double(2.423619178935841904839E-1);
        f1 := f1*w+Double(2.322724180937565842585E-2);
        f1 := f1*w+Double(2.372880440493179832059E-4);
        f1 := f1*w-Double(8.343219561192552752335E-5);
        f1 := f1*w+Double(1.363408795605250394881E-5);
        f1 := f1*w-Double(3.655412321999253963714E-7);
        f1 := f1*w+Double(1.464941733975961318456E-8);
        f1 := f1*w+Double(6.176407863710360207074E-10);
        f2 := Double(1.000000000000000000000);
        f2 := f2*w-Double(2.298062239901678075778E-1);
        f2 := f2*w+Double(1.105077041474037862347E-1);
        f2 := f2*w-Double(1.566542966630792353556E-2);
        f2 := f2*w+Double(2.761106850817352773874E-3);
        f2 := f2*w-Double(2.089148012284048449115E-4);
        f2 := f2*w+Double(1.708528938807675304186E-5);
        f2 := f2*w-Double(4.459311796356686423199E-7);
        f2 := f2*w+Double(1.394634930353847498145E-8);
        f2 := f2*w+Double(6.150865933977338354138E-10);
        f := f1/f2;
        Result := Exp(x)*w*(1+w*f);
        Exit;
    end;
    if AP_FP_Less(x,32) then
    begin
        w := 1/x;
        f1 := -Double(2.458119367674020323359E-1);
        f1 := f1*w-Double(1.483382253322077687183E-1);
        f1 := f1*w+Double(7.248291795735551591813E-2);
        f1 := f1*w-Double(1.348315687380940523823E-2);
        f1 := f1*w+Double(1.342775069788636972294E-3);
        f1 := f1*w-Double(7.942465637159712264564E-5);
        f1 := f1*w+Double(2.644179518984235952241E-6);
        f1 := f1*w-Double(4.239473659313765177195E-8);
        f2 := Double(1.000000000000000000000);
        f2 := f2*w-Double(1.044225908443871106315E-1);
        f2 := f2*w-Double(2.676453128101402655055E-1);
        f2 := f2*w+Double(9.695000254621984627876E-2);
        f2 := f2*w-Double(1.601745692712991078208E-2);
        f2 := f2*w+Double(1.496414899205908021882E-3);
        f2 := f2*w-Double(8.462452563778485013756E-5);
        f2 := f2*w+Double(2.728938403476726394024E-6);
        f2 := f2*w-Double(4.239462431819542051337E-8);
        f := f1/f2;
        Result := Exp(x)*w*(1+w*f);
        Exit;
    end;
    if AP_FP_Less(x,64) then
    begin
        w := 1/x;
        f1 := Double(1.212561118105456670844E-1);
        f1 := f1*w-Double(5.823133179043894485122E-1);
        f1 := f1*w+Double(2.348887314557016779211E-1);
        f1 := f1*w-Double(3.040034318113248237280E-2);
        f1 := f1*w+Double(1.510082146865190661777E-3);
        f1 := f1*w-Double(2.523137095499571377122E-5);
        f2 := Double(1.000000000000000000000);
        f2 := f2*w-Double(1.002252150365854016662);
        f2 := f2*w+Double(2.928709694872224144953E-1);
        f2 := f2*w-Double(3.337004338674007801307E-2);
        f2 := f2*w+Double(1.560544881127388842819E-3);
        f2 := f2*w-Double(2.523137093603234562648E-5);
        f := f1/f2;
        Result := Exp(x)*w*(1+w*f);
        Exit;
    end;
    w := 1/x;
    f1 := -Double(7.657847078286127362028E-1);
    f1 := f1*w+Double(6.886192415566705051750E-1);
    f1 := f1*w-Double(2.132598113545206124553E-1);
    f1 := f1*w+Double(3.346107552384193813594E-2);
    f1 := f1*w-Double(3.076541477344756050249E-3);
    f1 := f1*w+Double(1.747119316454907477380E-4);
    f1 := f1*w-Double(6.103711682274170530369E-6);
    f1 := f1*w+Double(1.218032765428652199087E-7);
    f1 := f1*w-Double(1.086076102793290233007E-9);
    f2 := Double(1.000000000000000000000);
    f2 := f2*w-Double(1.888802868662308731041);
    f2 := f2*w+Double(1.066691687211408896850);
    f2 := f2*w-Double(2.751915982306380647738E-1);
    f2 := f2*w+Double(3.930852688233823569726E-2);
    f2 := f2*w-Double(3.414684558602365085394E-3);
    f2 := f2*w+Double(1.866844370703555398195E-4);
    f2 := f2*w-Double(6.345146083130515357861E-6);
    f2 := f2*w+Double(1.239754287483206878024E-7);
    f2 := f2*w-Double(1.086076102793126632978E-9);
    f := f1/f2;
    Result := exp(x)*w*(1+w*f);
end;


(*************************************************************************
Exponential integral En(x)

Evaluates the exponential integral

                inf.
                  -
                 | |   -xt
                 |    e
     E (x)  =    |    ----  dt.
      n          |      n
               | |     t
                -
                 1


Both n and x must be nonnegative.

The routine employs either a power series, a continued
fraction, or an asymptotic formula depending on the
relative values of n and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       10000       1.7e-15     3.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 2000 by Stephen L. Moshier
*************************************************************************)
function ExponentialIntegralEN(X : Double; N : AlglibInteger):Double;
var
    r : Double;
    t : Double;
    yk : Double;
    xk : Double;
    pk : Double;
    pkm1 : Double;
    pkm2 : Double;
    qk : Double;
    qkm1 : Double;
    qkm2 : Double;
    psi : Double;
    z : Double;
    i : AlglibInteger;
    k : AlglibInteger;
    big : Double;
    EUL : Double;
begin
    EUL := Double(0.57721566490153286060);
    big := Double(1.44115188075855872)*Power(10, 17);
    if (n<0) or AP_FP_Less(x,0) or AP_FP_Greater(x,170) or AP_FP_Eq(x,0) and (n<2) then
    begin
        Result := -1;
        Exit;
    end;
    if AP_FP_Eq(x,0) then
    begin
        Result := AP_Double(1)/(n-1);
        Exit;
    end;
    if n=0 then
    begin
        Result := Exp(-x)/x;
        Exit;
    end;
    if n>5000 then
    begin
        xk := x+n;
        yk := 1/(xk*xk);
        t := n;
        Result := yk*t*(6*x*x-8*t*x+t*t);
        Result := yk*(Result+t*(t-Double(2.0)*x));
        Result := yk*(Result+t);
        Result := (Result+1)*exp(-x)/xk;
        Exit;
    end;
    if AP_FP_Less_Eq(x,1) then
    begin
        psi := -EUL-Ln(x);
        i:=1;
        while i<=n-1 do
        begin
            psi := psi+AP_Double(1)/i;
            Inc(i);
        end;
        z := -x;
        xk := 0;
        yk := 1;
        pk := 1-n;
        if n=1 then
        begin
            Result := Double(0.0);
        end
        else
        begin
            Result := Double(1.0)/pk;
        end;
        repeat
            xk := xk+1;
            yk := yk*z/xk;
            pk := pk+1;
            if AP_FP_Neq(pk,0) then
            begin
                Result := Result+yk/pk;
            end;
            if AP_FP_Neq(Result,0) then
            begin
                t := AbsReal(yk/Result);
            end
            else
            begin
                t := 1;
            end;
        until AP_FP_Less(t,MachineEpsilon);
        t := 1;
        I:=1;
        while I<=N-1 do
        begin
            t := t*z/I;
            Inc(I);
        end;
        Result := psi*t-Result;
        Exit;
    end
    else
    begin
        k := 1;
        pkm2 := 1;
        qkm2 := x;
        pkm1 := Double(1.0);
        qkm1 := x+n;
        Result := pkm1/qkm1;
        repeat
            k := k+1;
            if k mod 2=1 then
            begin
                yk := 1;
                xk := n+AP_Double((k-1))/2;
            end
            else
            begin
                yk := x;
                xk := AP_Double(k)/2;
            end;
            pk := pkm1*yk+pkm2*xk;
            qk := qkm1*yk+qkm2*xk;
            if AP_FP_Neq(qk,0) then
            begin
                r := pk/qk;
                t := AbsReal((Result-r)/r);
                Result := r;
            end
            else
            begin
                t := 1;
            end;
            pkm2 := pkm1;
            pkm1 := pk;
            qkm2 := qkm1;
            qkm1 := qk;
            if AP_FP_Greater(AbsReal(pk),big) then
            begin
                pkm2 := pkm2/big;
                pkm1 := pkm1/big;
                qkm2 := qkm2/big;
                qkm1 := qkm1/big;
            end;
        until AP_FP_Less(t,MachineEpsilon);
        Result := Result*exp(-x);
    end;
end;


end.