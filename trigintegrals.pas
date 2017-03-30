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
unit trigintegrals;
interface
uses Math, Sysutils, Ap;

procedure SineCosineIntegrals(X : Double; var SI : Double; var CI : Double);
procedure HyperbolicSineCosineIntegrals(X : Double;
     var Shi : Double;
     var Chi : Double);

implementation

procedure ChebIterationShiChi(x : Double;
     c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);forward;


(*************************************************************************
Sine and cosine integrals

Evaluates the integrals

                         x
                         -
                        |  cos t - 1
  Ci(x) = eul + ln x +  |  --------- dt,
                        |      t
                       -
                        0
            x
            -
           |  sin t
  Si(x) =  |  ----- dt
           |    t
          -
           0

where eul = 0.57721566490153286061 is Euler's constant.
The integrals are approximated by rational functions.
For x > 8 auxiliary functions f(x) and g(x) are employed
such that

Ci(x) = f(x) sin(x) - g(x) cos(x)
Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)


ACCURACY:
   Test interval = [0,50].
Absolute error, except relative when > 1:
arithmetic   function   # trials      peak         rms
   IEEE        Si        30000       4.4e-16     7.3e-17
   IEEE        Ci        30000       6.9e-16     5.1e-17

Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
*************************************************************************)
procedure SineCosineIntegrals(X : Double; var SI : Double; var CI : Double);
var
    Z : Double;
    C : Double;
    S : Double;
    F : Double;
    G : Double;
    Sg : AlglibInteger;
    SN : Double;
    SD : Double;
    CN : Double;
    CD : Double;
    FN : Double;
    FD : Double;
    GN : Double;
    GD : Double;
begin
    if AP_FP_Less(x,0) then
    begin
        Sg := -1;
        x := -x;
    end
    else
    begin
        Sg := 0;
    end;
    if AP_FP_Eq(x,0) then
    begin
        Si := 0;
        Ci := -MaxRealNumber;
        Exit;
    end;
    if AP_FP_Greater(x,Double(1.0E9)) then
    begin
        Si := Double(1.570796326794896619)-cos(x)/x;
        Ci := sin(x)/x;
        Exit;
    end;
    if AP_FP_Less_Eq(X,4) then
    begin
        z := x*x;
        SN := -Double(8.39167827910303881427E-11);
        SN := SN*z+Double(4.62591714427012837309E-8);
        SN := SN*z-Double(9.75759303843632795789E-6);
        SN := SN*z+Double(9.76945438170435310816E-4);
        SN := SN*z-Double(4.13470316229406538752E-2);
        SN := SN*z+Double(1.00000000000000000302E0);
        SD := Double(2.03269266195951942049E-12);
        SD := SD*z+Double(1.27997891179943299903E-9);
        SD := SD*z+Double(4.41827842801218905784E-7);
        SD := SD*z+Double(9.96412122043875552487E-5);
        SD := SD*z+Double(1.42085239326149893930E-2);
        SD := SD*z+Double(9.99999999999999996984E-1);
        s := x*SN/SD;
        CN := Double(2.02524002389102268789E-11);
        CN := CN*z-Double(1.35249504915790756375E-8);
        CN := CN*z+Double(3.59325051419993077021E-6);
        CN := CN*z-Double(4.74007206873407909465E-4);
        CN := CN*z+Double(2.89159652607555242092E-2);
        CN := CN*z-Double(1.00000000000000000080E0);
        CD := Double(4.07746040061880559506E-12);
        CD := CD*z+Double(3.06780997581887812692E-9);
        CD := CD*z+Double(1.23210355685883423679E-6);
        CD := CD*z+Double(3.17442024775032769882E-4);
        CD := CD*z+Double(5.10028056236446052392E-2);
        CD := CD*z+Double(4.00000000000000000080E0);
        c := z*CN/CD;
        if Sg<>0 then
        begin
            s := -s;
        end;
        si := s;
        ci := Double(0.57721566490153286061)+Ln(x)+c;
        Exit;
    end;
    s := sin(x);
    c := cos(x);
    z := Double(1.0)/(x*x);
    if AP_FP_Less(x,8) then
    begin
        FN := Double(4.23612862892216586994E0);
        FN := FN*z+Double(5.45937717161812843388E0);
        FN := FN*z+Double(1.62083287701538329132E0);
        FN := FN*z+Double(1.67006611831323023771E-1);
        FN := FN*z+Double(6.81020132472518137426E-3);
        FN := FN*z+Double(1.08936580650328664411E-4);
        FN := FN*z+Double(5.48900223421373614008E-7);
        FD := Double(1.00000000000000000000E0);
        FD := FD*z+Double(8.16496634205391016773E0);
        FD := FD*z+Double(7.30828822505564552187E0);
        FD := FD*z+Double(1.86792257950184183883E0);
        FD := FD*z+Double(1.78792052963149907262E-1);
        FD := FD*z+Double(7.01710668322789753610E-3);
        FD := FD*z+Double(1.10034357153915731354E-4);
        FD := FD*z+Double(5.48900252756255700982E-7);
        f := FN/(x*FD);
        GN := Double(8.71001698973114191777E-2);
        GN := GN*z+Double(6.11379109952219284151E-1);
        GN := GN*z+Double(3.97180296392337498885E-1);
        GN := GN*z+Double(7.48527737628469092119E-2);
        GN := GN*z+Double(5.38868681462177273157E-3);
        GN := GN*z+Double(1.61999794598934024525E-4);
        GN := GN*z+Double(1.97963874140963632189E-6);
        GN := GN*z+Double(7.82579040744090311069E-9);
        GD := Double(1.00000000000000000000E0);
        GD := GD*z+Double(1.64402202413355338886E0);
        GD := GD*z+Double(6.66296701268987968381E-1);
        GD := GD*z+Double(9.88771761277688796203E-2);
        GD := GD*z+Double(6.22396345441768420760E-3);
        GD := GD*z+Double(1.73221081474177119497E-4);
        GD := GD*z+Double(2.02659182086343991969E-6);
        GD := GD*z+Double(7.82579218933534490868E-9);
        g := z*GN/GD;
    end
    else
    begin
        FN := Double(4.55880873470465315206E-1);
        FN := FN*z+Double(7.13715274100146711374E-1);
        FN := FN*z+Double(1.60300158222319456320E-1);
        FN := FN*z+Double(1.16064229408124407915E-2);
        FN := FN*z+Double(3.49556442447859055605E-4);
        FN := FN*z+Double(4.86215430826454749482E-6);
        FN := FN*z+Double(3.20092790091004902806E-8);
        FN := FN*z+Double(9.41779576128512936592E-11);
        FN := FN*z+Double(9.70507110881952024631E-14);
        FD := Double(1.00000000000000000000E0);
        FD := FD*z+Double(9.17463611873684053703E-1);
        FD := FD*z+Double(1.78685545332074536321E-1);
        FD := FD*z+Double(1.22253594771971293032E-2);
        FD := FD*z+Double(3.58696481881851580297E-4);
        FD := FD*z+Double(4.92435064317881464393E-6);
        FD := FD*z+Double(3.21956939101046018377E-8);
        FD := FD*z+Double(9.43720590350276732376E-11);
        FD := FD*z+Double(9.70507110881952025725E-14);
        f := FN/(x*FD);
        GN := Double(6.97359953443276214934E-1);
        GN := GN*z+Double(3.30410979305632063225E-1);
        GN := GN*z+Double(3.84878767649974295920E-2);
        GN := GN*z+Double(1.71718239052347903558E-3);
        GN := GN*z+Double(3.48941165502279436777E-5);
        GN := GN*z+Double(3.47131167084116673800E-7);
        GN := GN*z+Double(1.70404452782044526189E-9);
        GN := GN*z+Double(3.85945925430276600453E-12);
        GN := GN*z+Double(3.14040098946363334640E-15);
        GD := Double(1.00000000000000000000E0);
        GD := GD*z+Double(1.68548898811011640017E0);
        GD := GD*z+Double(4.87852258695304967486E-1);
        GD := GD*z+Double(4.67913194259625806320E-2);
        GD := GD*z+Double(1.90284426674399523638E-3);
        GD := GD*z+Double(3.68475504442561108162E-5);
        GD := GD*z+Double(3.57043223443740838771E-7);
        GD := GD*z+Double(1.72693748966316146736E-9);
        GD := GD*z+Double(3.87830166023954706752E-12);
        GD := GD*z+Double(3.14040098946363335242E-15);
        g := z*GN/GD;
    end;
    si := Double(1.570796326794896619)-f*c-g*s;
    if sg<>0 then
    begin
        si := -si;
    end;
    ci := f*s-g*c;
end;


(*************************************************************************
Hyperbolic sine and cosine integrals

Approximates the integrals

                           x
                           -
                          | |   cosh t - 1
  Chi(x) = eul + ln x +   |    -----------  dt,
                        | |          t
                         -
                         0

              x
              -
             | |  sinh t
  Shi(x) =   |    ------  dt
           | |       t
            -
            0

where eul = 0.57721566490153286061 is Euler's constant.
The integrals are evaluated by power series for x < 8
and by Chebyshev expansions for x between 8 and 88.
For large x, both functions approach exp(x)/2x.
Arguments greater than 88 in magnitude return MAXNUM.


ACCURACY:

Test interval 0 to 88.
                     Relative error:
arithmetic   function  # trials      peak         rms
   IEEE         Shi      30000       6.9e-16     1.6e-16
       Absolute error, except relative when |Chi| > 1:
   IEEE         Chi      30000       8.4e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
procedure HyperbolicSineCosineIntegrals(X : Double;
     var Shi : Double;
     var Chi : Double);
var
    k : Double;
    z : Double;
    c : Double;
    s : Double;
    a : Double;
    sg : AlglibInteger;
    b0 : Double;
    b1 : Double;
    b2 : Double;
begin
    if AP_FP_Less(x,0) then
    begin
        sg := -1;
        x := -x;
    end
    else
    begin
        sg := 0;
    end;
    if AP_FP_Eq(x,0) then
    begin
        Shi := 0;
        Chi := -MaxRealNumber;
        Exit;
    end;
    if AP_FP_Less(x,Double(8.0)) then
    begin
        z := x*x;
        a := Double(1.0);
        s := Double(1.0);
        c := Double(0.0);
        k := Double(2.0);
        repeat
            a := a*z/k;
            c := c+a/k;
            k := k+Double(1.0);
            a := a/k;
            s := s+a/k;
            k := k+Double(1.0);
        until AP_FP_Less(AbsReal(a/s),MachineEpsilon);
        s := s*x;
    end
    else
    begin
        if AP_FP_Less(x,Double(18.0)) then
        begin
            a := (Double(576.0)/x-Double(52.0))/Double(10.0);
            k := exp(x)/x;
            b0 := Double(1.83889230173399459482E-17);
            b1 := Double(0.0);
            ChebIterationShiChi(a, -Double(9.55485532279655569575E-17), b0, b1, b2);
            ChebIterationShiChi(a, Double(2.04326105980879882648E-16), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.09896949074905343022E-15), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.31313534344092599234E-14), b0, b1, b2);
            ChebIterationShiChi(a, Double(5.93976226264314278932E-14), b0, b1, b2);
            ChebIterationShiChi(a, -Double(3.47197010497749154755E-14), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.40059764613117131000E-12), b0, b1, b2);
            ChebIterationShiChi(a, Double(9.49044626224223543299E-12), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.61596181145435454033E-11), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.77899784436430310321E-10), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.35455469767246947469E-9), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.03257121792819495123E-9), b0, b1, b2);
            ChebIterationShiChi(a, -Double(3.56699611114982536845E-8), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.44818877384267342057E-7), b0, b1, b2);
            ChebIterationShiChi(a, Double(7.82018215184051295296E-7), b0, b1, b2);
            ChebIterationShiChi(a, -Double(5.39919118403805073710E-6), b0, b1, b2);
            ChebIterationShiChi(a, -Double(3.12458202168959833422E-5), b0, b1, b2);
            ChebIterationShiChi(a, Double(8.90136741950727517826E-5), b0, b1, b2);
            ChebIterationShiChi(a, Double(2.02558474743846862168E-3), b0, b1, b2);
            ChebIterationShiChi(a, Double(2.96064440855633256972E-2), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.11847751047257036625E0), b0, b1, b2);
            s := k*Double(0.5)*(b0-b2);
            b0 := -Double(8.12435385225864036372E-18);
            b1 := Double(0.0);
            ChebIterationShiChi(a, Double(2.17586413290339214377E-17), b0, b1, b2);
            ChebIterationShiChi(a, Double(5.22624394924072204667E-17), b0, b1, b2);
            ChebIterationShiChi(a, -Double(9.48812110591690559363E-16), b0, b1, b2);
            ChebIterationShiChi(a, Double(5.35546311647465209166E-15), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.21009970113732918701E-14), b0, b1, b2);
            ChebIterationShiChi(a, -Double(6.00865178553447437951E-14), b0, b1, b2);
            ChebIterationShiChi(a, Double(7.16339649156028587775E-13), b0, b1, b2);
            ChebIterationShiChi(a, -Double(2.93496072607599856104E-12), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.40359438136491256904E-12), b0, b1, b2);
            ChebIterationShiChi(a, Double(8.76302288609054966081E-11), b0, b1, b2);
            ChebIterationShiChi(a, -Double(4.40092476213282340617E-10), b0, b1, b2);
            ChebIterationShiChi(a, -Double(1.87992075640569295479E-10), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.31458150989474594064E-8), b0, b1, b2);
            ChebIterationShiChi(a, -Double(4.75513930924765465590E-8), b0, b1, b2);
            ChebIterationShiChi(a, -Double(2.21775018801848880741E-7), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.94635531373272490962E-6), b0, b1, b2);
            ChebIterationShiChi(a, Double(4.33505889257316408893E-6), b0, b1, b2);
            ChebIterationShiChi(a, -Double(6.13387001076494349496E-5), b0, b1, b2);
            ChebIterationShiChi(a, -Double(3.13085477492997465138E-4), b0, b1, b2);
            ChebIterationShiChi(a, Double(4.97164789823116062801E-4), b0, b1, b2);
            ChebIterationShiChi(a, Double(2.64347496031374526641E-2), b0, b1, b2);
            ChebIterationShiChi(a, Double(1.11446150876699213025E0), b0, b1, b2);
            c := k*Double(0.5)*(b0-b2);
        end
        else
        begin
            if AP_FP_Less_Eq(x,Double(88.0)) then
            begin
                a := (Double(6336.0)/x-Double(212.0))/Double(70.0);
                k := exp(x)/x;
                b0 := -Double(1.05311574154850938805E-17);
                b1 := Double(0.0);
                ChebIterationShiChi(a, Double(2.62446095596355225821E-17), b0, b1, b2);
                ChebIterationShiChi(a, Double(8.82090135625368160657E-17), b0, b1, b2);
                ChebIterationShiChi(a, -Double(3.38459811878103047136E-16), b0, b1, b2);
                ChebIterationShiChi(a, -Double(8.30608026366935789136E-16), b0, b1, b2);
                ChebIterationShiChi(a, Double(3.93397875437050071776E-15), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.01765565969729044505E-14), b0, b1, b2);
                ChebIterationShiChi(a, -Double(4.21128170307640802703E-14), b0, b1, b2);
                ChebIterationShiChi(a, -Double(1.60818204519802480035E-13), b0, b1, b2);
                ChebIterationShiChi(a, Double(3.34714954175994481761E-13), b0, b1, b2);
                ChebIterationShiChi(a, Double(2.72600352129153073807E-12), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.66894954752839083608E-12), b0, b1, b2);
                ChebIterationShiChi(a, -Double(3.49278141024730899554E-11), b0, b1, b2);
                ChebIterationShiChi(a, -Double(1.58580661666482709598E-10), b0, b1, b2);
                ChebIterationShiChi(a, -Double(1.79289437183355633342E-10), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.76281629144264523277E-9), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.69050228879421288846E-8), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.25391771228487041649E-7), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.16229947068677338732E-6), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.61038260117376323993E-5), b0, b1, b2);
                ChebIterationShiChi(a, Double(3.49810375601053973070E-4), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.28478065259647610779E-2), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.03665722588798326712E0), b0, b1, b2);
                s := k*Double(0.5)*(b0-b2);
                b0 := Double(8.06913408255155572081E-18);
                b1 := Double(0.0);
                ChebIterationShiChi(a, -Double(2.08074168180148170312E-17), b0, b1, b2);
                ChebIterationShiChi(a, -Double(5.98111329658272336816E-17), b0, b1, b2);
                ChebIterationShiChi(a, Double(2.68533951085945765591E-16), b0, b1, b2);
                ChebIterationShiChi(a, Double(4.52313941698904694774E-16), b0, b1, b2);
                ChebIterationShiChi(a, -Double(3.10734917335299464535E-15), b0, b1, b2);
                ChebIterationShiChi(a, -Double(4.42823207332531972288E-15), b0, b1, b2);
                ChebIterationShiChi(a, Double(3.49639695410806959872E-14), b0, b1, b2);
                ChebIterationShiChi(a, Double(6.63406731718911586609E-14), b0, b1, b2);
                ChebIterationShiChi(a, -Double(3.71902448093119218395E-13), b0, b1, b2);
                ChebIterationShiChi(a, -Double(1.27135418132338309016E-12), b0, b1, b2);
                ChebIterationShiChi(a, Double(2.74851141935315395333E-12), b0, b1, b2);
                ChebIterationShiChi(a, Double(2.33781843985453438400E-11), b0, b1, b2);
                ChebIterationShiChi(a, Double(2.71436006377612442764E-11), b0, b1, b2);
                ChebIterationShiChi(a, -Double(2.56600180000355990529E-10), b0, b1, b2);
                ChebIterationShiChi(a, -Double(1.61021375163803438552E-9), b0, b1, b2);
                ChebIterationShiChi(a, -Double(4.72543064876271773512E-9), b0, b1, b2);
                ChebIterationShiChi(a, -Double(3.00095178028681682282E-9), b0, b1, b2);
                ChebIterationShiChi(a, Double(7.79387474390914922337E-8), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.06942765566401507066E-6), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.59503164802313196374E-5), b0, b1, b2);
                ChebIterationShiChi(a, Double(3.49592575153777996871E-4), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.28475387530065247392E-2), b0, b1, b2);
                ChebIterationShiChi(a, Double(1.03665693917934275131E0), b0, b1, b2);
                c := k*Double(0.5)*(b0-b2);
            end
            else
            begin
                if sg<>0 then
                begin
                    shi := -MaxRealNumber;
                end
                else
                begin
                    shi := MaxRealNumber;
                end;
                chi := MaxRealNumber;
                Exit;
            end;
        end;
    end;
    if sg<>0 then
    begin
        s := -s;
    end;
    shi := s;
    chi := Double(0.57721566490153286061)+Ln(x)+c;
end;


procedure ChebIterationShiChi(x : Double;
     c : Double;
     var b0 : Double;
     var b1 : Double;
     var b2 : Double);
begin
    b2 := b1;
    b1 := b0;
    b0 := x*b1-b2+c;
end;


end.