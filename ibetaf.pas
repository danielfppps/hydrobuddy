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
unit ibetaf;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr;

function IncompleteBeta(a : Double; b : Double; x : Double):Double;
function InvIncompleteBeta(a : Double; b : Double; y : Double):Double;

implementation

function IncompleteBetaFE(a : Double;
     b : Double;
     x : Double;
     big : Double;
     biginv : Double):Double;forward;
function IncompleteBetaFE2(a : Double;
     b : Double;
     x : Double;
     big : Double;
     biginv : Double):Double;forward;
function IncompleteBetaPS(a : Double;
     b : Double;
     x : Double;
     MAXGAM : Double):Double;forward;


(*************************************************************************
Incomplete beta integral

Returns incomplete beta integral of the arguments, evaluated
from zero to x.  The function is defined as

                 x
    -            -
   | (a+b)      | |  a-1     b-1
 -----------    |   t   (1-t)   dt.
  -     -     | |
 | (a) | (b)   -
                0

The domain of definition is 0 <= x <= 1.  In this
implementation a and b are restricted to positive values.
The integral from x to 1 may be obtained by the symmetry
relation

   1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).

The integral is evaluated by a continued fraction expansion
or, when b*x is small, by a power series.

ACCURACY:

Tested at uniformly distributed random points (a,b,x) with a and b
in "domain" and x between 0 and 1.
                                       Relative error
arithmetic   domain     # trials      peak         rms
   IEEE      0,5         10000       6.9e-15     4.5e-16
   IEEE      0,85       250000       2.2e-13     1.7e-14
   IEEE      0,1000      30000       5.3e-12     6.3e-13
   IEEE      0,10000    250000       9.3e-11     7.1e-12
   IEEE      0,100000    10000       8.7e-10     4.8e-11
Outputs smaller than the IEEE gradual underflow threshold
were excluded from these statistics.

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteBeta(a : Double; b : Double; x : Double):Double;
var
    t : Double;
    xc : Double;
    w : Double;
    y : Double;
    flag : AlglibInteger;
    sg : Double;
    big : Double;
    biginv : Double;
    MAXGAM : Double;
    MINLOG : Double;
    MAXLOG : Double;
begin
    big := Double(4.503599627370496e15);
    biginv := Double(2.22044604925031308085e-16);
    MAXGAM := Double(171.624376956302725);
    MINLOG := Ln(MinRealNumber);
    MAXLOG := Ln(MaxRealNumber);
    Assert(AP_FP_Greater(a,0) and AP_FP_Greater(b,0), 'Domain error in IncompleteBeta');
    Assert(AP_FP_Greater_Eq(x,0) and AP_FP_Less_Eq(x,1), 'Domain error in IncompleteBeta');
    if AP_FP_Eq(x,0) then
    begin
        Result := 0;
        Exit;
    end;
    if AP_FP_Eq(x,1) then
    begin
        Result := 1;
        Exit;
    end;
    flag := 0;
    if AP_FP_Less_Eq(b*x,Double(1.0)) and AP_FP_Less_Eq(x,Double(0.95)) then
    begin
        Result := IncompleteBetaPS(a, b, x, MAXGAM);
        Exit;
    end;
    w := Double(1.0)-x;
    if AP_FP_Greater(x,a/(a+b)) then
    begin
        flag := 1;
        t := a;
        a := b;
        b := t;
        xc := x;
        x := w;
    end
    else
    begin
        xc := w;
    end;
    if (flag=1) and AP_FP_Less_Eq(b*x,Double(1.0)) and AP_FP_Less_Eq(x,Double(0.95)) then
    begin
        t := IncompleteBetaPS(a, b, x, MAXGAM);
        if AP_FP_Less_Eq(t,MachineEpsilon) then
        begin
            Result := Double(1.0)-MachineEpsilon;
        end
        else
        begin
            Result := Double(1.0)-t;
        end;
        Exit;
    end;
    y := x*(a+b-Double(2.0))-(a-Double(1.0));
    if AP_FP_Less(y,Double(0.0)) then
    begin
        w := IncompleteBetaFE(a, b, x, big, biginv);
    end
    else
    begin
        w := IncompleteBetaFE2(a, b, x, big, biginv)/xc;
    end;
    y := a*ln(x);
    t := b*ln(xc);
    if AP_FP_Less(a+b,MAXGAM) and AP_FP_Less(absReal(y),MAXLOG) and AP_FP_Less(absReal(t),MAXLOG) then
    begin
        t := power(xc, b);
        t := t*power(x, a);
        t := t/a;
        t := t*w;
        t := t*(gamma(a+b)/(gamma(a)*gamma(b)));
        if flag=1 then
        begin
            if AP_FP_Less_Eq(t,MachineEpsilon) then
            begin
                Result := Double(1.0)-MachineEpsilon;
            end
            else
            begin
                Result := Double(1.0)-t;
            end;
        end
        else
        begin
            Result := t;
        end;
        Exit;
    end;
    y := y+t+lngamma(a+b, sg)-lngamma(a, sg)-lngamma(b, sg);
    y := y+ln(w/a);
    if AP_FP_Less(y,MINLOG) then
    begin
        t := Double(0.0);
    end
    else
    begin
        t := exp(y);
    end;
    if flag=1 then
    begin
        if AP_FP_Less_Eq(t,MachineEpsilon) then
        begin
            t := Double(1.0)-MachineEpsilon;
        end
        else
        begin
            t := Double(1.0)-t;
        end;
    end;
    Result := t;
end;


(*************************************************************************
Inverse of imcomplete beta integral

Given y, the function finds x such that

 incbet( a, b, x ) = y .

The routine performs interval halving or Newton iterations to find the
root of incbet(a,b,x) - y = 0.


ACCURACY:

                     Relative error:
               x     a,b
arithmetic   domain  domain  # trials    peak       rms
   IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
   IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
   IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
With a and b constrained to half-integer or integer values:
   IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
   IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
With a = .5, b constrained to half-integer or integer values:
   IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1996, 2000 by Stephen L. Moshier
*************************************************************************)
function InvIncompleteBeta(a : Double; b : Double; y : Double):Double;
var
    aaa : Double;
    bbb : Double;
    y0 : Double;
    d : Double;
    yyy : Double;
    x : Double;
    x0 : Double;
    x1 : Double;
    lgm : Double;
    yp : Double;
    di : Double;
    dithresh : Double;
    yl : Double;
    yh : Double;
    xt : Double;
    i : AlglibInteger;
    rflg : AlglibInteger;
    dir : AlglibInteger;
    nflg : AlglibInteger;
    s : Double;
    MainLoopPos : AlglibInteger;
    ihalve : AlglibInteger;
    ihalvecycle : AlglibInteger;
    newt : AlglibInteger;
    newtcycle : AlglibInteger;
    breaknewtcycle : AlglibInteger;
    breakihalvecycle : AlglibInteger;
begin
    i := 0;
    Assert(AP_FP_Greater_Eq(y,0) and AP_FP_Less_Eq(y,1), 'Domain error in InvIncompleteBeta');
    if AP_FP_Eq(y,0) then
    begin
        Result := 0;
        Exit;
    end;
    if AP_FP_Eq(y,Double(1.0)) then
    begin
        Result := 1;
        Exit;
    end;
    x0 := Double(0.0);
    yl := Double(0.0);
    x1 := Double(1.0);
    yh := Double(1.0);
    nflg := 0;
    MainLoopPos := 0;
    ihalve := 1;
    ihalvecycle := 2;
    newt := 3;
    newtcycle := 4;
    breaknewtcycle := 5;
    breakihalvecycle := 6;
    while True do
    begin
        
        //
        // start
        //
        if MainLoopPos=0 then
        begin
            if AP_FP_Less_Eq(a,Double(1.0)) or AP_FP_Less_Eq(b,Double(1.0)) then
            begin
                dithresh := Double(1.0e-6);
                rflg := 0;
                aaa := a;
                bbb := b;
                y0 := y;
                x := aaa/(aaa+bbb);
                yyy := IncompleteBeta(aaa, bbb, x);
                MainLoopPos := ihalve;
                Continue;
            end
            else
            begin
                dithresh := Double(1.0e-4);
            end;
            yp := -InvNormalDistribution(y);
            if AP_FP_Greater(y,Double(0.5)) then
            begin
                rflg := 1;
                aaa := b;
                bbb := a;
                y0 := Double(1.0)-y;
                yp := -yp;
            end
            else
            begin
                rflg := 0;
                aaa := a;
                bbb := b;
                y0 := y;
            end;
            lgm := (yp*yp-Double(3.0))/Double(6.0);
            x := Double(2.0)/(Double(1.0)/(Double(2.0)*aaa-Double(1.0))+Double(1.0)/(Double(2.0)*bbb-Double(1.0)));
            d := yp*sqrt(x+lgm)/x-(Double(1.0)/(Double(2.0)*bbb-Double(1.0))-Double(1.0)/(Double(2.0)*aaa-Double(1.0)))*(lgm+Double(5.0)/Double(6.0)-Double(2.0)/(Double(3.0)*x));
            d := Double(2.0)*d;
            if AP_FP_Less(d,Ln(MinRealNumber)) then
            begin
                x := 0;
                Break;
            end;
            x := aaa/(aaa+bbb*exp(d));
            yyy := IncompleteBeta(aaa, bbb, x);
            yp := (yyy-y0)/y0;
            if AP_FP_Less(absReal(yp),Double(0.2)) then
            begin
                MainLoopPos := newt;
                Continue;
            end;
            MainLoopPos := ihalve;
            Continue;
        end;
        
        //
        // ihalve
        //
        if MainLoopPos=ihalve then
        begin
            dir := 0;
            di := Double(0.5);
            i := 0;
            MainLoopPos := ihalvecycle;
            Continue;
        end;
        
        //
        // ihalvecycle
        //
        if MainLoopPos=ihalvecycle then
        begin
            if i<=99 then
            begin
                if i<>0 then
                begin
                    x := x0+di*(x1-x0);
                    if AP_FP_Eq(x,Double(1.0)) then
                    begin
                        x := Double(1.0)-MachineEpsilon;
                    end;
                    if AP_FP_Eq(x,Double(0.0)) then
                    begin
                        di := Double(0.5);
                        x := x0+di*(x1-x0);
                        if AP_FP_Eq(x,Double(0.0)) then
                        begin
                            Break;
                        end;
                    end;
                    yyy := IncompleteBeta(aaa, bbb, x);
                    yp := (x1-x0)/(x1+x0);
                    if AP_FP_Less(absReal(yp),dithresh) then
                    begin
                        MainLoopPos := newt;
                        Continue;
                    end;
                    yp := (yyy-y0)/y0;
                    if AP_FP_Less(absReal(yp),dithresh) then
                    begin
                        MainLoopPos := newt;
                        Continue;
                    end;
                end;
                if AP_FP_Less(yyy,y0) then
                begin
                    x0 := x;
                    yl := yyy;
                    if dir<0 then
                    begin
                        dir := 0;
                        di := Double(0.5);
                    end
                    else
                    begin
                        if dir>3 then
                        begin
                            di := Double(1.0)-(Double(1.0)-di)*(Double(1.0)-di);
                        end
                        else
                        begin
                            if dir>1 then
                            begin
                                di := Double(0.5)*di+Double(0.5);
                            end
                            else
                            begin
                                di := (y0-yyy)/(yh-yl);
                            end;
                        end;
                    end;
                    dir := dir+1;
                    if AP_FP_Greater(x0,Double(0.75)) then
                    begin
                        if rflg=1 then
                        begin
                            rflg := 0;
                            aaa := a;
                            bbb := b;
                            y0 := y;
                        end
                        else
                        begin
                            rflg := 1;
                            aaa := b;
                            bbb := a;
                            y0 := Double(1.0)-y;
                        end;
                        x := Double(1.0)-x;
                        yyy := IncompleteBeta(aaa, bbb, x);
                        x0 := Double(0.0);
                        yl := Double(0.0);
                        x1 := Double(1.0);
                        yh := Double(1.0);
                        MainLoopPos := ihalve;
                        Continue;
                    end;
                end
                else
                begin
                    x1 := x;
                    if (rflg=1) and AP_FP_Less(x1,MachineEpsilon) then
                    begin
                        x := Double(0.0);
                        Break;
                    end;
                    yh := yyy;
                    if dir>0 then
                    begin
                        dir := 0;
                        di := Double(0.5);
                    end
                    else
                    begin
                        if dir<-3 then
                        begin
                            di := di*di;
                        end
                        else
                        begin
                            if dir<-1 then
                            begin
                                di := Double(0.5)*di;
                            end
                            else
                            begin
                                di := (yyy-y0)/(yh-yl);
                            end;
                        end;
                    end;
                    dir := dir-1;
                end;
                i := i+1;
                MainLoopPos := ihalvecycle;
                Continue;
            end
            else
            begin
                MainLoopPos := breakihalvecycle;
                Continue;
            end;
        end;
        
        //
        // breakihalvecycle
        //
        if MainLoopPos=breakihalvecycle then
        begin
            if AP_FP_Greater_Eq(x0,Double(1.0)) then
            begin
                x := Double(1.0)-MachineEpsilon;
                Break;
            end;
            if AP_FP_Less_Eq(x,Double(0.0)) then
            begin
                x := Double(0.0);
                Break;
            end;
            MainLoopPos := newt;
            Continue;
        end;
        
        //
        // newt
        //
        if MainLoopPos=newt then
        begin
            if nflg<>0 then
            begin
                Break;
            end;
            nflg := 1;
            lgm := lngamma(aaa+bbb, s)-lngamma(aaa, s)-lngamma(bbb, s);
            i := 0;
            MainLoopPos := newtcycle;
            Continue;
        end;
        
        //
        // newtcycle
        //
        if MainLoopPos=newtcycle then
        begin
            if i<=7 then
            begin
                if i<>0 then
                begin
                    yyy := IncompleteBeta(aaa, bbb, x);
                end;
                if AP_FP_Less(yyy,yl) then
                begin
                    x := x0;
                    yyy := yl;
                end
                else
                begin
                    if AP_FP_Greater(yyy,yh) then
                    begin
                        x := x1;
                        yyy := yh;
                    end
                    else
                    begin
                        if AP_FP_Less(yyy,y0) then
                        begin
                            x0 := x;
                            yl := yyy;
                        end
                        else
                        begin
                            x1 := x;
                            yh := yyy;
                        end;
                    end;
                end;
                if AP_FP_Eq(x,Double(1.0)) or AP_FP_Eq(x,Double(0.0)) then
                begin
                    MainLoopPos := breaknewtcycle;
                    Continue;
                end;
                d := (aaa-Double(1.0))*ln(x)+(bbb-Double(1.0))*ln(Double(1.0)-x)+lgm;
                if AP_FP_Less(d,Ln(MinRealNumber)) then
                begin
                    Break;
                end;
                if AP_FP_Greater(d,Ln(MaxRealNumber)) then
                begin
                    MainLoopPos := breaknewtcycle;
                    Continue;
                end;
                d := exp(d);
                d := (yyy-y0)/d;
                xt := x-d;
                if AP_FP_Less_Eq(xt,x0) then
                begin
                    yyy := (x-x0)/(x1-x0);
                    xt := x0+Double(0.5)*yyy*(x-x0);
                    if AP_FP_Less_Eq(xt,Double(0.0)) then
                    begin
                        MainLoopPos := breaknewtcycle;
                        Continue;
                    end;
                end;
                if AP_FP_Greater_Eq(xt,x1) then
                begin
                    yyy := (x1-x)/(x1-x0);
                    xt := x1-Double(0.5)*yyy*(x1-x);
                    if AP_FP_Greater_Eq(xt,Double(1.0)) then
                    begin
                        MainLoopPos := breaknewtcycle;
                        Continue;
                    end;
                end;
                x := xt;
                if AP_FP_Less(absReal(d/x),Double(128.0)*MachineEpsilon) then
                begin
                    Break;
                end;
                i := i+1;
                MainLoopPos := newtcycle;
                Continue;
            end
            else
            begin
                MainLoopPos := breaknewtcycle;
                Continue;
            end;
        end;
        
        //
        // breaknewtcycle
        //
        if MainLoopPos=breaknewtcycle then
        begin
            dithresh := Double(256.0)*MachineEpsilon;
            MainLoopPos := ihalve;
            Continue;
        end;
    end;
    
    //
    // done
    //
    if rflg<>0 then
    begin
        if AP_FP_Less_Eq(x,MachineEpsilon) then
        begin
            x := Double(1.0)-MachineEpsilon;
        end
        else
        begin
            x := Double(1.0)-x;
        end;
    end;
    Result := x;
end;


(*************************************************************************
Continued fraction expansion #1 for incomplete beta integral

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteBetaFE(a : Double;
     b : Double;
     x : Double;
     big : Double;
     biginv : Double):Double;
var
    xk : Double;
    pk : Double;
    pkm1 : Double;
    pkm2 : Double;
    qk : Double;
    qkm1 : Double;
    qkm2 : Double;
    k1 : Double;
    k2 : Double;
    k3 : Double;
    k4 : Double;
    k5 : Double;
    k6 : Double;
    k7 : Double;
    k8 : Double;
    r : Double;
    t : Double;
    ans : Double;
    thresh : Double;
    n : AlglibInteger;
begin
    k1 := a;
    k2 := a+b;
    k3 := a;
    k4 := a+Double(1.0);
    k5 := Double(1.0);
    k6 := b-Double(1.0);
    k7 := k4;
    k8 := a+Double(2.0);
    pkm2 := Double(0.0);
    qkm2 := Double(1.0);
    pkm1 := Double(1.0);
    qkm1 := Double(1.0);
    ans := Double(1.0);
    r := Double(1.0);
    n := 0;
    thresh := Double(3.0)*MachineEpsilon;
    repeat
        xk := -x*k1*k2/(k3*k4);
        pk := pkm1+pkm2*xk;
        qk := qkm1+qkm2*xk;
        pkm2 := pkm1;
        pkm1 := pk;
        qkm2 := qkm1;
        qkm1 := qk;
        xk := x*k5*k6/(k7*k8);
        pk := pkm1+pkm2*xk;
        qk := qkm1+qkm2*xk;
        pkm2 := pkm1;
        pkm1 := pk;
        qkm2 := qkm1;
        qkm1 := qk;
        if AP_FP_Neq(qk,0) then
        begin
            r := pk/qk;
        end;
        if AP_FP_Neq(r,0) then
        begin
            t := absReal((ans-r)/r);
            ans := r;
        end
        else
        begin
            t := Double(1.0);
        end;
        if AP_FP_Less(t,thresh) then
        begin
            Break;
        end;
        k1 := k1+Double(1.0);
        k2 := k2+Double(1.0);
        k3 := k3+Double(2.0);
        k4 := k4+Double(2.0);
        k5 := k5+Double(1.0);
        k6 := k6-Double(1.0);
        k7 := k7+Double(2.0);
        k8 := k8+Double(2.0);
        if AP_FP_Greater(absReal(qk)+absReal(pk),big) then
        begin
            pkm2 := pkm2*biginv;
            pkm1 := pkm1*biginv;
            qkm2 := qkm2*biginv;
            qkm1 := qkm1*biginv;
        end;
        if AP_FP_Less(absReal(qk),biginv) or AP_FP_Less(absReal(pk),biginv) then
        begin
            pkm2 := pkm2*big;
            pkm1 := pkm1*big;
            qkm2 := qkm2*big;
            qkm1 := qkm1*big;
        end;
        n := n+1;
    until n=300;
    Result := ans;
end;


(*************************************************************************
Continued fraction expansion #2
for incomplete beta integral

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteBetaFE2(a : Double;
     b : Double;
     x : Double;
     big : Double;
     biginv : Double):Double;
var
    xk : Double;
    pk : Double;
    pkm1 : Double;
    pkm2 : Double;
    qk : Double;
    qkm1 : Double;
    qkm2 : Double;
    k1 : Double;
    k2 : Double;
    k3 : Double;
    k4 : Double;
    k5 : Double;
    k6 : Double;
    k7 : Double;
    k8 : Double;
    r : Double;
    t : Double;
    ans : Double;
    z : Double;
    thresh : Double;
    n : AlglibInteger;
begin
    k1 := a;
    k2 := b-Double(1.0);
    k3 := a;
    k4 := a+Double(1.0);
    k5 := Double(1.0);
    k6 := a+b;
    k7 := a+Double(1.0);
    k8 := a+Double(2.0);
    pkm2 := Double(0.0);
    qkm2 := Double(1.0);
    pkm1 := Double(1.0);
    qkm1 := Double(1.0);
    z := x/(Double(1.0)-x);
    ans := Double(1.0);
    r := Double(1.0);
    n := 0;
    thresh := Double(3.0)*MachineEpsilon;
    repeat
        xk := -z*k1*k2/(k3*k4);
        pk := pkm1+pkm2*xk;
        qk := qkm1+qkm2*xk;
        pkm2 := pkm1;
        pkm1 := pk;
        qkm2 := qkm1;
        qkm1 := qk;
        xk := z*k5*k6/(k7*k8);
        pk := pkm1+pkm2*xk;
        qk := qkm1+qkm2*xk;
        pkm2 := pkm1;
        pkm1 := pk;
        qkm2 := qkm1;
        qkm1 := qk;
        if AP_FP_Neq(qk,0) then
        begin
            r := pk/qk;
        end;
        if AP_FP_Neq(r,0) then
        begin
            t := absReal((ans-r)/r);
            ans := r;
        end
        else
        begin
            t := Double(1.0);
        end;
        if AP_FP_Less(t,thresh) then
        begin
            Break;
        end;
        k1 := k1+Double(1.0);
        k2 := k2-Double(1.0);
        k3 := k3+Double(2.0);
        k4 := k4+Double(2.0);
        k5 := k5+Double(1.0);
        k6 := k6+Double(1.0);
        k7 := k7+Double(2.0);
        k8 := k8+Double(2.0);
        if AP_FP_Greater(absReal(qk)+absReal(pk),big) then
        begin
            pkm2 := pkm2*biginv;
            pkm1 := pkm1*biginv;
            qkm2 := qkm2*biginv;
            qkm1 := qkm1*biginv;
        end;
        if AP_FP_Less(absReal(qk),biginv) or AP_FP_Less(absReal(pk),biginv) then
        begin
            pkm2 := pkm2*big;
            pkm1 := pkm1*big;
            qkm2 := qkm2*big;
            qkm1 := qkm1*big;
        end;
        n := n+1;
    until n=300;
    Result := ans;
end;


(*************************************************************************
Power series for incomplete beta integral.
Use when b*x is small and x not too close to 1.

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteBetaPS(a : Double;
     b : Double;
     x : Double;
     MAXGAM : Double):Double;
var
    s : Double;
    t : Double;
    u : Double;
    v : Double;
    n : Double;
    t1 : Double;
    z : Double;
    ai : Double;
    sg : Double;
begin
    ai := Double(1.0)/a;
    u := (Double(1.0)-b)*x;
    v := u/(a+Double(1.0));
    t1 := v;
    t := u;
    n := Double(2.0);
    s := Double(0.0);
    z := MACHinEePsilon*ai;
    while AP_FP_Greater(absReal(v),z) do
    begin
        u := (n-b)*x/n;
        t := t*u;
        v := t/(a+n);
        s := s+v;
        n := n+Double(1.0);
    end;
    s := s+t1;
    s := s+ai;
    u := a*ln(x);
    if AP_FP_Less(a+b,MAXGAM) and AP_FP_Less(absReal(u),Ln(MaxRealNumber)) then
    begin
        t := gamma(a+b)/(gamma(a)*gamma(b));
        s := s*t*power(x, a);
    end
    else
    begin
        t := lngamma(a+b, sg)-lngamma(a, sg)-lngamma(b, sg)+u+ln(s);
        if AP_FP_Less(t,Ln(MinRealNumber)) then
        begin
            s := Double(0.0);
        end
        else
        begin
            s := exp(t);
        end;
    end;
    Result := s;
end;


end.