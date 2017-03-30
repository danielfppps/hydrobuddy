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
unit igammaf;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr;

function IncompleteGamma(a : Double; x : Double):Double;
function IncompleteGammaC(a : Double; x : Double):Double;
function InvIncompleteGammaC(a : Double; y0 : Double):Double;

implementation

(*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteGamma(a : Double; x : Double):Double;
var
    IGammaEpsilon : Double;
    ans : Double;
    ax : Double;
    c : Double;
    r : Double;
    Tmp : Double;
begin
    IGammaEpsilon := Double(0.000000000000001);
    if AP_FP_Less_Eq(x,0) or AP_FP_Less_Eq(a,0) then
    begin
        Result := 0;
        Exit;
    end;
    if AP_FP_Greater(x,1) and AP_FP_Greater(x,a) then
    begin
        Result := 1-IncompleteGammaC(a, x);
        Exit;
    end;
    ax := a*Ln(x)-x-LnGamma(a, Tmp);
    if AP_FP_Less(ax,-Double(709.78271289338399)) then
    begin
        Result := 0;
        Exit;
    end;
    ax := Exp(ax);
    r := a;
    c := 1;
    ans := 1;
    repeat
        r := r+1;
        c := c*x/r;
        ans := ans+c;
    until AP_FP_Less_Eq(c/ans,IGammaEpsilon);
    Result := ans*ax/a;
end;


(*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function IncompleteGammaC(a : Double; x : Double):Double;
var
    IGammaEpsilon : Double;
    IGammaBigNumber : Double;
    IGammaBigNumberInv : Double;
    ans : Double;
    ax : Double;
    c : Double;
    yc : Double;
    r : Double;
    t : Double;
    y : Double;
    z : Double;
    pk : Double;
    pkm1 : Double;
    pkm2 : Double;
    qk : Double;
    qkm1 : Double;
    qkm2 : Double;
    Tmp : Double;
begin
    IGammaEpsilon := Double(0.000000000000001);
    IGammaBigNumber := Double(4503599627370496.0);
    IGammaBigNumberInv := Double(2.22044604925031308085)*Double(0.0000000000000001);
    if AP_FP_Less_Eq(x,0) or AP_FP_Less_Eq(a,0) then
    begin
        Result := 1;
        Exit;
    end;
    if AP_FP_Less(x,1) or AP_FP_Less(x,a) then
    begin
        Result := 1-IncompleteGamma(a, x);
        Exit;
    end;
    ax := a*Ln(x)-x-LnGamma(a, Tmp);
    if AP_FP_Less(ax,-Double(709.78271289338399)) then
    begin
        Result := 0;
        Exit;
    end;
    ax := Exp(ax);
    y := 1-a;
    z := x+y+1;
    c := 0;
    pkm2 := 1;
    qkm2 := x;
    pkm1 := x+1;
    qkm1 := z*x;
    ans := pkm1/qkm1;
    repeat
        c := c+1;
        y := y+1;
        z := z+2;
        yc := y*c;
        pk := pkm1*z-pkm2*yc;
        qk := qkm1*z-qkm2*yc;
        if AP_FP_Neq(qk,0) then
        begin
            r := pk/qk;
            t := AbsReal((ans-r)/r);
            ans := r;
        end
        else
        begin
            t := 1;
        end;
        pkm2 := pkm1;
        pkm1 := pk;
        qkm2 := qkm1;
        qkm1 := qk;
        if AP_FP_Greater(AbsReal(pk),IGammaBigNumber) then
        begin
            pkm2 := pkm2*IGammaBigNumberInv;
            pkm1 := pkm1*IGammaBigNumberInv;
            qkm2 := qkm2*IGammaBigNumberInv;
            qkm1 := qkm1*IGammaBigNumberInv;
        end;
    until AP_FP_Less_Eq(t,IGammaEpsilon);
    Result := ans*ax;
end;


(*************************************************************************
Inverse of complemented imcomplete gamma integral

Given p, the function finds x such that

 igamc( a, x ) = p.

Starting with the approximate value

        3
 x = a t

 where

 t = 1 - d - ndtri(p) sqrt(d)

and

 d = 1/9a,

the routine performs up to 10 Newton iterations to find the
root of igamc(a,x) - p = 0.

ACCURACY:

Tested at random a, p in the intervals indicated.

               a        p                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
   IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
   IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function InvIncompleteGammaC(a : Double; y0 : Double):Double;
var
    IGammaEpsilon : Double;
    IInvGammaBigNumber : Double;
    x0 : Double;
    x1 : Double;
    x : Double;
    yl : Double;
    yh : Double;
    y : Double;
    d : Double;
    lgm : Double;
    dithresh : Double;
    i : AlglibInteger;
    dir : AlglibInteger;
    Tmp : Double;
begin
    IGammaEpsilon := Double(0.000000000000001);
    IInvGammaBigNumber := Double(4503599627370496.0);
    x0 := IInvGammaBigNumber;
    yl := 0;
    x1 := 0;
    yh := 1;
    dithresh := 5*IGammaEpsilon;
    d := 1/(9*a);
    y := 1-d-InvNormalDistribution(y0)*Sqrt(d);
    x := a*y*y*y;
    lgm := LnGamma(a, Tmp);
    i := 0;
    while i<10 do
    begin
        if AP_FP_Greater(x,x0) or AP_FP_Less(x,x1) then
        begin
            d := Double(0.0625);
            Break;
        end;
        y := IncompleteGammaC(a, x);
        if AP_FP_Less(y,yl) or AP_FP_Greater(y,yh) then
        begin
            d := Double(0.0625);
            Break;
        end;
        if AP_FP_Less(y,y0) then
        begin
            x0 := x;
            yl := y;
        end
        else
        begin
            x1 := x;
            yh := y;
        end;
        d := (a-1)*Ln(x)-x-lgm;
        if AP_FP_Less(d,-Double(709.78271289338399)) then
        begin
            d := Double(0.0625);
            Break;
        end;
        d := -Exp(d);
        d := (y-y0)/d;
        if AP_FP_Less(AbsReal(d/x),IGammaEpsilon) then
        begin
            Result := x;
            Exit;
        end;
        x := x-d;
        i := i+1;
    end;
    if AP_FP_Eq(x0,IInvGammaBigNumber) then
    begin
        if AP_FP_Less_Eq(x,0) then
        begin
            x := 1;
        end;
        while AP_FP_Eq(x0,IInvGammaBigNumber) do
        begin
            x := (1+d)*x;
            y := IncompleteGammaC(a, x);
            if AP_FP_Less(y,y0) then
            begin
                x0 := x;
                yl := y;
                Break;
            end;
            d := d+d;
        end;
    end;
    d := Double(0.5);
    dir := 0;
    i := 0;
    while i<400 do
    begin
        x := x1+d*(x0-x1);
        y := IncompleteGammaC(a, x);
        lgm := (x0-x1)/(x1+x0);
        if AP_FP_Less(AbsReal(lgm),dithresh) then
        begin
            Break;
        end;
        lgm := (y-y0)/y0;
        if AP_FP_Less(AbsReal(lgm),dithresh) then
        begin
            Break;
        end;
        if AP_FP_Less_Eq(x,Double(0.0)) then
        begin
            Break;
        end;
        if AP_FP_Greater_Eq(y,y0) then
        begin
            x1 := x;
            yh := y;
            if dir<0 then
            begin
                dir := 0;
                d := Double(0.5);
            end
            else
            begin
                if dir>1 then
                begin
                    d := Double(0.5)*d+Double(0.5);
                end
                else
                begin
                    d := (y0-yl)/(yh-yl);
                end;
            end;
            dir := dir+1;
        end
        else
        begin
            x0 := x;
            yl := y;
            if dir>0 then
            begin
                dir := 0;
                d := Double(0.5);
            end
            else
            begin
                if dir<-1 then
                begin
                    d := Double(0.5)*d;
                end
                else
                begin
                    d := (y0-yl)/(yh-yl);
                end;
            end;
            dir := dir-1;
        end;
        i := i+1;
    end;
    Result := x;
end;


end.