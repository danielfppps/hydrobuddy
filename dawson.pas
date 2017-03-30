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
unit dawson;
interface
uses Math, Sysutils, Ap;

function DawsonIntegral(X : Double):Double;

implementation

(*************************************************************************
Dawson's Integral

Approximates the integral

                            x
                            -
                     2     | |        2
 dawsn(x)  =  exp( -x  )   |    exp( t  ) dt
                         | |
                          -
                          0

Three different rational approximations are employed, for
the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,10        10000       6.9e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************)
function DawsonIntegral(X : Double):Double;
var
    X2 : Double;
    Y : Double;
    Sg : AlglibInteger;
    AN : Double;
    AD : Double;
    BN : Double;
    BD : Double;
    CN : Double;
    CD : Double;
begin
    Sg := 1;
    if AP_FP_Less(X,0) then
    begin
        Sg := -1;
        X := -X;
    end;
    if AP_FP_Less(X,Double(3.25)) then
    begin
        x2 := x*x;
        AN := Double(1.13681498971755972054E-11);
        AN := AN*x2+Double(8.49262267667473811108E-10);
        AN := AN*x2+Double(1.94434204175553054283E-8);
        AN := AN*x2+Double(9.53151741254484363489E-7);
        AN := AN*x2+Double(3.07828309874913200438E-6);
        AN := AN*x2+Double(3.52513368520288738649E-4);
        AN := AN*x2+-Double(8.50149846724410912031E-4);
        AN := AN*x2+Double(4.22618223005546594270E-2);
        AN := AN*x2+-Double(9.17480371773452345351E-2);
        AN := AN*x2+Double(9.99999999999999994612E-1);
        AD := Double(2.40372073066762605484E-11);
        AD := AD*x2+Double(1.48864681368493396752E-9);
        AD := AD*x2+Double(5.21265281010541664570E-8);
        AD := AD*x2+Double(1.27258478273186970203E-6);
        AD := AD*x2+Double(2.32490249820789513991E-5);
        AD := AD*x2+Double(3.25524741826057911661E-4);
        AD := AD*x2+Double(3.48805814657162590916E-3);
        AD := AD*x2+Double(2.79448531198828973716E-2);
        AD := AD*x2+Double(1.58874241960120565368E-1);
        AD := AD*x2+Double(5.74918629489320327824E-1);
        AD := AD*x2+Double(1.00000000000000000539E0);
        y := x*AN/AD;
        Result := Sg*Y;
        Exit;
    end;
    x2 := Double(1.0)/(x*x);
    if AP_FP_Less(x,Double(6.25)) then
    begin
        BN := Double(5.08955156417900903354E-1);
        BN := BN*x2-Double(2.44754418142697847934E-1);
        BN := BN*x2+Double(9.41512335303534411857E-2);
        BN := BN*x2-Double(2.18711255142039025206E-2);
        BN := BN*x2+Double(3.66207612329569181322E-3);
        BN := BN*x2-Double(4.23209114460388756528E-4);
        BN := BN*x2+Double(3.59641304793896631888E-5);
        BN := BN*x2-Double(2.14640351719968974225E-6);
        BN := BN*x2+Double(9.10010780076391431042E-8);
        BN := BN*x2-Double(2.40274520828250956942E-9);
        BN := BN*x2+Double(3.59233385440928410398E-11);
        BD := Double(1.00000000000000000000E0);
        BD := BD*x2-Double(6.31839869873368190192E-1);
        BD := BD*x2+Double(2.36706788228248691528E-1);
        BD := BD*x2-Double(5.31806367003223277662E-2);
        BD := BD*x2+Double(8.48041718586295374409E-3);
        BD := BD*x2-Double(9.47996768486665330168E-4);
        BD := BD*x2+Double(7.81025592944552338085E-5);
        BD := BD*x2-Double(4.55875153252442634831E-6);
        BD := BD*x2+Double(1.89100358111421846170E-7);
        BD := BD*x2-Double(4.91324691331920606875E-9);
        BD := BD*x2+Double(7.18466403235734541950E-11);
        y := Double(1.0)/x+x2*BN/(BD*x);
        Result := sg*Double(0.5)*y;
        Exit;
    end;
    if AP_FP_Greater(x,Double(1.0E9)) then
    begin
        Result := sg*Double(0.5)/x;
        Exit;
    end;
    CN := -Double(5.90592860534773254987E-1);
    CN := CN*x2+Double(6.29235242724368800674E-1);
    CN := CN*x2-Double(1.72858975380388136411E-1);
    CN := CN*x2+Double(1.64837047825189632310E-2);
    CN := CN*x2-Double(4.86827613020462700845E-4);
    CD := Double(1.00000000000000000000E0);
    CD := CD*x2-Double(2.69820057197544900361E0);
    CD := CD*x2+Double(1.73270799045947845857E0);
    CD := CD*x2-Double(3.93708582281939493482E-1);
    CD := CD*x2+Double(3.44278924041233391079E-2);
    CD := CD*x2-Double(9.73655226040941223894E-4);
    y := Double(1.0)/x+x2*CN/(CD*x);
    Result := sg*Double(0.5)*y;
end;


end.