{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier

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
unit jacobianelliptic;
interface
uses Math, Sysutils, Ap;

procedure JacobianEllipticFunctions(u : Double;
     m : Double;
     var sn : Double;
     var cn : Double;
     var dn : Double;
     var ph : Double);

implementation

(*************************************************************************
Jacobian Elliptic Functions

Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
and dn(u|m) of parameter m between 0 and 1, and real
argument u.

These functions are periodic, with quarter-period on the
real axis equal to the complete elliptic integral
ellpk(1.0-m).

Relation to incomplete elliptic integral:
If u = ellik(phi,m), then sn(u|m) = sin(phi),
and cn(u|m) = cos(phi).  Phi is called the amplitude of u.

Computation is by means of the arithmetic-geometric mean
algorithm, except when m is within 1e-9 of 0 or 1.  In the
latter case with m close to 1, the approximation applies
only for phi < pi/2.

ACCURACY:

Tested at random points with u between 0 and 10, m between
0 and 1.

           Absolute error (* = relative error):
arithmetic   function   # trials      peak         rms
   IEEE      phi         10000       9.2e-16*    1.4e-16*
   IEEE      sn          50000       4.1e-15     4.6e-16
   IEEE      cn          40000       3.6e-15     4.4e-16
   IEEE      dn          10000       1.3e-12     1.8e-14

 Peak error observed in consistency check using addition
theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
the above relation to the incomplete elliptic integral.
Accuracy deteriorates when u is large.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
procedure JacobianEllipticFunctions(u : Double;
     m : Double;
     var sn : Double;
     var cn : Double;
     var dn : Double;
     var ph : Double);
var
    ai : Double;
    b : Double;
    phi : Double;
    t : Double;
    twon : Double;
    a : TReal1DArray;
    c : TReal1DArray;
    i : AlglibInteger;
begin
    Assert(AP_FP_Greater_Eq(m,0) and AP_FP_Less_Eq(m,1), 'Domain error in JacobianEllipticFunctions: m<0 or m>1');
    SetLength(a, 8+1);
    SetLength(c, 8+1);
    if AP_FP_Less(m,Double(1.0e-9)) then
    begin
        t := sin(u);
        b := cos(u);
        ai := Double(0.25)*m*(u-t*b);
        sn := t-ai*b;
        cn := b+ai*t;
        ph := u-ai;
        dn := Double(1.0)-Double(0.5)*m*t*t;
        Exit;
    end;
    if AP_FP_Greater_Eq(m,Double(0.9999999999)) then
    begin
        ai := Double(0.25)*(Double(1.0)-m);
        b := cosh(u);
        t := tanh(u);
        phi := Double(1.0)/b;
        twon := b*sinh(u);
        sn := t+ai*(twon-u)/(b*b);
        ph := Double(2.0)*arctan(exp(u))-Double(1.57079632679489661923)+ai*(twon-u)/b;
        ai := ai*t*phi;
        cn := phi-ai*(twon-u);
        dn := phi+ai*(twon+u);
        Exit;
    end;
    a[0] := Double(1.0);
    b := sqrt(Double(1.0)-m);
    c[0] := sqrt(m);
    twon := Double(1.0);
    i := 0;
    while AP_FP_Greater(AbsReal(c[i]/a[i]),MachineEpsilon) do
    begin
        if i>7 then
        begin
            Assert(False, 'Overflow in JacobianEllipticFunctions');
            Break;
        end;
        ai := a[i];
        i := i+1;
        c[i] := Double(0.5)*(ai-b);
        t := sqrt(ai*b);
        a[i] := Double(0.5)*(ai+b);
        b := t;
        twon := twon*Double(2.0);
    end;
    phi := twon*a[i]*u;
    repeat
        t := c[i]*sin(phi)/a[i];
        b := phi;
        phi := (arcsin(t)+phi)/Double(2.0);
        i := i-1;
    until i=0;
    sn := sin(phi);
    t := cos(phi);
    cn := t;
    dn := t/cos(phi-b);
    ph := phi;
end;


end.