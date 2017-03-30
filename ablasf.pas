{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
unit ablasf;
interface
uses Math, Sysutils, Ap;

function CMatrixRank1F(M : AlglibInteger;
     N : AlglibInteger;
     var A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     var U : TComplex1DArray;
     IU : AlglibInteger;
     var V : TComplex1DArray;
     IV : AlglibInteger):Boolean;
function RMatrixRank1F(M : AlglibInteger;
     N : AlglibInteger;
     var A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     var U : TReal1DArray;
     IU : AlglibInteger;
     var V : TReal1DArray;
     IV : AlglibInteger):Boolean;
function CMatrixMVF(M : AlglibInteger;
     N : AlglibInteger;
     var A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpA : AlglibInteger;
     var X : TComplex1DArray;
     IX : AlglibInteger;
     var Y : TComplex1DArray;
     IY : AlglibInteger):Boolean;
function RMatrixMVF(M : AlglibInteger;
     N : AlglibInteger;
     var A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpA : AlglibInteger;
     var X : TReal1DArray;
     IX : AlglibInteger;
     var Y : TReal1DArray;
     IY : AlglibInteger):Boolean;
function CMatrixRightTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
function CMatrixLeftTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
function RMatrixRightTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
function RMatrixLeftTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
function CMatrixSYRKF(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean):Boolean;
function RMatrixSYRKF(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean):Boolean;
function RMatrixGEMMF(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TReal2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger):Boolean;
function CMatrixGEMMF(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Complex;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TComplex2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Complex;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger):Boolean;

implementation

(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixRank1F(M : AlglibInteger;
     N : AlglibInteger;
     var A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     var U : TComplex1DArray;
     IU : AlglibInteger;
     var V : TComplex1DArray;
     IV : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixRank1F(M : AlglibInteger;
     N : AlglibInteger;
     var A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     var U : TReal1DArray;
     IU : AlglibInteger;
     var V : TReal1DArray;
     IV : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixMVF(M : AlglibInteger;
     N : AlglibInteger;
     var A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpA : AlglibInteger;
     var X : TComplex1DArray;
     IX : AlglibInteger;
     var Y : TComplex1DArray;
     IY : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixMVF(M : AlglibInteger;
     N : AlglibInteger;
     var A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpA : AlglibInteger;
     var X : TReal1DArray;
     IX : AlglibInteger;
     var Y : TReal1DArray;
     IY : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixRightTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixLeftTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixRightTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixLeftTRSMF(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixSYRKF(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixSYRKF(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixGEMMF(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TReal2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger):Boolean;
begin
    Result := False;
end;


(*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixGEMMF(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Complex;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TComplex2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Complex;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger):Boolean;
begin
    Result := False;
end;


end.