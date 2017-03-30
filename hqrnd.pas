{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c)
    2007, Sergey Bochkanov (ALGLIB project).
    1988, Pierre L'Ecuyer

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
unit hqrnd;
interface
uses Math, Sysutils, Ap;

type
(*************************************************************************
Portable high quality random number generator state.
Initialized with HQRNDRandomize() or HQRNDSeed().

Fields:
    S1, S2      -   seed values
    V           -   precomputed value
    MagicV      -   'magic' value used to determine whether State structure
                    was correctly initialized.
*************************************************************************)
HQRNDState = record
    S1 : AlglibInteger;
    S2 : AlglibInteger;
    V : Double;
    MagicV : AlglibInteger;
end;



procedure HQRNDRandomize(var State : HQRNDState);
procedure HQRNDSeed(S1 : AlglibInteger;
     S2 : AlglibInteger;
     var State : HQRNDState);
function HQRNDUniformR(var State : HQRNDState):Double;
function HQRNDUniformI(N : AlglibInteger;
     var State : HQRNDState):AlglibInteger;
function HQRNDNormal(var State : HQRNDState):Double;
procedure HQRNDUnit2(var State : HQRNDState; var X : Double; var Y : Double);
procedure HQRNDNormal2(var State : HQRNDState;
     var X1 : Double;
     var X2 : Double);
function HQRNDExponential(Lambda : Double; var State : HQRNDState):Double;

implementation

const
    HQRNDMax = 2147483563;
    HQRNDM1 = 2147483563;
    HQRNDM2 = 2147483399;
    HQRNDMagic = 1634357784;

function HQRNDIntegerBase(var State : HQRNDState):AlglibInteger;forward;


(*************************************************************************
HQRNDState  initialization  with  random  values  which come from standard
RNG.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure HQRNDRandomize(var State : HQRNDState);
begin
    HQRNDSeed(RandomInteger(HQRNDM1), RandomInteger(HQRNDM2), State);
end;


(*************************************************************************
HQRNDState initialization with seed values

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure HQRNDSeed(S1 : AlglibInteger;
     S2 : AlglibInteger;
     var State : HQRNDState);
begin
    State.S1 := S1 mod (HQRNDM1-1)+1;
    State.S2 := S2 mod (HQRNDM2-1)+1;
    State.V := AP_Double(1)/HQRNDMax;
    State.MagicV := HQRNDMagic;
end;


(*************************************************************************
This function generates random real number in (0,1),
not including interval boundaries

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function HQRNDUniformR(var State : HQRNDState):Double;
begin
    Result := State.V*HQRNDIntegerBase(State);
end;


(*************************************************************************
This function generates random integer number in [0, N)

1. N must be less than HQRNDMax-1.
2. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function HQRNDUniformI(N : AlglibInteger;
     var State : HQRNDState):AlglibInteger;
var
    MX : AlglibInteger;
begin
    
    //
    // Correct handling of N's close to RNDBaseMax
    // (avoiding skewed distributions for RNDBaseMax<>K*N)
    //
    Assert(N>0, 'HQRNDUniformI: N<=0!');
    Assert(N<HQRNDMax-1, 'HQRNDUniformI: N>=RNDBaseMax-1!');
    MX := HQRNDMax-1-(HQRNDMax-1) mod N;
    repeat
        Result := HQRNDIntegerBase(State)-1;
    until Result<MX;
    Result := Result mod N;
end;


(*************************************************************************
Random number generator: normal numbers

This function generates one random number from normal distribution.
Its performance is equal to that of HQRNDNormal2()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
function HQRNDNormal(var State : HQRNDState):Double;
var
    V1 : Double;
    V2 : Double;
begin
    HQRNDNormal2(State, V1, V2);
    Result := V1;
end;


(*************************************************************************
Random number generator: random X and Y such that X^2+Y^2=1

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure HQRNDUnit2(var State : HQRNDState; var X : Double; var Y : Double);
var
    V : Double;
    MX : Double;
    MN : Double;
begin
    repeat
        HQRNDNormal2(State, X, Y);
    until AP_FP_Neq(X,0) or AP_FP_Neq(Y,0);
    MX := Max(AbsReal(X), AbsReal(Y));
    MN := Min(AbsReal(X), AbsReal(Y));
    V := MX*Sqrt(1+AP_Sqr(MN/MX));
    X := X/V;
    Y := Y/V;
end;


(*************************************************************************
Random number generator: normal numbers

This function generates two independent random numbers from normal
distribution. Its performance is equal to that of HQRNDNormal()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure HQRNDNormal2(var State : HQRNDState;
     var X1 : Double;
     var X2 : Double);
var
    U : Double;
    V : Double;
    S : Double;
begin
    while True do
    begin
        U := 2*HQRNDUniformR(State)-1;
        V := 2*HQRNDUniformR(State)-1;
        S := AP_Sqr(u)+AP_Sqr(v);
        if AP_FP_Greater(S,0) and AP_FP_Less(S,1) then
        begin
            
            //
            // two Sqrt's instead of one to
            // avoid overflow when S is too small
            //
            S := Sqrt(-2*Ln(S))/Sqrt(S);
            X1 := U*S;
            X2 := V*S;
            Exit;
        end;
    end;
end;


(*************************************************************************
Random number generator: exponential distribution

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 11.08.2007 by Bochkanov Sergey
*************************************************************************)
function HQRNDExponential(Lambda : Double; var State : HQRNDState):Double;
begin
    Assert(AP_FP_Greater(Lambda,0), 'HQRNDExponential: Lambda<=0!');
    Result := -Ln(HQRNDUniformR(State))/Lambda;
end;


(*************************************************************************

L'Ecuyer, Efficient and portable combined random number generators
*************************************************************************)
function HQRNDIntegerBase(var State : HQRNDState):AlglibInteger;
var
    K : AlglibInteger;
begin
    Assert(State.MagicV=HQRNDMagic, 'HQRNDIntegerBase: State is not correctly initialized!');
    K := State.S1 div 53668;
    State.S1 := 40014*(State.S1-K*53668)-K*12211;
    if State.S1<0 then
    begin
        State.S1 := State.S1+2147483563;
    end;
    K := State.S2 div 52774;
    State.S2 := 40692*(State.S2-K*52774)-K*3791;
    if State.S2<0 then
    begin
        State.S2 := State.S2+2147483399;
    end;
    
    //
    // Result
    //
    Result := State.S1-State.S2;
    if Result<1 then
    begin
        Result := Result+2147483562;
    end;
end;


end.