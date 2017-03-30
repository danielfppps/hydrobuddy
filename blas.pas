{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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
unit blas;
interface
uses Math, Sysutils, Ap;

function VectorNorm2(const X : TReal1DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger):Double;
function VectorIdxAbsMax(const X : TReal1DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger):AlglibInteger;
function ColumnIdxAbsMax(const X : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J : AlglibInteger):AlglibInteger;
function RowIdxAbsMax(const X : TReal2DArray;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     I : AlglibInteger):AlglibInteger;
function UpperHessenberg1Norm(const A : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     var WORK : TReal1DArray):Double;
procedure CopyMatrix(const A : TReal2DArray;
     IS1 : AlglibInteger;
     IS2 : AlglibInteger;
     JS1 : AlglibInteger;
     JS2 : AlglibInteger;
     var B : TReal2DArray;
     ID1 : AlglibInteger;
     ID2 : AlglibInteger;
     JD1 : AlglibInteger;
     JD2 : AlglibInteger);
procedure InplaceTranspose(var A : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     var WORK : TReal1DArray);
procedure CopyAndTranspose(const A : TReal2DArray;
     IS1 : AlglibInteger;
     IS2 : AlglibInteger;
     JS1 : AlglibInteger;
     JS2 : AlglibInteger;
     var B : TReal2DArray;
     ID1 : AlglibInteger;
     ID2 : AlglibInteger;
     JD1 : AlglibInteger;
     JD2 : AlglibInteger);
procedure MatrixVectorMultiply(const A : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     Trans : Boolean;
     const X : TReal1DArray;
     IX1 : AlglibInteger;
     IX2 : AlglibInteger;
     Alpha : Double;
     var Y : TReal1DArray;
     IY1 : AlglibInteger;
     IY2 : AlglibInteger;
     Beta : Double);
function Pythag2(X : Double; Y : Double):Double;
procedure MatrixMatrixMultiply(const A : TReal2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TReal2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     Alpha : Double;
     var C : TReal2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger;
     Beta : Double;
     var WORK : TReal1DArray);

implementation

function VectorNorm2(const X : TReal1DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger):Double;
var
    N : AlglibInteger;
    IX : AlglibInteger;
    ABSXI : Double;
    SCL : Double;
    SSQ : Double;
begin
    N := I2-I1+1;
    if N<1 then
    begin
        Result := 0;
        Exit;
    end;
    if N=1 then
    begin
        Result := AbsReal(X[I1]);
        Exit;
    end;
    SCL := 0;
    SSQ := 1;
    IX:=I1;
    while IX<=I2 do
    begin
        if AP_FP_Neq(X[IX],0) then
        begin
            ABSXI := ABSReal(X[IX]);
            if AP_FP_Less(SCL,ABSXI) then
            begin
                SSQ := 1+SSQ*AP_Sqr(SCL/ABSXI);
                SCL := ABSXI;
            end
            else
            begin
                SSQ := SSQ+AP_Sqr(ABSXI/SCL);
            end;
        end;
        Inc(IX);
    end;
    Result := SCL*SQRT(SSQ);
end;


function VectorIdxAbsMax(const X : TReal1DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger):AlglibInteger;
var
    I : AlglibInteger;
    A : Double;
begin
    Result := I1;
    A := AbsReal(X[Result]);
    I:=I1+1;
    while I<=I2 do
    begin
        if AP_FP_Greater(AbsReal(X[I]),AbsReal(X[Result])) then
        begin
            Result := I;
        end;
        Inc(I);
    end;
end;


function ColumnIdxAbsMax(const X : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J : AlglibInteger):AlglibInteger;
var
    I : AlglibInteger;
    A : Double;
begin
    Result := I1;
    A := AbsReal(X[Result,J]);
    I:=I1+1;
    while I<=I2 do
    begin
        if AP_FP_Greater(AbsReal(X[I,J]),AbsReal(X[Result,J])) then
        begin
            Result := I;
        end;
        Inc(I);
    end;
end;


function RowIdxAbsMax(const X : TReal2DArray;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     I : AlglibInteger):AlglibInteger;
var
    J : AlglibInteger;
    A : Double;
begin
    Result := J1;
    A := AbsReal(X[I,Result]);
    J:=J1+1;
    while J<=J2 do
    begin
        if AP_FP_Greater(AbsReal(X[I,J]),AbsReal(X[I,Result])) then
        begin
            Result := J;
        end;
        Inc(J);
    end;
end;


function UpperHessenberg1Norm(const A : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     var WORK : TReal1DArray):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(I2-I1=J2-J1, 'UpperHessenberg1Norm: I2-I1<>J2-J1!');
    J:=J1;
    while J<=J2 do
    begin
        WORK[J] := 0;
        Inc(J);
    end;
    I:=I1;
    while I<=I2 do
    begin
        J:=Max(J1, J1+I-I1-1);
        while J<=J2 do
        begin
            WORK[J] := WORK[J]+AbsReal(A[I,J]);
            Inc(J);
        end;
        Inc(I);
    end;
    Result := 0;
    J:=J1;
    while J<=J2 do
    begin
        Result := Max(Result, WORK[J]);
        Inc(J);
    end;
end;


procedure CopyMatrix(const A : TReal2DArray;
     IS1 : AlglibInteger;
     IS2 : AlglibInteger;
     JS1 : AlglibInteger;
     JS2 : AlglibInteger;
     var B : TReal2DArray;
     ID1 : AlglibInteger;
     ID2 : AlglibInteger;
     JD1 : AlglibInteger;
     JD2 : AlglibInteger);
var
    ISRC : AlglibInteger;
    IDST : AlglibInteger;
begin
    if (IS1>IS2) or (JS1>JS2) then
    begin
        Exit;
    end;
    Assert(IS2-IS1=ID2-ID1, 'CopyMatrix: different sizes!');
    Assert(JS2-JS1=JD2-JD1, 'CopyMatrix: different sizes!');
    ISRC:=IS1;
    while ISRC<=IS2 do
    begin
        IDST := ISRC-IS1+ID1;
        APVMove(@B[IDST][0], JD1, JD2, @A[ISRC][0], JS1, JS2);
        Inc(ISRC);
    end;
end;


procedure InplaceTranspose(var A : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     var WORK : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    IPS : AlglibInteger;
    JPS : AlglibInteger;
    L : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (I1>I2) or (J1>J2) then
    begin
        Exit;
    end;
    Assert(I1-I2=J1-J2, 'InplaceTranspose error: incorrect array size!');
    I:=I1;
    while I<=I2-1 do
    begin
        J := J1+I-I1;
        IPS := I+1;
        JPS := J1+IPS-I1;
        L := I2-I;
        i1_ := (IPS) - (1);
        for i_ := 1 to L do
        begin
            WORK[i_] := A[i_+i1_,J];
        end;
        i1_ := (JPS) - (IPS);
        for i_ := IPS to I2 do
        begin
            A[i_,J] := A[I,i_+i1_];
        end;
        APVMove(@A[I][0], JPS, J2, @WORK[0], 1, L);
        Inc(I);
    end;
end;


procedure CopyAndTranspose(const A : TReal2DArray;
     IS1 : AlglibInteger;
     IS2 : AlglibInteger;
     JS1 : AlglibInteger;
     JS2 : AlglibInteger;
     var B : TReal2DArray;
     ID1 : AlglibInteger;
     ID2 : AlglibInteger;
     JD1 : AlglibInteger;
     JD2 : AlglibInteger);
var
    ISRC : AlglibInteger;
    JDST : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (IS1>IS2) or (JS1>JS2) then
    begin
        Exit;
    end;
    Assert(IS2-IS1=JD2-JD1, 'CopyAndTranspose: different sizes!');
    Assert(JS2-JS1=ID2-ID1, 'CopyAndTranspose: different sizes!');
    ISRC:=IS1;
    while ISRC<=IS2 do
    begin
        JDST := ISRC-IS1+JD1;
        i1_ := (JS1) - (ID1);
        for i_ := ID1 to ID2 do
        begin
            B[i_,JDST] := A[ISRC,i_+i1_];
        end;
        Inc(ISRC);
    end;
end;


procedure MatrixVectorMultiply(const A : TReal2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     Trans : Boolean;
     const X : TReal1DArray;
     IX1 : AlglibInteger;
     IX2 : AlglibInteger;
     Alpha : Double;
     var Y : TReal1DArray;
     IY1 : AlglibInteger;
     IY2 : AlglibInteger;
     Beta : Double);
var
    I : AlglibInteger;
    V : Double;
begin
    if  not Trans then
    begin
        
        //
        // y := alpha*A*x + beta*y;
        //
        if (I1>I2) or (J1>J2) then
        begin
            Exit;
        end;
        Assert(J2-J1=IX2-IX1, 'MatrixVectorMultiply: A and X dont match!');
        Assert(I2-I1=IY2-IY1, 'MatrixVectorMultiply: A and Y dont match!');
        
        //
        // beta*y
        //
        if AP_FP_Eq(Beta,0) then
        begin
            I:=IY1;
            while I<=IY2 do
            begin
                Y[I] := 0;
                Inc(I);
            end;
        end
        else
        begin
            APVMul(@Y[0], IY1, IY2, Beta);
        end;
        
        //
        // alpha*A*x
        //
        I:=I1;
        while I<=I2 do
        begin
            V := APVDotProduct(@A[I][0], J1, J2, @X[0], IX1, IX2);
            Y[IY1+I-I1] := Y[IY1+I-I1]+Alpha*V;
            Inc(I);
        end;
    end
    else
    begin
        
        //
        // y := alpha*A'*x + beta*y;
        //
        if (I1>I2) or (J1>J2) then
        begin
            Exit;
        end;
        Assert(I2-I1=IX2-IX1, 'MatrixVectorMultiply: A and X dont match!');
        Assert(J2-J1=IY2-IY1, 'MatrixVectorMultiply: A and Y dont match!');
        
        //
        // beta*y
        //
        if AP_FP_Eq(Beta,0) then
        begin
            I:=IY1;
            while I<=IY2 do
            begin
                Y[I] := 0;
                Inc(I);
            end;
        end
        else
        begin
            APVMul(@Y[0], IY1, IY2, Beta);
        end;
        
        //
        // alpha*A'*x
        //
        I:=I1;
        while I<=I2 do
        begin
            V := Alpha*X[IX1+I-I1];
            APVAdd(@Y[0], IY1, IY2, @A[I][0], J1, J2, V);
            Inc(I);
        end;
    end;
end;


function Pythag2(X : Double; Y : Double):Double;
var
    W : Double;
    XABS : Double;
    YABS : Double;
    Z : Double;
begin
    XABS := AbsReal(X);
    YABS := AbsReal(Y);
    W := Max(XABS, YABS);
    Z := Min(XABS, YABS);
    if AP_FP_Eq(Z,0) then
    begin
        Result := W;
    end
    else
    begin
        Result := W*SQRT(1+AP_Sqr(Z/W));
    end;
end;


procedure MatrixMatrixMultiply(const A : TReal2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TReal2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     Alpha : Double;
     var C : TReal2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger;
     Beta : Double;
     var WORK : TReal1DArray);
var
    ARows : AlglibInteger;
    ACols : AlglibInteger;
    BRows : AlglibInteger;
    BCols : AlglibInteger;
    CRows : AlglibInteger;
    CCols : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Setup
    //
    if  not TransA then
    begin
        ARows := AI2-AI1+1;
        ACols := AJ2-AJ1+1;
    end
    else
    begin
        ARows := AJ2-AJ1+1;
        ACols := AI2-AI1+1;
    end;
    if  not TransB then
    begin
        BRows := BI2-BI1+1;
        BCols := BJ2-BJ1+1;
    end
    else
    begin
        BRows := BJ2-BJ1+1;
        BCols := BI2-BI1+1;
    end;
    Assert(ACols=BRows, 'MatrixMatrixMultiply: incorrect matrix sizes!');
    if (ARows<=0) or (ACols<=0) or (BRows<=0) or (BCols<=0) then
    begin
        Exit;
    end;
    CRows := ARows;
    CCols := BCols;
    
    //
    // Test WORK
    //
    I := Max(ARows, ACols);
    I := Max(BRows, I);
    I := Max(I, BCols);
    Work[1] := 0;
    Work[I] := 0;
    
    //
    // Prepare C
    //
    if AP_FP_Eq(Beta,0) then
    begin
        I:=CI1;
        while I<=CI2 do
        begin
            J:=CJ1;
            while J<=CJ2 do
            begin
                C[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=CI1;
        while I<=CI2 do
        begin
            APVMul(@C[I][0], CJ1, CJ2, Beta);
            Inc(I);
        end;
    end;
    
    //
    // A*B
    //
    if  not TransA and  not TransB then
    begin
        L:=AI1;
        while L<=AI2 do
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                V := Alpha*A[L,AJ1+R-BI1];
                K := CI1+L-AI1;
                APVAdd(@C[K][0], CJ1, CJ2, @B[R][0], BJ1, BJ2, V);
                Inc(R);
            end;
            Inc(L);
        end;
        Exit;
    end;
    
    //
    // A*B'
    //
    if  not TransA and TransB then
    begin
        if ARows*ACols<BRows*BCols then
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                L:=AI1;
                while L<=AI2 do
                begin
                    V := APVDotProduct(@A[L][0], AJ1, AJ2, @B[R][0], BJ1, BJ2);
                    C[CI1+L-AI1,CJ1+R-BI1] := C[CI1+L-AI1,CJ1+R-BI1]+Alpha*V;
                    Inc(L);
                end;
                Inc(R);
            end;
            Exit;
        end
        else
        begin
            L:=AI1;
            while L<=AI2 do
            begin
                R:=BI1;
                while R<=BI2 do
                begin
                    V := APVDotProduct(@A[L][0], AJ1, AJ2, @B[R][0], BJ1, BJ2);
                    C[CI1+L-AI1,CJ1+R-BI1] := C[CI1+L-AI1,CJ1+R-BI1]+Alpha*V;
                    Inc(R);
                end;
                Inc(L);
            end;
            Exit;
        end;
    end;
    
    //
    // A'*B
    //
    if TransA and  not TransB then
    begin
        L:=AJ1;
        while L<=AJ2 do
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                V := Alpha*A[AI1+R-BI1,L];
                K := CI1+L-AJ1;
                APVAdd(@C[K][0], CJ1, CJ2, @B[R][0], BJ1, BJ2, V);
                Inc(R);
            end;
            Inc(L);
        end;
        Exit;
    end;
    
    //
    // A'*B'
    //
    if TransA and TransB then
    begin
        if ARows*ACols<BRows*BCols then
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                I:=1;
                while I<=CRows do
                begin
                    WORK[I] := Double(0.0);
                    Inc(I);
                end;
                L:=AI1;
                while L<=AI2 do
                begin
                    V := Alpha*B[R,BJ1+L-AI1];
                    K := CJ1+R-BI1;
                    APVAdd(@WORK[0], 1, CRows, @A[L][0], AJ1, AJ2, V);
                    Inc(L);
                end;
                i1_ := (1) - (CI1);
                for i_ := CI1 to CI2 do
                begin
                    C[i_,K] := C[i_,K] + WORK[i_+i1_];
                end;
                Inc(R);
            end;
            Exit;
        end
        else
        begin
            L:=AJ1;
            while L<=AJ2 do
            begin
                K := AI2-AI1+1;
                i1_ := (AI1) - (1);
                for i_ := 1 to K do
                begin
                    WORK[i_] := A[i_+i1_,L];
                end;
                R:=BI1;
                while R<=BI2 do
                begin
                    V := APVDotProduct(@WORK[0], 1, K, @B[R][0], BJ1, BJ2);
                    C[CI1+L-AJ1,CJ1+R-BI1] := C[CI1+L-AJ1,CJ1+R-BI1]+Alpha*V;
                    Inc(R);
                end;
                Inc(L);
            end;
            Exit;
        end;
    end;
end;


end.