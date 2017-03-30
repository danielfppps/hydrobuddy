{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
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
unit rotations;
interface
uses Math, Sysutils, Ap;

procedure ApplyRotationsFromTheLeft(IsForward : Boolean;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
procedure ApplyRotationsFromTheRight(IsForward : Boolean;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
procedure GenerateRotation(F : Double;
     G : Double;
     var CS : Double;
     var SN : Double;
     var R : Double);

implementation

(*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [N1..N2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************)
procedure ApplyRotationsFromTheLeft(IsForward : Boolean;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
var
    J : AlglibInteger;
    JP1 : AlglibInteger;
    CTEMP : Double;
    STEMP : Double;
    TEMP : Double;
begin
    if (M1>M2) or (N1>N2) then
    begin
        Exit;
    end;
    
    //
    // Form  P * A
    //
    if IsForward then
    begin
        if N1<>N2 then
        begin
            
            //
            // Common case: N1<>N2
            //
            J:=M1;
            while J<=M2-1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    JP1 := J+1;
                    APVMove(@WORK[0], N1, N2, @A[JP1][0], N1, N2, CTEMP);
                    APVSub(@WORK[0], N1, N2, @A[J][0], N1, N2, STEMP);
                    APVMul(@A[J][0], N1, N2, CTEMP);
                    APVAdd(@A[J][0], N1, N2, @A[JP1][0], N1, N2, STEMP);
                    APVMove(@A[JP1][0], N1, N2, @WORK[0], N1, N2);
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // Special case: N1=N2
            //
            J:=M1;
            while J<=M2-1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    TEMP := A[J+1,N1];
                    A[J+1,N1] := CTEMP*TEMP-STEMP*A[J,N1];
                    A[J,N1] := STEMP*TEMP+CTEMP*A[J,N1];
                end;
                Inc(J);
            end;
        end;
    end
    else
    begin
        if N1<>N2 then
        begin
            
            //
            // Common case: N1<>N2
            //
            J:=M2-1;
            while J>=M1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    JP1 := J+1;
                    APVMove(@WORK[0], N1, N2, @A[JP1][0], N1, N2, CTEMP);
                    APVSub(@WORK[0], N1, N2, @A[J][0], N1, N2, STEMP);
                    APVMul(@A[J][0], N1, N2, CTEMP);
                    APVAdd(@A[J][0], N1, N2, @A[JP1][0], N1, N2, STEMP);
                    APVMove(@A[JP1][0], N1, N2, @WORK[0], N1, N2);
                end;
                Dec(J);
            end;
        end
        else
        begin
            
            //
            // Special case: N1=N2
            //
            J:=M2-1;
            while J>=M1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    TEMP := A[J+1,N1];
                    A[J+1,N1] := CTEMP*TEMP-STEMP*A[J,N1];
                    A[J,N1] := STEMP*TEMP+CTEMP*A[J,N1];
                end;
                Dec(J);
            end;
        end;
    end;
end;


(*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [M1..M2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************)
procedure ApplyRotationsFromTheRight(IsForward : Boolean;
     M1 : AlglibInteger;
     M2 : AlglibInteger;
     N1 : AlglibInteger;
     N2 : AlglibInteger;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
var
    J : AlglibInteger;
    JP1 : AlglibInteger;
    CTEMP : Double;
    STEMP : Double;
    TEMP : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Form A * P'
    //
    if IsForward then
    begin
        if M1<>M2 then
        begin
            
            //
            // Common case: M1<>M2
            //
            J:=N1;
            while J<=N2-1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    JP1 := J+1;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := CTEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := WORK[i_] - STEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := CTEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := A[i_,J] + STEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,JP1] := WORK[i_];
                    end;
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // Special case: M1=M2
            //
            J:=N1;
            while J<=N2-1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    TEMP := A[M1,J+1];
                    A[M1,J+1] := CTEMP*TEMP-STEMP*A[M1,J];
                    A[M1,J] := STEMP*TEMP+CTEMP*A[M1,J];
                end;
                Inc(J);
            end;
        end;
    end
    else
    begin
        if M1<>M2 then
        begin
            
            //
            // Common case: M1<>M2
            //
            J:=N2-1;
            while J>=N1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    JP1 := J+1;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := CTEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := WORK[i_] - STEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := CTEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := A[i_,J] + STEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,JP1] := WORK[i_];
                    end;
                end;
                Dec(J);
            end;
        end
        else
        begin
            
            //
            // Special case: M1=M2
            //
            J:=N2-1;
            while J>=N1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if AP_FP_Neq(CTEMP,1) or AP_FP_Neq(STEMP,0) then
                begin
                    TEMP := A[M1,J+1];
                    A[M1,J+1] := CTEMP*TEMP-STEMP*A[M1,J];
                    A[M1,J] := STEMP*TEMP+CTEMP*A[M1,J];
                end;
                Dec(J);
            end;
        end;
    end;
end;


(*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************)
procedure GenerateRotation(F : Double;
     G : Double;
     var CS : Double;
     var SN : Double;
     var R : Double);
var
    F1 : Double;
    G1 : Double;
begin
    if AP_FP_Eq(G,0) then
    begin
        CS := 1;
        SN := 0;
        R := F;
    end
    else
    begin
        if AP_FP_Eq(F,0) then
        begin
            CS := 0;
            SN := 1;
            R := G;
        end
        else
        begin
            F1 := F;
            G1 := G;
            if AP_FP_Greater(AbsReal(F1),AbsReal(G1)) then
            begin
                R := AbsReal(F1)*SQRT(1+AP_Sqr(G1/F1));
            end
            else
            begin
                R := AbsReal(G1)*SQRT(1+AP_Sqr(F1/G1));
            end;
            CS := F1/R;
            SN := G1/R;
            if AP_FP_Greater(ABSReal(F),ABSReal(G)) and AP_FP_Less(CS,0) then
            begin
                CS := -CS;
                SN := -SN;
                R := -R;
            end;
        end;
    end;
end;


end.