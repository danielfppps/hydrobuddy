{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright 2008 by Sergey Bochkanov (ALGLIB project).

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
unit bdss;
interface
uses Math, Sysutils, Ap, tsort, descriptivestatistics;

type
CVReport = record
    RelCLSError : Double;
    AvgCE : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
end;



procedure DSErrAllocate(NClasses : AlglibInteger; var Buf : TReal1DArray);
procedure DSErrAccumulate(var Buf : TReal1DArray;
     const Y : TReal1DArray;
     const DesiredY : TReal1DArray);
procedure DSErrFinish(var Buf : TReal1DArray);
procedure DSNormalize(var XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var Means : TReal1DArray;
     var Sigmas : TReal1DArray);
procedure DSNormalizeC(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var Means : TReal1DArray;
     var Sigmas : TReal1DArray);
function DSGetMeanMinDistance(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger):Double;
procedure DSTie(var A : TReal1DArray;
     N : AlglibInteger;
     var Ties : TInteger1DArray;
     var TieCount : AlglibInteger;
     var P1 : TInteger1DArray;
     var P2 : TInteger1DArray);
procedure DSTieFastI(var A : TReal1DArray;
     var B : TInteger1DArray;
     N : AlglibInteger;
     var Ties : TInteger1DArray;
     var TieCount : AlglibInteger);
procedure DSOptimalSplit2(A : TReal1DArray;
     C : TInteger1DArray;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var PAL : Double;
     var PBL : Double;
     var PAR : Double;
     var PBR : Double;
     var CVE : Double);
procedure DSOptimalSplit2Fast(var A : TReal1DArray;
     var C : TInteger1DArray;
     var TiesBuf : TInteger1DArray;
     var CntBuf : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     Alpha : Double;
     var Info : AlglibInteger;
     var Threshold : Double;
     var RMS : Double;
     var CVRMS : Double);
procedure DSSplitK(A : TReal1DArray;
     C : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     KMax : AlglibInteger;
     var Info : AlglibInteger;
     var Thresholds : TReal1DArray;
     var NI : AlglibInteger;
     var CVE : Double);
procedure DSOptimalSplitK(A : TReal1DArray;
     C : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     KMax : AlglibInteger;
     var Info : AlglibInteger;
     var Thresholds : TReal1DArray;
     var NI : AlglibInteger;
     var CVE : Double);

implementation

procedure DSKFoldSplit(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NClasses : AlglibInteger;
     FoldsCount : AlglibInteger;
     StratifiedSplits : Boolean;
     var Folds : TInteger1DArray);forward;
function XLNY(X : Double; Y : Double):Double;forward;
function GetCV(const Cnt : TInteger1DArray; NC : AlglibInteger):Double;forward;
procedure TieAddC(const C : TInteger1DArray;
     const Ties : TInteger1DArray;
     NTie : AlglibInteger;
     NC : AlglibInteger;
     var Cnt : TInteger1DArray);forward;
procedure TieSubC(const C : TInteger1DArray;
     const Ties : TInteger1DArray;
     NTie : AlglibInteger;
     NC : AlglibInteger;
     var Cnt : TInteger1DArray);forward;
procedure TieGetC(const C : TInteger1DArray;
     const Ties : TInteger1DArray;
     NTie : AlglibInteger;
     NC : AlglibInteger;
     var Cnt : TInteger1DArray);forward;


(*************************************************************************
This set of routines (DSErrAllocate, DSErrAccumulate, DSErrFinish)
calculates different error functions (classification error, cross-entropy,
rms, avg, avg.rel errors).

1. DSErrAllocate prepares buffer.
2. DSErrAccumulate accumulates individual errors:
    * Y contains predicted output (posterior probabilities for classification)
    * DesiredY contains desired output (class number for classification)
3. DSErrFinish outputs results:
   * Buf[0] contains relative classification error (zero for regression tasks)
   * Buf[1] contains avg. cross-entropy (zero for regression tasks)
   * Buf[2] contains rms error (regression, classification)
   * Buf[3] contains average error (regression, classification)
   * Buf[4] contains average relative error (regression, classification)
   
NOTES(1):
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.

NOTES(2):
    rms. avg, avg.rel errors for classification tasks are interpreted as
    errors in posterior probabilities with respect to probabilities given
    by training/test set.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************)
procedure DSErrAllocate(NClasses : AlglibInteger; var Buf : TReal1DArray);
begin
    SetLength(Buf, 7+1);
    Buf[0] := 0;
    Buf[1] := 0;
    Buf[2] := 0;
    Buf[3] := 0;
    Buf[4] := 0;
    Buf[5] := NClasses;
    Buf[6] := 0;
    Buf[7] := 0;
end;


(*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************)
procedure DSErrAccumulate(var Buf : TReal1DArray;
     const Y : TReal1DArray;
     const DesiredY : TReal1DArray);
var
    NClasses : AlglibInteger;
    NOut : AlglibInteger;
    Offs : AlglibInteger;
    MMax : AlglibInteger;
    RMax : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    EV : Double;
begin
    Offs := 5;
    NClasses := Round(Buf[Offs]);
    if NClasses>0 then
    begin
        
        //
        // Classification
        //
        RMax := Round(DesiredY[0]);
        MMax := 0;
        J:=1;
        while J<=NClasses-1 do
        begin
            if AP_FP_Greater(Y[J],Y[MMax]) then
            begin
                MMax := J;
            end;
            Inc(J);
        end;
        if MMax<>RMax then
        begin
            Buf[0] := Buf[0]+1;
        end;
        if AP_FP_Greater(Y[RMax],0) then
        begin
            Buf[1] := Buf[1]-Ln(Y[RMax]);
        end
        else
        begin
            Buf[1] := Buf[1]+Ln(MaxRealNumber);
        end;
        J:=0;
        while J<=NClasses-1 do
        begin
            V := Y[J];
            if J=RMax then
            begin
                EV := 1;
            end
            else
            begin
                EV := 0;
            end;
            Buf[2] := Buf[2]+AP_Sqr(V-EV);
            Buf[3] := Buf[3]+AbsReal(V-EV);
            if AP_FP_Neq(EV,0) then
            begin
                Buf[4] := Buf[4]+AbsReal((V-EV)/EV);
                Buf[Offs+2] := Buf[Offs+2]+1;
            end;
            Inc(J);
        end;
        Buf[Offs+1] := Buf[Offs+1]+1;
    end
    else
    begin
        
        //
        // Regression
        //
        NOut := -NClasses;
        RMax := 0;
        J:=1;
        while J<=NOut-1 do
        begin
            if AP_FP_Greater(DesiredY[J],DesiredY[RMax]) then
            begin
                RMax := J;
            end;
            Inc(J);
        end;
        MMax := 0;
        J:=1;
        while J<=NOut-1 do
        begin
            if AP_FP_Greater(Y[J],Y[MMax]) then
            begin
                MMax := J;
            end;
            Inc(J);
        end;
        if MMax<>RMax then
        begin
            Buf[0] := Buf[0]+1;
        end;
        J:=0;
        while J<=NOut-1 do
        begin
            V := Y[J];
            EV := DesiredY[J];
            Buf[2] := Buf[2]+AP_Sqr(V-EV);
            Buf[3] := Buf[3]+AbsReal(V-EV);
            if AP_FP_Neq(EV,0) then
            begin
                Buf[4] := Buf[4]+AbsReal((V-EV)/EV);
                Buf[Offs+2] := Buf[Offs+2]+1;
            end;
            Inc(J);
        end;
        Buf[Offs+1] := Buf[Offs+1]+1;
    end;
end;


(*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************)
procedure DSErrFinish(var Buf : TReal1DArray);
var
    NOut : AlglibInteger;
    Offs : AlglibInteger;
begin
    Offs := 5;
    NOut := AbsInt(Round(Buf[Offs]));
    if AP_FP_Neq(Buf[Offs+1],0) then
    begin
        Buf[0] := Buf[0]/Buf[Offs+1];
        Buf[1] := Buf[1]/Buf[Offs+1];
        Buf[2] := Sqrt(Buf[2]/(NOut*Buf[Offs+1]));
        Buf[3] := Buf[3]/(NOut*Buf[Offs+1]);
    end;
    if AP_FP_Neq(Buf[Offs+2],0) then
    begin
        Buf[4] := Buf[4]/Buf[Offs+2];
    end;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSNormalize(var XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var Means : TReal1DArray;
     var Sigmas : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    Tmp : TReal1DArray;
    Mean : Double;
    Variance : Double;
    Skewness : Double;
    Kurtosis : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Test parameters
    //
    if (NPoints<=0) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    
    //
    // Standartization
    //
    SetLength(Means, NVars-1+1);
    SetLength(Sigmas, NVars-1+1);
    SetLength(Tmp, NPoints-1+1);
    J:=0;
    while J<=NVars-1 do
    begin
        for i_ := 0 to NPoints-1 do
        begin
            Tmp[i_] := XY[i_,J];
        end;
        CalculateMoments(Tmp, NPoints, Mean, Variance, Skewness, Kurtosis);
        Means[J] := Mean;
        Sigmas[J] := Sqrt(Variance);
        if AP_FP_Eq(Sigmas[J],0) then
        begin
            Sigmas[J] := 1;
        end;
        I:=0;
        while I<=NPoints-1 do
        begin
            XY[I,J] := (XY[I,J]-Means[J])/Sigmas[J];
            Inc(I);
        end;
        Inc(J);
    end;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSNormalizeC(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var Means : TReal1DArray;
     var Sigmas : TReal1DArray);
var
    J : AlglibInteger;
    Tmp : TReal1DArray;
    Mean : Double;
    Variance : Double;
    Skewness : Double;
    Kurtosis : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Test parameters
    //
    if (NPoints<=0) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    
    //
    // Standartization
    //
    SetLength(Means, NVars-1+1);
    SetLength(Sigmas, NVars-1+1);
    SetLength(Tmp, NPoints-1+1);
    J:=0;
    while J<=NVars-1 do
    begin
        for i_ := 0 to NPoints-1 do
        begin
            Tmp[i_] := XY[i_,J];
        end;
        CalculateMoments(Tmp, NPoints, Mean, Variance, Skewness, Kurtosis);
        Means[J] := Mean;
        Sigmas[J] := Sqrt(Variance);
        if AP_FP_Eq(Sigmas[J],0) then
        begin
            Sigmas[J] := 1;
        end;
        Inc(J);
    end;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************)
function DSGetMeanMinDistance(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    Tmp : TReal1DArray;
    Tmp2 : TReal1DArray;
    V : Double;
begin
    
    //
    // Test parameters
    //
    if (NPoints<=0) or (NVars<1) then
    begin
        Result := 0;
        Exit;
    end;
    
    //
    // Process
    //
    SetLength(Tmp, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        Tmp[I] := MaxRealNumber;
        Inc(I);
    end;
    SetLength(Tmp2, NVars-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        J:=I+1;
        while J<=NPoints-1 do
        begin
            APVMove(@Tmp2[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
            APVSub(@Tmp2[0], 0, NVars-1, @XY[J][0], 0, NVars-1);
            V := APVDotProduct(@Tmp2[0], 0, NVars-1, @Tmp2[0], 0, NVars-1);
            V := Sqrt(V);
            Tmp[I] := Min(Tmp[I], V);
            Tmp[J] := Min(Tmp[J], V);
            Inc(J);
        end;
        Inc(I);
    end;
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        Result := Result+Tmp[I]/NPoints;
        Inc(I);
    end;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSTie(var A : TReal1DArray;
     N : AlglibInteger;
     var Ties : TInteger1DArray;
     var TieCount : AlglibInteger;
     var P1 : TInteger1DArray;
     var P2 : TInteger1DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    Tmp : TInteger1DArray;
begin
    
    //
    // Special case
    //
    if N<=0 then
    begin
        TieCount := 0;
        Exit;
    end;
    
    //
    // Sort A
    //
    TagSort(A, N, P1, P2);
    
    //
    // Process ties
    //
    TieCount := 1;
    I:=1;
    while I<=N-1 do
    begin
        if AP_FP_Neq(A[I],A[I-1]) then
        begin
            TieCount := TieCount+1;
        end;
        Inc(I);
    end;
    SetLength(Ties, TieCount+1);
    Ties[0] := 0;
    K := 1;
    I:=1;
    while I<=N-1 do
    begin
        if AP_FP_Neq(A[I],A[I-1]) then
        begin
            Ties[K] := I;
            K := K+1;
        end;
        Inc(I);
    end;
    Ties[TieCount] := N;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSTieFastI(var A : TReal1DArray;
     var B : TInteger1DArray;
     N : AlglibInteger;
     var Ties : TInteger1DArray;
     var TieCount : AlglibInteger);
var
    I : AlglibInteger;
    K : AlglibInteger;
    Tmp : TInteger1DArray;
begin
    
    //
    // Special case
    //
    if N<=0 then
    begin
        TieCount := 0;
        Exit;
    end;
    
    //
    // Sort A
    //
    TagSortFastI(A, B, N);
    
    //
    // Process ties
    //
    Ties[0] := 0;
    K := 1;
    I:=1;
    while I<=N-1 do
    begin
        if AP_FP_Neq(A[I],A[I-1]) then
        begin
            Ties[K] := I;
            K := K+1;
        end;
        Inc(I);
    end;
    Ties[K] := N;
    TieCount := K;
end;


(*************************************************************************
Optimal partition, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSOptimalSplit2(A : TReal1DArray;
     C : TInteger1DArray;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var PAL : Double;
     var PBL : Double;
     var PAR : Double;
     var PBR : Double;
     var CVE : Double);
var
    I : AlglibInteger;
    T : AlglibInteger;
    S : Double;
    Ties : TInteger1DArray;
    TieCount : AlglibInteger;
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
    K : AlglibInteger;
    KOptimal : AlglibInteger;
    PAK : Double;
    PBK : Double;
    CVOptimal : Double;
    CV : Double;
begin
    A := DynamicArrayCopy(A);
    C := DynamicArrayCopy(C);
    
    //
    // Test for errors in inputs
    //
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if (C[I]<>0) and (C[I]<>1) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Tie
    //
    DSTie(A, N, Ties, TieCount, P1, P2);
    I:=0;
    while I<=N-1 do
    begin
        if P2[I]<>I then
        begin
            T := C[I];
            C[I] := C[P2[I]];
            C[P2[I]] := T;
        end;
        Inc(I);
    end;
    
    //
    // Special case: number of ties is 1.
    //
    // NOTE: we assume that P[i,j] equals to 0 or 1,
    //       intermediate values are not allowed.
    //
    if TieCount=1 then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // General case, number of ties > 1
    //
    // NOTE: we assume that P[i,j] equals to 0 or 1,
    //       intermediate values are not allowed.
    //
    PAL := 0;
    PBL := 0;
    PAR := 0;
    PBR := 0;
    I:=0;
    while I<=N-1 do
    begin
        if C[I]=0 then
        begin
            PAR := PAR+1;
        end;
        if C[I]=1 then
        begin
            PBR := PBR+1;
        end;
        Inc(I);
    end;
    KOptimal := -1;
    CVOptimal := MaxRealNumber;
    K:=0;
    while K<=TieCount-2 do
    begin
        
        //
        // first, obtain information about K-th tie which is
        // moved from R-part to L-part
        //
        PAK := 0;
        PBK := 0;
        I:=Ties[K];
        while I<=Ties[K+1]-1 do
        begin
            if C[I]=0 then
            begin
                PAK := PAK+1;
            end;
            if C[I]=1 then
            begin
                PBK := PBK+1;
            end;
            Inc(I);
        end;
        
        //
        // Calculate cross-validation CE
        //
        CV := 0;
        CV := CV-XLNY(PAL+PAK, (PAL+PAK)/(PAL+PAK+PBL+PBK+1));
        CV := CV-XLNY(PBL+PBK, (PBL+PBK)/(PAL+PAK+1+PBL+PBK));
        CV := CV-XLNY(PAR-PAK, (PAR-PAK)/(PAR-PAK+PBR-PBK+1));
        CV := CV-XLNY(PBR-PBK, (PBR-PBK)/(PAR-PAK+1+PBR-PBK));
        
        //
        // Compare with best
        //
        if AP_FP_Less(CV,CVOptimal) then
        begin
            CVOptimal := CV;
            KOptimal := K;
        end;
        
        //
        // update
        //
        PAL := PAL+PAK;
        PBL := PBL+PBK;
        PAR := PAR-PAK;
        PBR := PBR-PBK;
        Inc(K);
    end;
    CVE := CVOptimal;
    Threshold := Double(0.5)*(A[Ties[KOptimal]]+A[Ties[KOptimal+1]]);
    PAL := 0;
    PBL := 0;
    PAR := 0;
    PBR := 0;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Less(A[I],Threshold) then
        begin
            if C[I]=0 then
            begin
                PAL := PAL+1;
            end
            else
            begin
                PBL := PBL+1;
            end;
        end
        else
        begin
            if C[I]=0 then
            begin
                PAR := PAR+1;
            end
            else
            begin
                PBR := PBR+1;
            end;
        end;
        Inc(I);
    end;
    S := PAL+PBL;
    PAL := PAL/S;
    PBL := PBL/S;
    S := PAR+PBR;
    PAR := PAR/S;
    PBR := PBR/S;
end;


(*************************************************************************
Optimal partition, internal subroutine. Fast version.

Accepts:
    A       array[0..N-1]       array of attributes     array[0..N-1]
    C       array[0..N-1]       array of class labels
    TiesBuf array[0..N]         temporaries (ties)
    CntBuf  array[0..2*NC-1]    temporaries (counts)
    Alpha                       centering factor (0<=alpha<=1, recommended value - 0.05)
    
Output:
    Info    error code (">0"=OK, "<0"=bad)
    RMS     training set RMS error
    CVRMS   leave-one-out RMS error
    
Note:
    content of all arrays is changed by subroutine

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSOptimalSplit2Fast(var A : TReal1DArray;
     var C : TInteger1DArray;
     var TiesBuf : TInteger1DArray;
     var CntBuf : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     Alpha : Double;
     var Info : AlglibInteger;
     var Threshold : Double;
     var RMS : Double;
     var CVRMS : Double);
var
    I : AlglibInteger;
    K : AlglibInteger;
    CL : AlglibInteger;
    TieCount : AlglibInteger;
    CBest : Double;
    CC : Double;
    KOptimal : AlglibInteger;
    SL : AlglibInteger;
    SR : AlglibInteger;
    V : Double;
    W : Double;
    X : Double;
begin
    
    //
    // Test for errors in inputs
    //
    if (N<=0) or (NC<2) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if (C[I]<0) or (C[I]>=NC) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Tie
    //
    DSTieFastI(A, C, N, TiesBuf, TieCount);
    
    //
    // Special case: number of ties is 1.
    //
    if TieCount=1 then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // General case, number of ties > 1
    //
    I:=0;
    while I<=2*NC-1 do
    begin
        CntBuf[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        CntBuf[NC+C[I]] := CntBuf[NC+C[I]]+1;
        Inc(I);
    end;
    KOptimal := -1;
    Threshold := A[N-1];
    CBest := MaxRealNumber;
    SL := 0;
    SR := N;
    K:=0;
    while K<=TieCount-2 do
    begin
        
        //
        // first, move Kth tie from right to left
        //
        I:=TiesBuf[K];
        while I<=TiesBuf[K+1]-1 do
        begin
            CL := C[I];
            CntBuf[CL] := CntBuf[CL]+1;
            CntBuf[NC+CL] := CntBuf[NC+CL]-1;
            Inc(I);
        end;
        SL := SL+(TiesBuf[K+1]-TiesBuf[K]);
        SR := SR-(TiesBuf[K+1]-TiesBuf[K]);
        
        //
        // Calculate RMS error
        //
        V := 0;
        I:=0;
        while I<=NC-1 do
        begin
            W := CntBuf[I];
            V := V+W*AP_Sqr(W/SL-1);
            V := V+(SL-W)*AP_Sqr(W/SL);
            W := CntBuf[NC+I];
            V := V+W*AP_Sqr(W/SR-1);
            V := V+(SR-W)*AP_Sqr(W/SR);
            Inc(I);
        end;
        V := Sqrt(V/(NC*N));
        
        //
        // Compare with best
        //
        X := AP_Double(2*SL)/(SL+SR)-1;
        CC := V*(1-Alpha+Alpha*AP_Sqr(X));
        if AP_FP_Less(CC,CBest) then
        begin
            
            //
            // store split
            //
            RMS := V;
            KOptimal := K;
            CBest := CC;
            
            //
            // calculate CVRMS error
            //
            CVRMS := 0;
            I:=0;
            while I<=NC-1 do
            begin
                if SL>1 then
                begin
                    W := CntBuf[I];
                    CVRMS := CVRMS+W*AP_Sqr((W-1)/(SL-1)-1);
                    CVRMS := CVRMS+(SL-W)*AP_Sqr(W/(SL-1));
                end
                else
                begin
                    W := CntBuf[I];
                    CVRMS := CVRMS+W*AP_Sqr(AP_Double(1)/NC-1);
                    CVRMS := CVRMS+(SL-W)*AP_Sqr(AP_Double(1)/NC);
                end;
                if SR>1 then
                begin
                    W := CntBuf[NC+I];
                    CVRMS := CVRMS+W*AP_Sqr((W-1)/(SR-1)-1);
                    CVRMS := CVRMS+(SR-W)*AP_Sqr(W/(SR-1));
                end
                else
                begin
                    W := CntBuf[NC+I];
                    CVRMS := CVRMS+W*AP_Sqr(AP_Double(1)/NC-1);
                    CVRMS := CVRMS+(SR-W)*AP_Sqr(AP_Double(1)/NC);
                end;
                Inc(I);
            end;
            CVRMS := Sqrt(CVRMS/(NC*N));
        end;
        Inc(K);
    end;
    
    //
    // Calculate threshold.
    // Code is a bit complicated because there can be such
    // numbers that 0.5(A+B) equals to A or B (if A-B=epsilon)
    //
    Threshold := Double(0.5)*(A[TiesBuf[KOptimal]]+A[TiesBuf[KOptimal+1]]);
    if AP_FP_Less_Eq(Threshold,A[TiesBuf[KOptimal]]) then
    begin
        Threshold := A[TiesBuf[KOptimal+1]];
    end;
end;


(*************************************************************************
Automatic non-optimal discretization, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSSplitK(A : TReal1DArray;
     C : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     KMax : AlglibInteger;
     var Info : AlglibInteger;
     var Thresholds : TReal1DArray;
     var NI : AlglibInteger;
     var CVE : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    K : AlglibInteger;
    Ties : TInteger1DArray;
    TieCount : AlglibInteger;
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
    Cnt : TInteger1DArray;
    V2 : Double;
    BestK : AlglibInteger;
    BestCVE : Double;
    BestSizes : TInteger1DArray;
    CurCVE : Double;
    CurSizes : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    C := DynamicArrayCopy(C);
    
    //
    // Test for errors in inputs
    //
    if (N<=0) or (NC<2) or (KMax<2) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if (C[I]<0) or (C[I]>=NC) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Tie
    //
    DSTie(A, N, Ties, TieCount, P1, P2);
    I:=0;
    while I<=N-1 do
    begin
        if P2[I]<>I then
        begin
            K := C[I];
            C[I] := C[P2[I]];
            C[P2[I]] := K;
        end;
        Inc(I);
    end;
    
    //
    // Special cases
    //
    if TieCount=1 then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // General case:
    // 0. allocate arrays
    //
    KMax := Min(KMax, TieCount);
    SetLength(BestSizes, KMax-1+1);
    SetLength(CurSizes, KMax-1+1);
    SetLength(Cnt, NC-1+1);
    
    //
    // General case:
    // 1. prepare "weak" solution (two subintervals, divided at median)
    //
    V2 := MaxRealNumber;
    J := -1;
    I:=1;
    while I<=TieCount-1 do
    begin
        if AP_FP_Less(AbsReal(Ties[I]-Double(0.5)*(N-1)),V2) then
        begin
            V2 := AbsReal(Ties[I]-Double(0.5)*N);
            J := I;
        end;
        Inc(I);
    end;
    Assert(J>0, 'DSSplitK: internal error #1!');
    BestK := 2;
    BestSizes[0] := Ties[J];
    BestSizes[1] := N-J;
    BestCVE := 0;
    I:=0;
    while I<=NC-1 do
    begin
        Cnt[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=J-1 do
    begin
        TieAddC(C, Ties, I, NC, Cnt);
        Inc(I);
    end;
    BestCVE := BestCVE+GetCV(Cnt, NC);
    I:=0;
    while I<=NC-1 do
    begin
        Cnt[I] := 0;
        Inc(I);
    end;
    I:=J;
    while I<=TieCount-1 do
    begin
        TieAddC(C, Ties, I, NC, Cnt);
        Inc(I);
    end;
    BestCVE := BestCVE+GetCV(Cnt, NC);
    
    //
    // General case:
    // 2. Use greedy algorithm to find sub-optimal split in O(KMax*N) time
    //
    K:=2;
    while K<=KMax do
    begin
        
        //
        // Prepare greedy K-interval split
        //
        I:=0;
        while I<=K-1 do
        begin
            CurSizes[I] := 0;
            Inc(I);
        end;
        I := 0;
        J := 0;
        while (J<=TieCount-1) and (I<=K-1) do
        begin
            
            //
            // Rule: I-th bin is empty, fill it
            //
            if CurSizes[I]=0 then
            begin
                CurSizes[I] := Ties[J+1]-Ties[J];
                J := J+1;
                Continue;
            end;
            
            //
            // Rule: (K-1-I) bins left, (K-1-I) ties left (1 tie per bin); next bin
            //
            if TieCount-J=K-1-I then
            begin
                I := I+1;
                Continue;
            end;
            
            //
            // Rule: last bin, always place in current
            //
            if I=K-1 then
            begin
                CurSizes[I] := CurSizes[I]+Ties[J+1]-Ties[J];
                J := J+1;
                Continue;
            end;
            
            //
            // Place J-th tie in I-th bin, or leave for I+1-th bin.
            //
            if AP_FP_Less(AbsReal(CurSizes[I]+Ties[J+1]-Ties[J]-AP_Double(N)/K),AbsReal(CurSizes[I]-AP_Double(N)/K)) then
            begin
                CurSizes[I] := CurSizes[I]+Ties[J+1]-Ties[J];
                J := J+1;
            end
            else
            begin
                I := I+1;
            end;
        end;
        Assert((CurSizes[K-1]<>0) and (J=TieCount), 'DSSplitK: internal error #1');
        
        //
        // Calculate CVE
        //
        CurCVE := 0;
        J := 0;
        I:=0;
        while I<=K-1 do
        begin
            J1:=0;
            while J1<=NC-1 do
            begin
                Cnt[J1] := 0;
                Inc(J1);
            end;
            J1:=J;
            while J1<=J+CurSizes[I]-1 do
            begin
                Cnt[C[J1]] := Cnt[C[J1]]+1;
                Inc(J1);
            end;
            CurCVE := CurCVE+GetCV(Cnt, NC);
            J := J+CurSizes[I];
            Inc(I);
        end;
        
        //
        // Choose best variant
        //
        if AP_FP_Less(CurCVE,BestCVE) then
        begin
            I:=0;
            while I<=K-1 do
            begin
                BestSizes[I] := CurSizes[I];
                Inc(I);
            end;
            BestCVE := CurCVE;
            BestK := K;
        end;
        Inc(K);
    end;
    
    //
    // Transform from sizes to thresholds
    //
    CVE := BestCVE;
    NI := BestK;
    SetLength(Thresholds, NI-2+1);
    J := BestSizes[0];
    I:=1;
    while I<=BestK-1 do
    begin
        Thresholds[I-1] := Double(0.5)*(A[J-1]+A[J]);
        J := J+BestSizes[I];
        Inc(I);
    end;
end;


(*************************************************************************
Automatic optimal discretization, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure DSOptimalSplitK(A : TReal1DArray;
     C : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     KMax : AlglibInteger;
     var Info : AlglibInteger;
     var Thresholds : TReal1DArray;
     var NI : AlglibInteger;
     var CVE : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    S : AlglibInteger;
    JL : AlglibInteger;
    JR : AlglibInteger;
    V2 : Double;
    Ties : TInteger1DArray;
    TieCount : AlglibInteger;
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
    CVTemp : Double;
    Cnt : TInteger1DArray;
    Cnt2 : TInteger1DArray;
    CV : TReal2DArray;
    Splits : TInteger2DArray;
    K : AlglibInteger;
    KOptimal : AlglibInteger;
    CVOptimal : Double;
begin
    A := DynamicArrayCopy(A);
    C := DynamicArrayCopy(C);
    
    //
    // Test for errors in inputs
    //
    if (N<=0) or (NC<2) or (KMax<2) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if (C[I]<0) or (C[I]>=NC) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Tie
    //
    DSTie(A, N, Ties, TieCount, P1, P2);
    I:=0;
    while I<=N-1 do
    begin
        if P2[I]<>I then
        begin
            K := C[I];
            C[I] := C[P2[I]];
            C[P2[I]] := K;
        end;
        Inc(I);
    end;
    
    //
    // Special cases
    //
    if TieCount=1 then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // General case
    // Use dynamic programming to find best split in O(KMax*NC*TieCount^2) time
    //
    KMax := Min(KMax, TieCount);
    SetLength(CV, KMax-1+1, TieCount-1+1);
    SetLength(Splits, KMax-1+1, TieCount-1+1);
    SetLength(Cnt, NC-1+1);
    SetLength(Cnt2, NC-1+1);
    J:=0;
    while J<=NC-1 do
    begin
        Cnt[J] := 0;
        Inc(J);
    end;
    J:=0;
    while J<=TieCount-1 do
    begin
        TieAddC(C, Ties, J, NC, Cnt);
        Splits[0,J] := 0;
        CV[0,J] := GetCV(Cnt, NC);
        Inc(J);
    end;
    K:=1;
    while K<=KMax-1 do
    begin
        J:=0;
        while J<=NC-1 do
        begin
            Cnt[J] := 0;
            Inc(J);
        end;
        
        //
        // Subtask size J in [K..TieCount-1]:
        // optimal K-splitting on ties from 0-th to J-th.
        //
        J:=K;
        while J<=TieCount-1 do
        begin
            
            //
            // Update Cnt - let it contain classes of ties from K-th to J-th
            //
            TieAddC(C, Ties, J, NC, Cnt);
            
            //
            // Search for optimal split point S in [K..J]
            //
            I:=0;
            while I<=NC-1 do
            begin
                Cnt2[I] := Cnt[I];
                Inc(I);
            end;
            CV[K,J] := CV[K-1,J-1]+GetCV(Cnt2, NC);
            Splits[K,J] := J;
            S:=K+1;
            while S<=J do
            begin
                
                //
                // Update Cnt2 - let it contain classes of ties from S-th to J-th
                //
                TieSubC(C, Ties, S-1, NC, Cnt2);
                
                //
                // Calculate CVE
                //
                CVTemp := CV[K-1,S-1]+GetCV(Cnt2, NC);
                if AP_FP_Less(CVTemp,CV[K,J]) then
                begin
                    CV[K,J] := CVTemp;
                    Splits[K,J] := S;
                end;
                Inc(S);
            end;
            Inc(J);
        end;
        Inc(K);
    end;
    
    //
    // Choose best partition, output result
    //
    KOptimal := -1;
    CVOptimal := MaxRealNumber;
    K:=0;
    while K<=KMax-1 do
    begin
        if AP_FP_Less(CV[K,TieCount-1],CVOptimal) then
        begin
            CVOptimal := CV[K,TieCount-1];
            KOptimal := K;
        end;
        Inc(K);
    end;
    Assert(KOptimal>=0, 'DSOptimalSplitK: internal error #1!');
    if KOptimal=0 then
    begin
        
        //
        // Special case: best partition is one big interval.
        // Even 2-partition is not better.
        // This is possible when dealing with "weak" predictor variables.
        //
        // Make binary split as close to the median as possible.
        //
        V2 := MaxRealNumber;
        J := -1;
        I:=1;
        while I<=TieCount-1 do
        begin
            if AP_FP_Less(AbsReal(Ties[I]-Double(0.5)*(N-1)),V2) then
            begin
                V2 := AbsReal(Ties[I]-Double(0.5)*(N-1));
                J := I;
            end;
            Inc(I);
        end;
        Assert(J>0, 'DSOptimalSplitK: internal error #2!');
        SetLength(Thresholds, 0+1);
        Thresholds[0] := Double(0.5)*(A[Ties[J-1]]+A[Ties[J]]);
        NI := 2;
        CVE := 0;
        I:=0;
        while I<=NC-1 do
        begin
            Cnt[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=J-1 do
        begin
            TieAddC(C, Ties, I, NC, Cnt);
            Inc(I);
        end;
        CVE := CVE+GetCV(Cnt, NC);
        I:=0;
        while I<=NC-1 do
        begin
            Cnt[I] := 0;
            Inc(I);
        end;
        I:=J;
        while I<=TieCount-1 do
        begin
            TieAddC(C, Ties, I, NC, Cnt);
            Inc(I);
        end;
        CVE := CVE+GetCV(Cnt, NC);
    end
    else
    begin
        
        //
        // General case: 2 or more intervals
        //
        SetLength(Thresholds, KOptimal-1+1);
        NI := KOptimal+1;
        CVE := CV[KOptimal,TieCount-1];
        JL := Splits[KOptimal,TieCount-1];
        JR := TieCount-1;
        K:=KOptimal;
        while K>=1 do
        begin
            Thresholds[K-1] := Double(0.5)*(A[Ties[JL-1]]+A[Ties[JL]]);
            JR := JL-1;
            JL := Splits[K-1,JL-1];
            Dec(K);
        end;
    end;
end;


(*************************************************************************
Subroutine prepares K-fold split of the training set.

NOTES:
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************)
procedure DSKFoldSplit(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NClasses : AlglibInteger;
     FoldsCount : AlglibInteger;
     StratifiedSplits : Boolean;
     var Folds : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
begin
    
    //
    // test parameters
    //
    Assert(NPoints>0, 'DSKFoldSplit: wrong NPoints!');
    Assert((NClasses>1) or (NClasses<0), 'DSKFoldSplit: wrong NClasses!');
    Assert((FoldsCount>=2) and (FoldsCount<=NPoints), 'DSKFoldSplit: wrong FoldsCount!');
    Assert( not StratifiedSplits, 'DSKFoldSplit: stratified splits are not supported!');
    
    //
    // Folds
    //
    SetLength(Folds, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        Folds[I] := I*FoldsCount div NPoints;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-2 do
    begin
        J := I+RandomInteger(NPoints-I);
        if J<>I then
        begin
            K := Folds[I];
            Folds[I] := Folds[J];
            Folds[J] := K;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal function
*************************************************************************)
function XLNY(X : Double; Y : Double):Double;
begin
    if AP_FP_Eq(X,0) then
    begin
        Result := 0;
    end
    else
    begin
        Result := X*Ln(Y);
    end;
end;


(*************************************************************************
Internal function,
returns number of samples of class I in Cnt[I]
*************************************************************************)
function GetCV(const Cnt : TInteger1DArray; NC : AlglibInteger):Double;
var
    I : AlglibInteger;
    S : Double;
begin
    S := 0;
    I:=0;
    while I<=NC-1 do
    begin
        S := S+Cnt[I];
        Inc(I);
    end;
    Result := 0;
    I:=0;
    while I<=NC-1 do
    begin
        Result := Result-XLNY(Cnt[I], Cnt[I]/(S+NC-1));
        Inc(I);
    end;
end;


(*************************************************************************
Internal function, adds number of samples of class I in tie NTie to Cnt[I]
*************************************************************************)
procedure TieAddC(const C : TInteger1DArray;
     const Ties : TInteger1DArray;
     NTie : AlglibInteger;
     NC : AlglibInteger;
     var Cnt : TInteger1DArray);
var
    I : AlglibInteger;
begin
    I:=Ties[NTie];
    while I<=Ties[NTie+1]-1 do
    begin
        Cnt[C[I]] := Cnt[C[I]]+1;
        Inc(I);
    end;
end;


(*************************************************************************
Internal function, subtracts number of samples of class I in tie NTie to Cnt[I]
*************************************************************************)
procedure TieSubC(const C : TInteger1DArray;
     const Ties : TInteger1DArray;
     NTie : AlglibInteger;
     NC : AlglibInteger;
     var Cnt : TInteger1DArray);
var
    I : AlglibInteger;
begin
    I:=Ties[NTie];
    while I<=Ties[NTie+1]-1 do
    begin
        Cnt[C[I]] := Cnt[C[I]]-1;
        Inc(I);
    end;
end;


(*************************************************************************
Internal function,
returns number of samples of class I in Cnt[I]
*************************************************************************)
procedure TieGetC(const C : TInteger1DArray;
     const Ties : TInteger1DArray;
     NTie : AlglibInteger;
     NC : AlglibInteger;
     var Cnt : TInteger1DArray);
var
    I : AlglibInteger;
begin
    I:=0;
    while I<=NC-1 do
    begin
        Cnt[I] := 0;
        Inc(I);
    end;
    I:=Ties[NTie];
    while I<=Ties[NTie+1]-1 do
    begin
        Cnt[C[I]] := Cnt[C[I]]+1;
        Inc(I);
    end;
end;


end.