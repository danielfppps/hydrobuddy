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
unit dforest;
interface
uses Math, Sysutils, Ap, tsort, descriptivestatistics, bdss;

type
DecisionForest = record
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    NTrees : AlglibInteger;
    BufSize : AlglibInteger;
    Trees : TReal1DArray;
end;


DFReport = record
    RelClsError : Double;
    AvgCE : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    OOBRelClsError : Double;
    OOBAvgCE : Double;
    OOBRMSError : Double;
    OOBAvgError : Double;
    OOBAvgRelError : Double;
end;


DFInternalBuffers = record
    TreeBuf : TReal1DArray;
    IdxBuf : TInteger1DArray;
    TmpBufR : TReal1DArray;
    TmpBufR2 : TReal1DArray;
    TmpBufI : TInteger1DArray;
    ClassIBuf : TInteger1DArray;
    VarPool : TInteger1DArray;
    EVSBin : TBoolean1DArray;
    EVSSplits : TReal1DArray;
end;



procedure DFBuildRandomDecisionForest(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NTrees : AlglibInteger;
     R : Double;
     var Info : AlglibInteger;
     var DF : DecisionForest;
     var Rep : DFReport);
procedure DFBuildInternal(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NTrees : AlglibInteger;
     SampleSize : AlglibInteger;
     NFeatures : AlglibInteger;
     Flags : AlglibInteger;
     var Info : AlglibInteger;
     var DF : DecisionForest;
     var Rep : DFReport);
procedure DFProcess(const DF : DecisionForest;
     const X : TReal1DArray;
     var Y : TReal1DArray);
function DFRelClsError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function DFAvgCE(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function DFRMSError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function DFAvgError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function DFAvgRelError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
procedure DFCopy(const DF1 : DecisionForest; var DF2 : DecisionForest);
procedure DFSerialize(const DF : DecisionForest;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
procedure DFUnserialize(const RA : TReal1DArray; var DF : DecisionForest);

implementation

const
    DFVNum = 8;
    InnerNodeWidth = 3;
    LeafNodeWidth = 2;
    DFUseStrongSplits = 1;
    DFUseEVS = 2;

function DFClsError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):AlglibInteger;forward;
procedure DFProcessInternal(const DF : DecisionForest;
     Offs : AlglibInteger;
     const X : TReal1DArray;
     var Y : TReal1DArray);forward;
procedure DFBuildTree(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NFeatures : AlglibInteger;
     NVarsInPool : AlglibInteger;
     Flags : AlglibInteger;
     var Bufs : DFInternalBuffers);forward;
procedure DFBuildTreeRec(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NFeatures : AlglibInteger;
     NVarsInPool : AlglibInteger;
     Flags : AlglibInteger;
     var NumProcessed : AlglibInteger;
     Idx1 : AlglibInteger;
     Idx2 : AlglibInteger;
     var Bufs : DFInternalBuffers);forward;
procedure DFWeakSplitI(var X : TReal1DArray;
     var Y : TInteger1DArray;
     N : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var E : Double);forward;
procedure DFSplitC(var X : TReal1DArray;
     var C : TInteger1DArray;
     var CntBuf : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     Flags : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var E : Double);forward;
procedure DFSplitR(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger;
     Flags : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var E : Double);forward;


(*************************************************************************
This subroutine builds random decision forest.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure DFBuildRandomDecisionForest(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NTrees : AlglibInteger;
     R : Double;
     var Info : AlglibInteger;
     var DF : DecisionForest;
     var Rep : DFReport);
var
    SampleSize : AlglibInteger;
begin
    if AP_FP_Less_Eq(R,0) or AP_FP_Greater(R,1) then
    begin
        Info := -1;
        Exit;
    end;
    SampleSize := Max(Round(R*NPoints), 1);
    DFBuildInternal(XY, NPoints, NVars, NClasses, NTrees, SampleSize, Max(NVars div 2, 1), DFUseStrongSplits+DFUseEVS, Info, DF, Rep);
end;


procedure DFBuildInternal(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NTrees : AlglibInteger;
     SampleSize : AlglibInteger;
     NFeatures : AlglibInteger;
     Flags : AlglibInteger;
     var Info : AlglibInteger;
     var DF : DecisionForest;
     var Rep : DFReport);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TmpI : AlglibInteger;
    LastTreeOffs : AlglibInteger;
    Offs : AlglibInteger;
    OOBOffs : AlglibInteger;
    TreeSize : AlglibInteger;
    NVarsInPool : AlglibInteger;
    UseEVS : Boolean;
    Bufs : DFInternalBuffers;
    PermBuf : TInteger1DArray;
    OOBBuf : TReal1DArray;
    OOBCntBuf : TInteger1DArray;
    XYS : TReal2DArray;
    X : TReal1DArray;
    Y : TReal1DArray;
    OOBCnt : AlglibInteger;
    OOBRelCnt : AlglibInteger;
    V : Double;
    VMin : Double;
    VMax : Double;
    BFlag : Boolean;
begin
    
    //
    // Test for inputs
    //
    if (NPoints<1) or (SampleSize<1) or (SampleSize>NPoints) or (NVars<1) or (NClasses<1) or (NTrees<1) or (NFeatures<1) then
    begin
        Info := -1;
        Exit;
    end;
    if NClasses>1 then
    begin
        I:=0;
        while I<=NPoints-1 do
        begin
            if (Round(XY[I,NVars])<0) or (Round(XY[I,NVars])>=NClasses) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
    end;
    Info := 1;
    
    //
    // Flags
    //
    UseEVS := Flags div DFUseEVS mod 2<>0;
    
    //
    // Allocate data, prepare header
    //
    TreeSize := 1+InnerNodeWidth*(SampleSize-1)+LeafNodeWidth*SampleSize;
    SetLength(PermBuf, NPoints-1+1);
    SetLength(Bufs.TreeBuf, TreeSize-1+1);
    SetLength(Bufs.IdxBuf, NPoints-1+1);
    SetLength(Bufs.TmpBufR, NPoints-1+1);
    SetLength(Bufs.TmpBufR2, NPoints-1+1);
    SetLength(Bufs.TmpBufI, NPoints-1+1);
    SetLength(Bufs.VarPool, NVars-1+1);
    SetLength(Bufs.EVSBin, NVars-1+1);
    SetLength(Bufs.EVSSplits, NVars-1+1);
    SetLength(Bufs.ClassIBuf, 2*NClasses-1+1);
    SetLength(OOBBuf, NClasses*NPoints-1+1);
    SetLength(OOBCntBuf, NPoints-1+1);
    SetLength(DF.Trees, NTrees*TreeSize-1+1);
    SetLength(XYS, SampleSize-1+1, NVars+1);
    SetLength(X, NVars-1+1);
    SetLength(Y, NClasses-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        PermBuf[I] := I;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints*NClasses-1 do
    begin
        OOBBuf[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        OOBCntBuf[I] := 0;
        Inc(I);
    end;
    
    //
    // Prepare variable pool and EVS (extended variable selection/splitting) buffers
    // (whether EVS is turned on or not):
    // 1. detect binary variables and pre-calculate splits for them
    // 2. detect variables with non-distinct values and exclude them from pool
    //
    I:=0;
    while I<=NVars-1 do
    begin
        Bufs.VarPool[I] := I;
        Inc(I);
    end;
    NVarsInPool := NVars;
    if UseEVS then
    begin
        J:=0;
        while J<=NVars-1 do
        begin
            VMin := XY[0,J];
            VMax := VMin;
            I:=0;
            while I<=NPoints-1 do
            begin
                V := XY[I,J];
                VMin := Min(VMin, V);
                VMax := Max(VMax, V);
                Inc(I);
            end;
            if AP_FP_Eq(VMin,VMax) then
            begin
                
                //
                // exclude variable from pool
                //
                Bufs.VarPool[J] := Bufs.VarPool[NVarsInPool-1];
                Bufs.VarPool[NVarsInPool-1] := -1;
                NVarsInPool := NVarsInPool-1;
                Inc(J);
                Continue;
            end;
            BFlag := False;
            I:=0;
            while I<=NPoints-1 do
            begin
                V := XY[I,J];
                if AP_FP_Neq(V,VMin) and AP_FP_Neq(V,VMax) then
                begin
                    BFlag := True;
                    Break;
                end;
                Inc(I);
            end;
            if BFlag then
            begin
                
                //
                // non-binary variable
                //
                Bufs.EVSBin[J] := False;
            end
            else
            begin
                
                //
                // Prepare
                //
                Bufs.EVSBin[J] := True;
                Bufs.EVSSplits[J] := Double(0.5)*(VMin+VMax);
                if AP_FP_Less_Eq(Bufs.EVSSplits[J],VMin) then
                begin
                    Bufs.EVSSplits[J] := VMax;
                end;
            end;
            Inc(J);
        end;
    end;
    
    //
    // RANDOM FOREST FORMAT
    // W[0]         -   size of array
    // W[1]         -   version number
    // W[2]         -   NVars
    // W[3]         -   NClasses (1 for regression)
    // W[4]         -   NTrees
    // W[5]         -   trees offset
    //
    //
    // TREE FORMAT
    // W[Offs]      -   size of sub-array
    //     node info:
    // W[K+0]       -   variable number        (-1 for leaf mode)
    // W[K+1]       -   threshold              (class/value for leaf node)
    // W[K+2]       -   ">=" branch index      (absent for leaf node)
    //
    //
    DF.NVars := NVars;
    DF.NClasses := NClasses;
    DF.NTrees := NTrees;
    
    //
    // Build forest
    //
    Offs := 0;
    I:=0;
    while I<=NTrees-1 do
    begin
        
        //
        // Prepare sample
        //
        K:=0;
        while K<=SampleSize-1 do
        begin
            J := K+RandomInteger(NPoints-K);
            TmpI := PermBuf[K];
            PermBuf[K] := PermBuf[J];
            PermBuf[J] := TmpI;
            J := PermBuf[K];
            APVMove(@XYS[K][0], 0, NVars, @XY[J][0], 0, NVars);
            Inc(K);
        end;
        
        //
        // build tree, copy
        //
        DFBuildTree(XYS, SampleSize, NVars, NClasses, NFeatures, NVarsInPool, Flags, Bufs);
        J := Round(Bufs.TreeBuf[0]);
        APVMove(@DF.Trees[0], Offs, Offs+J-1, @Bufs.TreeBuf[0], 0, J-1);
        LastTreeOffs := Offs;
        Offs := Offs+J;
        
        //
        // OOB estimates
        //
        K:=SampleSize;
        while K<=NPoints-1 do
        begin
            J:=0;
            while J<=NClasses-1 do
            begin
                Y[J] := 0;
                Inc(J);
            end;
            J := PermBuf[K];
            APVMove(@X[0], 0, NVars-1, @XY[J][0], 0, NVars-1);
            DFProcessInternal(DF, LastTreeOffs, X, Y);
            APVAdd(@OOBBuf[0], J*NClasses, (J+1)*NClasses-1, @Y[0], 0, NClasses-1);
            OOBCntBuf[J] := OOBCntBuf[J]+1;
            Inc(K);
        end;
        Inc(I);
    end;
    DF.BufSize := Offs;
    
    //
    // Normalize OOB results
    //
    I:=0;
    while I<=NPoints-1 do
    begin
        if OOBCntBuf[I]<>0 then
        begin
            V := AP_Double(1)/OOBCntBuf[I];
            APVMul(@OOBBuf[0], I*NClasses, I*NClasses+NClasses-1, V);
        end;
        Inc(I);
    end;
    
    //
    // Calculate training set estimates
    //
    Rep.RelClsError := DFRelClsError(DF, XY, NPoints);
    Rep.AvgCE := DFAvgCE(DF, XY, NPoints);
    Rep.RMSError := DFRMSError(DF, XY, NPoints);
    Rep.AvgError := DFAvgError(DF, XY, NPoints);
    Rep.AvgRelError := DFAvgRelError(DF, XY, NPoints);
    
    //
    // Calculate OOB estimates.
    //
    Rep.OOBRelClsError := 0;
    Rep.OOBAvgCE := 0;
    Rep.OOBRMSError := 0;
    Rep.OOBAvgError := 0;
    Rep.OOBAvgRelError := 0;
    OOBCnt := 0;
    OOBRelCnt := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        if OOBCntBuf[I]<>0 then
        begin
            OOBOffs := I*NClasses;
            if NClasses>1 then
            begin
                
                //
                // classification-specific code
                //
                K := Round(XY[I,NVars]);
                TmpI := 0;
                J:=1;
                while J<=NClasses-1 do
                begin
                    if AP_FP_Greater(OOBBuf[OOBOffs+J],OOBBuf[OOBOffs+TmpI]) then
                    begin
                        TmpI := J;
                    end;
                    Inc(J);
                end;
                if TmpI<>K then
                begin
                    Rep.OOBRelClsError := Rep.OOBRelClsError+1;
                end;
                if AP_FP_Neq(OOBBuf[OOBOffs+K],0) then
                begin
                    Rep.OOBAvgCE := Rep.OOBAvgCE-Ln(OOBBuf[OOBOffs+K]);
                end
                else
                begin
                    Rep.OOBAvgCE := Rep.OOBAvgCE-Ln(MinRealNumber);
                end;
                J:=0;
                while J<=NClasses-1 do
                begin
                    if J=K then
                    begin
                        Rep.OOBRMSError := Rep.OOBRMSError+AP_Sqr(OOBBuf[OOBOffs+J]-1);
                        Rep.OOBAvgError := Rep.OOBAvgError+AbsReal(OOBBuf[OOBOffs+J]-1);
                        Rep.OOBAvgRelError := Rep.OOBAvgRelError+AbsReal(OOBBuf[OOBOffs+J]-1);
                        OOBRelCnt := OOBRelCnt+1;
                    end
                    else
                    begin
                        Rep.OOBRMSError := Rep.OOBRMSError+AP_Sqr(OOBBuf[OOBOffs+J]);
                        Rep.OOBAvgError := Rep.OOBAvgError+AbsReal(OOBBuf[OOBOffs+J]);
                    end;
                    Inc(J);
                end;
            end
            else
            begin
                
                //
                // regression-specific code
                //
                Rep.OOBRMSError := Rep.OOBRMSError+AP_Sqr(OOBBuf[OOBOffs]-XY[I,NVars]);
                Rep.OOBAvgError := Rep.OOBAvgError+AbsReal(OOBBuf[OOBOffs]-XY[I,NVars]);
                if AP_FP_Neq(XY[I,NVars],0) then
                begin
                    Rep.OOBAvgRelError := Rep.OOBAvgRelError+AbsReal((OOBBuf[OOBOffs]-XY[I,NVars])/XY[I,NVars]);
                    OOBRelCnt := OOBRelCnt+1;
                end;
            end;
            
            //
            // update OOB estimates count.
            //
            OOBCnt := OOBCnt+1;
        end;
        Inc(I);
    end;
    if OOBCnt>0 then
    begin
        Rep.OOBRelClsError := Rep.OOBRelClsError/OOBCnt;
        Rep.OOBAvgCE := Rep.OOBAvgCE/OOBCnt;
        Rep.OOBRMSError := Sqrt(Rep.OOBRMSError/(OOBCnt*NClasses));
        Rep.OOBAvgError := Rep.OOBAvgError/(OOBCnt*NClasses);
        if OOBRelCnt>0 then
        begin
            Rep.OOBAvgRelError := Rep.OOBAvgRelError/OOBRelCnt;
        end;
    end;
end;


(*************************************************************************
Procesing

INPUT PARAMETERS:
    DF      -   decision forest model
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NClasses-1].

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure DFProcess(const DF : DecisionForest;
     const X : TReal1DArray;
     var Y : TReal1DArray);
var
    Offs : AlglibInteger;
    I : AlglibInteger;
    V : Double;
begin
    
    //
    // Proceed
    //
    Offs := 0;
    I:=0;
    while I<=DF.NClasses-1 do
    begin
        Y[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=DF.NTrees-1 do
    begin
        
        //
        // Process basic tree
        //
        DFProcessInternal(DF, Offs, X, Y);
        
        //
        // Next tree
        //
        Offs := Offs+Round(DF.Trees[Offs]);
        Inc(I);
    end;
    V := AP_Double(1)/DF.NTrees;
    APVMul(@Y[0], 0, DF.NClasses-1, V);
end;


(*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************)
function DFRelClsError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
begin
    Result := AP_Double(DFClsError(DF, XY, NPoints))/NPoints;
end;


(*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************)
function DFAvgCE(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    X : TReal1DArray;
    Y : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TmpI : AlglibInteger;
begin
    SetLength(X, DF.NVars-1+1);
    SetLength(Y, DF.NClasses-1+1);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, DF.NVars-1, @XY[I][0], 0, DF.NVars-1);
        DFProcess(DF, X, Y);
        if DF.NClasses>1 then
        begin
            
            //
            // classification-specific code
            //
            K := Round(XY[I,DF.NVars]);
            TmpI := 0;
            J:=1;
            while J<=DF.NClasses-1 do
            begin
                if AP_FP_Greater(Y[J],Y[TmpI]) then
                begin
                    TmpI := J;
                end;
                Inc(J);
            end;
            if AP_FP_Neq(Y[K],0) then
            begin
                Result := Result-Ln(Y[K]);
            end
            else
            begin
                Result := Result-Ln(MinRealNumber);
            end;
        end;
        Inc(I);
    end;
    Result := Result/NPoints;
end;


(*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************)
function DFRMSError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    X : TReal1DArray;
    Y : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TmpI : AlglibInteger;
begin
    SetLength(X, DF.NVars-1+1);
    SetLength(Y, DF.NClasses-1+1);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, DF.NVars-1, @XY[I][0], 0, DF.NVars-1);
        DFProcess(DF, X, Y);
        if DF.NClasses>1 then
        begin
            
            //
            // classification-specific code
            //
            K := Round(XY[I,DF.NVars]);
            TmpI := 0;
            J:=1;
            while J<=DF.NClasses-1 do
            begin
                if AP_FP_Greater(Y[J],Y[TmpI]) then
                begin
                    TmpI := J;
                end;
                Inc(J);
            end;
            J:=0;
            while J<=DF.NClasses-1 do
            begin
                if J=K then
                begin
                    Result := Result+AP_Sqr(Y[J]-1);
                end
                else
                begin
                    Result := Result+AP_Sqr(Y[J]);
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // regression-specific code
            //
            Result := Result+AP_Sqr(Y[0]-XY[I,DF.NVars]);
        end;
        Inc(I);
    end;
    Result := Sqrt(Result/(NPoints*DF.NClasses));
end;


(*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************)
function DFAvgError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    X : TReal1DArray;
    Y : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
begin
    SetLength(X, DF.NVars-1+1);
    SetLength(Y, DF.NClasses-1+1);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, DF.NVars-1, @XY[I][0], 0, DF.NVars-1);
        DFProcess(DF, X, Y);
        if DF.NClasses>1 then
        begin
            
            //
            // classification-specific code
            //
            K := Round(XY[I,DF.NVars]);
            J:=0;
            while J<=DF.NClasses-1 do
            begin
                if J=K then
                begin
                    Result := Result+AbsReal(Y[J]-1);
                end
                else
                begin
                    Result := Result+AbsReal(Y[J]);
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // regression-specific code
            //
            Result := Result+AbsReal(Y[0]-XY[I,DF.NVars]);
        end;
        Inc(I);
    end;
    Result := Result/(NPoints*DF.NClasses);
end;


(*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************)
function DFAvgRelError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    X : TReal1DArray;
    Y : TReal1DArray;
    RelCnt : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
begin
    SetLength(X, DF.NVars-1+1);
    SetLength(Y, DF.NClasses-1+1);
    Result := 0;
    RelCnt := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, DF.NVars-1, @XY[I][0], 0, DF.NVars-1);
        DFProcess(DF, X, Y);
        if DF.NClasses>1 then
        begin
            
            //
            // classification-specific code
            //
            K := Round(XY[I,DF.NVars]);
            J:=0;
            while J<=DF.NClasses-1 do
            begin
                if J=K then
                begin
                    Result := Result+AbsReal(Y[J]-1);
                    RelCnt := RelCnt+1;
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // regression-specific code
            //
            if AP_FP_Neq(XY[I,DF.NVars],0) then
            begin
                Result := Result+AbsReal((Y[0]-XY[I,DF.NVars])/XY[I,DF.NVars]);
                RelCnt := RelCnt+1;
            end;
        end;
        Inc(I);
    end;
    if RelCnt>0 then
    begin
        Result := Result/RelCnt;
    end;
end;


(*************************************************************************
Copying of DecisionForest strucure

INPUT PARAMETERS:
    DF1 -   original

OUTPUT PARAMETERS:
    DF2 -   copy

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure DFCopy(const DF1 : DecisionForest; var DF2 : DecisionForest);
begin
    DF2.NVars := DF1.NVars;
    DF2.NClasses := DF1.NClasses;
    DF2.NTrees := DF1.NTrees;
    DF2.BufSize := DF1.BufSize;
    SetLength(DF2.Trees, DF1.BufSize-1+1);
    APVMove(@DF2.Trees[0], 0, DF1.BufSize-1, @DF1.Trees[0], 0, DF1.BufSize-1);
end;


(*************************************************************************
Serialization of DecisionForest strucure

INPUT PARAMETERS:
    DF      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores decision forest,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure DFSerialize(const DF : DecisionForest;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
begin
    SetLength(RA, DF.BufSize+5-1+1);
    RA[0] := DFVNum;
    RA[1] := DF.NVars;
    RA[2] := DF.NClasses;
    RA[3] := DF.NTrees;
    RA[4] := DF.BufSize;
    APVMove(@RA[0], 5, 5+DF.BufSize-1, @DF.Trees[0], 0, DF.BufSize-1);
    RLen := 5+DF.BufSize;
end;


(*************************************************************************
Unserialization of DecisionForest strucure

INPUT PARAMETERS:
    RA      -   real array which stores decision forest

OUTPUT PARAMETERS:
    DF      -   restored structure

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure DFUnserialize(const RA : TReal1DArray; var DF : DecisionForest);
begin
    Assert(Round(RA[0])=DFVNum, 'DFUnserialize: incorrect array!');
    DF.NVars := Round(RA[1]);
    DF.NClasses := Round(RA[2]);
    DF.NTrees := Round(RA[3]);
    DF.BufSize := Round(RA[4]);
    SetLength(DF.Trees, DF.BufSize-1+1);
    APVMove(@DF.Trees[0], 0, DF.BufSize-1, @RA[0], 5, 5+DF.BufSize-1);
end;


(*************************************************************************
Classification error
*************************************************************************)
function DFClsError(const DF : DecisionForest;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):AlglibInteger;
var
    X : TReal1DArray;
    Y : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TmpI : AlglibInteger;
begin
    if DF.NClasses<=1 then
    begin
        Result := 0;
        Exit;
    end;
    SetLength(X, DF.NVars-1+1);
    SetLength(Y, DF.NClasses-1+1);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, DF.NVars-1, @XY[I][0], 0, DF.NVars-1);
        DFProcess(DF, X, Y);
        K := Round(XY[I,DF.NVars]);
        TmpI := 0;
        J:=1;
        while J<=DF.NClasses-1 do
        begin
            if AP_FP_Greater(Y[J],Y[TmpI]) then
            begin
                TmpI := J;
            end;
            Inc(J);
        end;
        if TmpI<>K then
        begin
            Result := Result+1;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal subroutine for processing one decision tree starting at Offs
*************************************************************************)
procedure DFProcessInternal(const DF : DecisionForest;
     Offs : AlglibInteger;
     const X : TReal1DArray;
     var Y : TReal1DArray);
var
    K : AlglibInteger;
    Idx : AlglibInteger;
begin
    
    //
    // Set pointer to the root
    //
    K := Offs+1;
    
    //
    // Navigate through the tree
    //
    while True do
    begin
        if AP_FP_Eq(DF.Trees[K],-1) then
        begin
            if DF.NClasses=1 then
            begin
                Y[0] := Y[0]+DF.Trees[K+1];
            end
            else
            begin
                Idx := Round(DF.Trees[K+1]);
                Y[Idx] := Y[Idx]+1;
            end;
            Break;
        end;
        if AP_FP_Less(X[Round(DF.Trees[K])],DF.Trees[K+1]) then
        begin
            K := K+InnerNodeWidth;
        end
        else
        begin
            K := Offs+Round(DF.Trees[K+2]);
        end;
    end;
end;


(*************************************************************************
Builds one decision tree. Just a wrapper for the DFBuildTreeRec.
*************************************************************************)
procedure DFBuildTree(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NFeatures : AlglibInteger;
     NVarsInPool : AlglibInteger;
     Flags : AlglibInteger;
     var Bufs : DFInternalBuffers);
var
    NumProcessed : AlglibInteger;
    I : AlglibInteger;
begin
    Assert(NPoints>0);
    
    //
    // Prepare IdxBuf. It stores indices of the training set elements.
    // When training set is being split, contents of IdxBuf is
    // correspondingly reordered so we can know which elements belong
    // to which branch of decision tree.
    //
    I:=0;
    while I<=NPoints-1 do
    begin
        Bufs.IdxBuf[I] := I;
        Inc(I);
    end;
    
    //
    // Recursive procedure
    //
    NumProcessed := 1;
    DFBuildTreeRec(XY, NPoints, NVars, NClasses, NFeatures, NVarsInPool, Flags, NumProcessed, 0, NPoints-1, Bufs);
    Bufs.TreeBuf[0] := NumProcessed;
end;


(*************************************************************************
Builds one decision tree (internal recursive subroutine)

Parameters:
    TreeBuf     -   large enough array, at least TreeSize
    IdxBuf      -   at least NPoints elements
    TmpBufR     -   at least NPoints
    TmpBufR2    -   at least NPoints
    TmpBufI     -   at least NPoints
    TmpBufI2    -   at least NPoints+1
*************************************************************************)
procedure DFBuildTreeRec(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     NFeatures : AlglibInteger;
     NVarsInPool : AlglibInteger;
     Flags : AlglibInteger;
     var NumProcessed : AlglibInteger;
     Idx1 : AlglibInteger;
     Idx2 : AlglibInteger;
     var Bufs : DFInternalBuffers);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    BFlag : Boolean;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    Info : AlglibInteger;
    SL : Double;
    SR : Double;
    W : Double;
    IdxBest : AlglibInteger;
    EBest : Double;
    TBest : Double;
    VarCur : AlglibInteger;
    S : Double;
    V : Double;
    V1 : Double;
    V2 : Double;
    Threshold : Double;
    OldNP : AlglibInteger;
    CurRMS : Double;
    UseEVS : Boolean;
begin
    Assert(NPoints>0);
    Assert(Idx2>=Idx1);
    UseEVS := Flags div DFUseEVS mod 2<>0;
    
    //
    // Leaf node
    //
    if Idx2=Idx1 then
    begin
        Bufs.TreeBuf[NumProcessed] := -1;
        Bufs.TreeBuf[NumProcessed+1] := XY[Bufs.IdxBuf[Idx1],NVars];
        NumProcessed := NumProcessed+LeafNodeWidth;
        Exit;
    end;
    
    //
    // Non-leaf node.
    // Select random variable, prepare split:
    // 1. prepare default solution - no splitting, class at random
    // 2. investigate possible splits, compare with default/best
    //
    IdxBest := -1;
    if NClasses>1 then
    begin
        
        //
        // default solution for classification
        //
        I:=0;
        while I<=NClasses-1 do
        begin
            Bufs.ClassIBuf[I] := 0;
            Inc(I);
        end;
        S := Idx2-Idx1+1;
        I:=Idx1;
        while I<=Idx2 do
        begin
            J := Round(XY[Bufs.IdxBuf[I],NVars]);
            Bufs.ClassIBuf[J] := Bufs.ClassIBuf[J]+1;
            Inc(I);
        end;
        EBest := 0;
        I:=0;
        while I<=NClasses-1 do
        begin
            EBest := EBest+Bufs.ClassIBuf[I]*AP_Sqr(1-Bufs.ClassIBuf[I]/S)+(S-Bufs.ClassIBuf[I])*AP_Sqr(Bufs.ClassIBuf[I]/S);
            Inc(I);
        end;
        EBest := Sqrt(EBest/(NClasses*(Idx2-Idx1+1)));
    end
    else
    begin
        
        //
        // default solution for regression
        //
        V := 0;
        I:=Idx1;
        while I<=Idx2 do
        begin
            V := V+XY[Bufs.IdxBuf[I],NVars];
            Inc(I);
        end;
        V := V/(Idx2-Idx1+1);
        EBest := 0;
        I:=Idx1;
        while I<=Idx2 do
        begin
            EBest := EBest+AP_Sqr(XY[Bufs.IdxBuf[I],NVars]-V);
            Inc(I);
        end;
        EBest := Sqrt(EBest/(Idx2-Idx1+1));
    end;
    I := 0;
    while I<=Min(NFeatures, NVarsInPool)-1 do
    begin
        
        //
        // select variables from pool
        //
        J := I+RandomInteger(NVarsInPool-I);
        K := Bufs.VarPool[I];
        Bufs.VarPool[I] := Bufs.VarPool[J];
        Bufs.VarPool[J] := K;
        VarCur := Bufs.VarPool[I];
        
        //
        // load variable values to working array
        //
        // apply EVS preprocessing: if all variable values are same,
        // variable is excluded from pool.
        //
        // This is necessary for binary pre-splits (see later) to work.
        //
        J:=Idx1;
        while J<=Idx2 do
        begin
            Bufs.TmpBufR[J-Idx1] := XY[Bufs.IdxBuf[J],VarCur];
            Inc(J);
        end;
        if UseEVS then
        begin
            BFlag := False;
            V := Bufs.TmpBufR[0];
            J:=0;
            while J<=Idx2-Idx1 do
            begin
                if AP_FP_Neq(Bufs.TmpBufR[J],V) then
                begin
                    BFlag := True;
                    Break;
                end;
                Inc(J);
            end;
            if  not BFlag then
            begin
                
                //
                // exclude variable from pool,
                // go to the next iteration.
                // I is not increased.
                //
                K := Bufs.VarPool[I];
                Bufs.VarPool[I] := Bufs.VarPool[NVarsInPool-1];
                Bufs.VarPool[NVarsInPool-1] := K;
                NVarsInPool := NVarsInPool-1;
                Continue;
            end;
        end;
        
        //
        // load labels to working array
        //
        if NClasses>1 then
        begin
            J:=Idx1;
            while J<=Idx2 do
            begin
                Bufs.TmpBufI[J-Idx1] := Round(XY[Bufs.IdxBuf[J],NVars]);
                Inc(J);
            end;
        end
        else
        begin
            J:=Idx1;
            while J<=Idx2 do
            begin
                Bufs.TmpBufR2[J-Idx1] := XY[Bufs.IdxBuf[J],NVars];
                Inc(J);
            end;
        end;
        
        //
        // calculate split
        //
        if UseEVS and Bufs.EVSBin[VarCur] then
        begin
            
            //
            // Pre-calculated splits for binary variables.
            // Threshold is already known, just calculate RMS error
            //
            Threshold := Bufs.EVSSplits[VarCur];
            if NClasses>1 then
            begin
                
                //
                // classification-specific code
                //
                J:=0;
                while J<=2*NClasses-1 do
                begin
                    Bufs.ClassIBuf[J] := 0;
                    Inc(J);
                end;
                SL := 0;
                SR := 0;
                J:=0;
                while J<=Idx2-Idx1 do
                begin
                    K := Bufs.TmpBufI[J];
                    if AP_FP_Less(Bufs.TmpBufR[J],Threshold) then
                    begin
                        Bufs.ClassIBuf[K] := Bufs.ClassIBuf[K]+1;
                        SL := SL+1;
                    end
                    else
                    begin
                        Bufs.ClassIBuf[K+NClasses] := Bufs.ClassIBuf[K+NClasses]+1;
                        SR := SR+1;
                    end;
                    Inc(J);
                end;
                Assert(AP_FP_Neq(SL,0) and AP_FP_Neq(SR,0), 'DFBuildTreeRec: something strange!');
                CurRMS := 0;
                J:=0;
                while J<=NClasses-1 do
                begin
                    W := Bufs.ClassIBuf[J];
                    CurRMS := CurRMS+W*AP_Sqr(W/SL-1);
                    CurRMS := CurRMS+(SL-W)*AP_Sqr(W/SL);
                    W := Bufs.ClassIBuf[NClasses+J];
                    CurRMS := CurRMS+W*AP_Sqr(W/SR-1);
                    CurRMS := CurRMS+(SR-W)*AP_Sqr(W/SR);
                    Inc(J);
                end;
                CurRMS := Sqrt(CurRMS/(NClasses*(Idx2-Idx1+1)));
            end
            else
            begin
                
                //
                // regression-specific code
                //
                SL := 0;
                SR := 0;
                V1 := 0;
                V2 := 0;
                J:=0;
                while J<=Idx2-Idx1 do
                begin
                    if AP_FP_Less(Bufs.TmpBufR[J],Threshold) then
                    begin
                        V1 := V1+Bufs.TmpBufR2[J];
                        SL := SL+1;
                    end
                    else
                    begin
                        V2 := V2+Bufs.TmpBufR2[J];
                        SR := SR+1;
                    end;
                    Inc(J);
                end;
                Assert(AP_FP_Neq(SL,0) and AP_FP_Neq(SR,0), 'DFBuildTreeRec: something strange!');
                V1 := V1/SL;
                V2 := V2/SR;
                CurRMS := 0;
                J:=0;
                while J<=Idx2-Idx1 do
                begin
                    if AP_FP_Less(Bufs.TmpBufR[J],Threshold) then
                    begin
                        CurRMS := CurRMS+AP_Sqr(V1-Bufs.TmpBufR2[J]);
                    end
                    else
                    begin
                        CurRMS := CurRMS+AP_Sqr(V2-Bufs.TmpBufR2[J]);
                    end;
                    Inc(J);
                end;
                CurRMS := Sqrt(CurRMS/(Idx2-Idx1+1));
            end;
            Info := 1;
        end
        else
        begin
            
            //
            // Generic splits
            //
            if NClasses>1 then
            begin
                DFSplitC(Bufs.TmpBufR, Bufs.TmpBufI, Bufs.ClassIBuf, Idx2-Idx1+1, NClasses, DFUseStrongSplits, Info, Threshold, CurRMS);
            end
            else
            begin
                DFSplitR(Bufs.TmpBufR, Bufs.TmpBufR2, Idx2-Idx1+1, DFUseStrongSplits, Info, Threshold, CurRMS);
            end;
        end;
        if Info>0 then
        begin
            if AP_FP_Less_Eq(CurRMS,EBest) then
            begin
                EBest := CurRMS;
                IdxBest := VarCur;
                TBest := Threshold;
            end;
        end;
        
        //
        // Next iteration
        //
        I := I+1;
    end;
    
    //
    // to split or not to split
    //
    if IdxBest<0 then
    begin
        
        //
        // All values are same, cannot split.
        //
        Bufs.TreeBuf[NumProcessed] := -1;
        if NClasses>1 then
        begin
            
            //
            // Select random class label (randomness allows us to
            // approximate distribution of the classes)
            //
            Bufs.TreeBuf[NumProcessed+1] := Round(XY[Bufs.IdxBuf[Idx1+RandomInteger(Idx2-Idx1+1)],NVars]);
        end
        else
        begin
            
            //
            // Select average (for regression task).
            //
            V := 0;
            I:=Idx1;
            while I<=Idx2 do
            begin
                V := V+XY[Bufs.IdxBuf[I],NVars]/(Idx2-Idx1+1);
                Inc(I);
            end;
            Bufs.TreeBuf[NumProcessed+1] := V;
        end;
        NumProcessed := NumProcessed+LeafNodeWidth;
    end
    else
    begin
        
        //
        // we can split
        //
        Bufs.TreeBuf[NumProcessed] := IdxBest;
        Bufs.TreeBuf[NumProcessed+1] := TBest;
        I1 := Idx1;
        I2 := Idx2;
        while I1<=I2 do
        begin
            
            //
            // Reorder indices so that left partition is in [Idx1..I1-1],
            // and right partition is in [I2+1..Idx2]
            //
            if AP_FP_Less(XY[Bufs.IdxBuf[I1],IdxBest],TBest) then
            begin
                I1 := I1+1;
                Continue;
            end;
            if AP_FP_Greater_Eq(XY[Bufs.IdxBuf[I2],IdxBest],TBest) then
            begin
                I2 := I2-1;
                Continue;
            end;
            J := Bufs.IdxBuf[I1];
            Bufs.IdxBuf[I1] := Bufs.IdxBuf[I2];
            Bufs.IdxBuf[I2] := J;
            I1 := I1+1;
            I2 := I2-1;
        end;
        OldNP := NumProcessed;
        NumProcessed := NumProcessed+InnerNodeWidth;
        DFBuildTreeRec(XY, NPoints, NVars, NClasses, NFeatures, NVarsInPool, Flags, NumProcessed, Idx1, I1-1, Bufs);
        Bufs.TreeBuf[OldNP+2] := NumProcessed;
        DFBuildTreeRec(XY, NPoints, NVars, NClasses, NFeatures, NVarsInPool, Flags, NumProcessed, I2+1, Idx2, Bufs);
    end;
end;


(*************************************************************************
Makes weak split on attribute
*************************************************************************)
procedure DFWeakSplitI(var X : TReal1DArray;
     var Y : TInteger1DArray;
     N : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var E : Double);
var
    I : AlglibInteger;
    NEq : AlglibInteger;
    NLess : AlglibInteger;
    NGreater : AlglibInteger;
begin
    TagSortFastI(X, Y, N);
    if N mod 2=1 then
    begin
        
        //
        // odd number of elements
        //
        Threshold := X[N div 2];
    end
    else
    begin
        
        //
        // even number of elements.
        //
        // if two closest to the middle of the array are equal,
        // we will select one of them (to avoid possible problems with
        // floating point errors).
        // we will select halfsum otherwise.
        //
        if AP_FP_Eq(X[N div 2-1],X[N div 2]) then
        begin
            Threshold := X[N div 2-1];
        end
        else
        begin
            Threshold := Double(0.5)*(X[N div 2-1]+X[N div 2]);
        end;
    end;
    NEq := 0;
    NLess := 0;
    NGreater := 0;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Less(X[I],Threshold) then
        begin
            NLess := NLess+1;
        end;
        if AP_FP_Eq(X[I],Threshold) then
        begin
            NEq := NEq+1;
        end;
        if AP_FP_Greater(X[I],Threshold) then
        begin
            NGreater := NGreater+1;
        end;
        Inc(I);
    end;
    if (NLess=0) and (NGreater=0) then
    begin
        Info := -3;
    end
    else
    begin
        if NEq<>0 then
        begin
            if NLess<NGreater then
            begin
                Threshold := Double(0.5)*(X[NLess+NEq-1]+X[NLess+NEq]);
            end
            else
            begin
                Threshold := Double(0.5)*(X[NLess-1]+X[NLess]);
            end;
        end;
        Info := 1;
        E := 0;
    end;
end;


(*************************************************************************
Makes split on attribute
*************************************************************************)
procedure DFSplitC(var X : TReal1DArray;
     var C : TInteger1DArray;
     var CntBuf : TInteger1DArray;
     N : AlglibInteger;
     NC : AlglibInteger;
     Flags : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var E : Double);
var
    I : AlglibInteger;
    NEq : AlglibInteger;
    NLess : AlglibInteger;
    NGreater : AlglibInteger;
    Q : AlglibInteger;
    QMin : AlglibInteger;
    QMax : AlglibInteger;
    QCnt : AlglibInteger;
    CurSplit : Double;
    NLeft : AlglibInteger;
    V : Double;
    CurE : Double;
    W : Double;
    SL : Double;
    SR : Double;
begin
    TagSortFastI(X, C, N);
    E := MaxRealNumber;
    Threshold := Double(0.5)*(X[0]+X[N-1]);
    Info := -3;
    if Flags div DFUseStrongSplits mod 2=0 then
    begin
        
        //
        // weak splits, split at half
        //
        QCnt := 2;
        QMin := 1;
        QMax := 1;
    end
    else
    begin
        
        //
        // strong splits: choose best quartile
        //
        QCnt := 4;
        QMin := 1;
        QMax := 3;
    end;
    Q:=QMin;
    while Q<=QMax do
    begin
        CurSplit := X[N*Q div QCnt];
        NEq := 0;
        NLess := 0;
        NGreater := 0;
        I:=0;
        while I<=N-1 do
        begin
            if AP_FP_Less(X[I],CurSplit) then
            begin
                NLess := NLess+1;
            end;
            if AP_FP_Eq(X[I],CurSplit) then
            begin
                NEq := NEq+1;
            end;
            if AP_FP_Greater(X[I],CurSplit) then
            begin
                NGreater := NGreater+1;
            end;
            Inc(I);
        end;
        Assert(NEq<>0, 'DFSplitR: NEq=0, something strange!!!');
        if (NLess<>0) or (NGreater<>0) then
        begin
            
            //
            // set threshold between two partitions, with
            // some tweaking to avoid problems with floating point
            // arithmetics.
            //
            // The problem is that when you calculates C = 0.5*(A+B) there
            // can be no C which lies strictly between A and B (for example,
            // there is no floating point number which is
            // greater than 1 and less than 1+eps). In such situations
            // we choose right side as theshold (remember that
            // points which lie on threshold falls to the right side).
            //
            if NLess<NGreater then
            begin
                CurSplit := Double(0.5)*(X[NLess+NEq-1]+X[NLess+NEq]);
                NLeft := NLess+NEq;
                if AP_FP_Less_Eq(CurSplit,X[NLess+NEq-1]) then
                begin
                    CurSplit := X[NLess+NEq];
                end;
            end
            else
            begin
                CurSplit := Double(0.5)*(X[NLess-1]+X[NLess]);
                NLeft := NLess;
                if AP_FP_Less_Eq(CurSplit,X[NLess-1]) then
                begin
                    CurSplit := X[NLess];
                end;
            end;
            Info := 1;
            CurE := 0;
            I:=0;
            while I<=2*NC-1 do
            begin
                CntBuf[I] := 0;
                Inc(I);
            end;
            I:=0;
            while I<=NLeft-1 do
            begin
                CntBuf[C[I]] := CntBuf[C[I]]+1;
                Inc(I);
            end;
            I:=NLeft;
            while I<=N-1 do
            begin
                CntBuf[NC+C[I]] := CntBuf[NC+C[I]]+1;
                Inc(I);
            end;
            SL := NLeft;
            SR := N-NLeft;
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
            CurE := Sqrt(V/(NC*N));
            if AP_FP_Less(CurE,E) then
            begin
                Threshold := CurSplit;
                E := CurE;
            end;
        end;
        Inc(Q);
    end;
end;


(*************************************************************************
Makes split on attribute
*************************************************************************)
procedure DFSplitR(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger;
     Flags : AlglibInteger;
     var Info : AlglibInteger;
     var Threshold : Double;
     var E : Double);
var
    I : AlglibInteger;
    NEq : AlglibInteger;
    NLess : AlglibInteger;
    NGreater : AlglibInteger;
    Q : AlglibInteger;
    QMin : AlglibInteger;
    QMax : AlglibInteger;
    QCnt : AlglibInteger;
    CurSplit : Double;
    NLeft : AlglibInteger;
    V : Double;
    CurE : Double;
begin
    TagSortFastR(X, Y, N);
    E := MaxRealNumber;
    Threshold := Double(0.5)*(X[0]+X[N-1]);
    Info := -3;
    if Flags div DFUseStrongSplits mod 2=0 then
    begin
        
        //
        // weak splits, split at half
        //
        QCnt := 2;
        QMin := 1;
        QMax := 1;
    end
    else
    begin
        
        //
        // strong splits: choose best quartile
        //
        QCnt := 4;
        QMin := 1;
        QMax := 3;
    end;
    Q:=QMin;
    while Q<=QMax do
    begin
        CurSplit := X[N*Q div QCnt];
        NEq := 0;
        NLess := 0;
        NGreater := 0;
        I:=0;
        while I<=N-1 do
        begin
            if AP_FP_Less(X[I],CurSplit) then
            begin
                NLess := NLess+1;
            end;
            if AP_FP_Eq(X[I],CurSplit) then
            begin
                NEq := NEq+1;
            end;
            if AP_FP_Greater(X[I],CurSplit) then
            begin
                NGreater := NGreater+1;
            end;
            Inc(I);
        end;
        Assert(NEq<>0, 'DFSplitR: NEq=0, something strange!!!');
        if (NLess<>0) or (NGreater<>0) then
        begin
            
            //
            // set threshold between two partitions, with
            // some tweaking to avoid problems with floating point
            // arithmetics.
            //
            // The problem is that when you calculates C = 0.5*(A+B) there
            // can be no C which lies strictly between A and B (for example,
            // there is no floating point number which is
            // greater than 1 and less than 1+eps). In such situations
            // we choose right side as theshold (remember that
            // points which lie on threshold falls to the right side).
            //
            if NLess<NGreater then
            begin
                CurSplit := Double(0.5)*(X[NLess+NEq-1]+X[NLess+NEq]);
                NLeft := NLess+NEq;
                if AP_FP_Less_Eq(CurSplit,X[NLess+NEq-1]) then
                begin
                    CurSplit := X[NLess+NEq];
                end;
            end
            else
            begin
                CurSplit := Double(0.5)*(X[NLess-1]+X[NLess]);
                NLeft := NLess;
                if AP_FP_Less_Eq(CurSplit,X[NLess-1]) then
                begin
                    CurSplit := X[NLess];
                end;
            end;
            Info := 1;
            CurE := 0;
            V := 0;
            I:=0;
            while I<=NLeft-1 do
            begin
                V := V+Y[I];
                Inc(I);
            end;
            V := V/NLeft;
            I:=0;
            while I<=NLeft-1 do
            begin
                CurE := CurE+AP_Sqr(Y[I]-V);
                Inc(I);
            end;
            V := 0;
            I:=NLeft;
            while I<=N-1 do
            begin
                V := V+Y[I];
                Inc(I);
            end;
            V := V/(N-NLeft);
            I:=NLeft;
            while I<=N-1 do
            begin
                CurE := CurE+AP_Sqr(Y[I]-V);
                Inc(I);
            end;
            CurE := Sqrt(CurE/N);
            if AP_FP_Less(CurE,E) then
            begin
                Threshold := CurSplit;
                E := CurE;
            end;
        end;
        Inc(Q);
    end;
end;


end.