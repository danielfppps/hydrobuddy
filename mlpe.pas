{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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
unit mlpe;
interface
uses Math, Sysutils, Ap, mlpbase, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, linmin, minlbfgs, hblas, sblas, ortfac, blas, rotations, bdsvd, svd, xblas, densesolver, mlptrain, tsort, descriptivestatistics, bdss;

type
(*************************************************************************
Neural networks ensemble
*************************************************************************)
MLPEnsemble = record
    StructInfo : TInteger1DArray;
    EnsembleSize : AlglibInteger;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    WCount : AlglibInteger;
    IsSoftmax : Boolean;
    PostProcessing : Boolean;
    Weights : TReal1DArray;
    ColumnMeans : TReal1DArray;
    ColumnSigmas : TReal1DArray;
    SerializedLen : AlglibInteger;
    SerializedMLP : TReal1DArray;
    TmpWeights : TReal1DArray;
    TmpMeans : TReal1DArray;
    TmpSigmas : TReal1DArray;
    Neurons : TReal1DArray;
    DFDNET : TReal1DArray;
    Y : TReal1DArray;
end;



procedure MLPECreate0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreate1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreate2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateB0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     B : Double;
     D : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateB1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     B : Double;
     D : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateB2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     B : Double;
     D : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateR0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     A : Double;
     B : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateR1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     A : Double;
     B : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateR2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     A : Double;
     B : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateC0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateC1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateC2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECreateFromNetwork(const Network : MultiLayerPerceptron;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
procedure MLPECopy(const Ensemble1 : MLPEnsemble; var Ensemble2 : MLPEnsemble);
procedure MLPESerialize(var Ensemble : MLPEnsemble;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
procedure MLPEUnserialize(const RA : TReal1DArray; var Ensemble : MLPEnsemble);
procedure MLPERandomize(var Ensemble : MLPEnsemble);
procedure MLPEProperties(const Ensemble : MLPEnsemble;
     var NIn : AlglibInteger;
     var NOut : AlglibInteger);
function MLPEIsSoftmax(const Ensemble : MLPEnsemble):Boolean;
procedure MLPEProcess(var Ensemble : MLPEnsemble;
     const X : TReal1DArray;
     var Y : TReal1DArray);
function MLPERelClsError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MLPEAvgCE(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MLPERMSError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MLPEAvgError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MLPEAvgRelError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
procedure MLPEBaggingLM(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var OOBErrors : MLPCVReport);
procedure MLPEBaggingLBFGS(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var OOBErrors : MLPCVReport);
procedure MLPETrainES(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);

implementation

const
    MLPNTotalOffset = 3;
    MLPEVNum = 9;

procedure MLPEAllErrors(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     var RelCls : Double;
     var AvgCE : Double;
     var RMS : Double;
     var Avg : Double;
     var AvgRel : Double);forward;
procedure MLPEBaggingInternal(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     LMAlgorithm : Boolean;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var OOBErrors : MLPCVReport);forward;


(*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreate0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreate0(NIn, NOut, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreate1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreate1(NIn, NHid, NOut, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreate2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreate2(NIn, NHid1, NHid2, NOut, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateB0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     B : Double;
     D : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateB0(NIn, NOut, B, D, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateB1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     B : Double;
     D : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateB1(NIn, NHid, NOut, B, D, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateB2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     B : Double;
     D : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateB2(NIn, NHid1, NHid2, NOut, B, D, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateR0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     A : Double;
     B : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateR0(NIn, NOut, A, B, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateR1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     A : Double;
     B : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateR1(NIn, NHid, NOut, A, B, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateR2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     A : Double;
     B : Double;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateR2(NIn, NHid1, NHid2, NOut, A, B, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateC0(NIn : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateC0(NIn, NOut, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateC1(NIn : AlglibInteger;
     NHid : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateC1(NIn, NHid, NOut, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateC2(NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    Net : MultiLayerPerceptron;
begin
    MLPCreateC2(NIn, NHid1, NHid2, NOut, Net);
    MLPECreateFromNetwork(Net, EnsembleSize, Ensemble);
end;


(*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECreateFromNetwork(const Network : MultiLayerPerceptron;
     EnsembleSize : AlglibInteger;
     var Ensemble : MLPEnsemble);
var
    I : AlglibInteger;
    CCount : AlglibInteger;
begin
    Assert(EnsembleSize>0, 'MLPECreate: incorrect ensemble size!');
    
    //
    // network properties
    //
    MLPProperties(Network, Ensemble.NIn, Ensemble.NOut, Ensemble.WCount);
    if MLPIsSoftmax(Network) then
    begin
        CCount := Ensemble.NIn;
    end
    else
    begin
        CCount := Ensemble.NIn+Ensemble.NOut;
    end;
    Ensemble.PostProcessing := False;
    Ensemble.IsSoftmax := MLPIsSoftmax(Network);
    Ensemble.EnsembleSize := EnsembleSize;
    
    //
    // structure information
    //
    SetLength(Ensemble.StructInfo, Network.StructInfo[0]-1+1);
    I:=0;
    while I<=Network.StructInfo[0]-1 do
    begin
        Ensemble.StructInfo[I] := Network.StructInfo[I];
        Inc(I);
    end;
    
    //
    // weights, means, sigmas
    //
    SetLength(Ensemble.Weights, EnsembleSize*Ensemble.WCount-1+1);
    SetLength(Ensemble.ColumnMeans, EnsembleSize*CCount-1+1);
    SetLength(Ensemble.ColumnSigmas, EnsembleSize*CCount-1+1);
    I:=0;
    while I<=EnsembleSize*Ensemble.WCount-1 do
    begin
        Ensemble.Weights[I] := RandomReal-Double(0.5);
        Inc(I);
    end;
    I:=0;
    while I<=EnsembleSize-1 do
    begin
        APVMove(@Ensemble.ColumnMeans[0], I*CCount, (I+1)*CCount-1, @Network.ColumnMeans[0], 0, CCount-1);
        APVMove(@Ensemble.ColumnSigmas[0], I*CCount, (I+1)*CCount-1, @Network.ColumnSigmas[0], 0, CCount-1);
        Inc(I);
    end;
    
    //
    // serialized part
    //
    MLPSerialize(Network, Ensemble.SerializedMLP, Ensemble.SerializedLen);
    
    //
    // temporaries, internal buffers
    //
    SetLength(Ensemble.TmpWeights, Ensemble.WCount-1+1);
    SetLength(Ensemble.TmpMeans, CCount-1+1);
    SetLength(Ensemble.TmpSigmas, CCount-1+1);
    SetLength(Ensemble.Neurons, Ensemble.StructInfo[MLPNTotalOffset]-1+1);
    SetLength(Ensemble.DFDNET, Ensemble.StructInfo[MLPNTotalOffset]-1+1);
    SetLength(Ensemble.Y, Ensemble.NOut-1+1);
end;


(*************************************************************************
Copying of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble1 -   original

OUTPUT PARAMETERS:
    Ensemble2 -   copy

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPECopy(const Ensemble1 : MLPEnsemble; var Ensemble2 : MLPEnsemble);
var
    I : AlglibInteger;
    SSize : AlglibInteger;
    CCount : AlglibInteger;
    NTotal : AlglibInteger;
begin
    
    //
    // Unload info
    //
    SSize := Ensemble1.StructInfo[0];
    if Ensemble1.IsSoftmax then
    begin
        CCount := Ensemble1.NIn;
    end
    else
    begin
        CCount := Ensemble1.NIn+Ensemble1.NOut;
    end;
    NTotal := Ensemble1.StructInfo[MLPNTotalOffset];
    
    //
    // Allocate space
    //
    SetLength(Ensemble2.StructInfo, SSize-1+1);
    SetLength(Ensemble2.Weights, Ensemble1.EnsembleSize*Ensemble1.WCount-1+1);
    SetLength(Ensemble2.ColumnMeans, Ensemble1.EnsembleSize*CCount-1+1);
    SetLength(Ensemble2.ColumnSigmas, Ensemble1.EnsembleSize*CCount-1+1);
    SetLength(Ensemble2.TmpWeights, Ensemble1.WCount-1+1);
    SetLength(Ensemble2.TmpMeans, CCount-1+1);
    SetLength(Ensemble2.TmpSigmas, CCount-1+1);
    SetLength(Ensemble2.SerializedMLP, Ensemble1.SerializedLen-1+1);
    SetLength(Ensemble2.Neurons, NTotal-1+1);
    SetLength(Ensemble2.DFDNET, NTotal-1+1);
    SetLength(Ensemble2.Y, Ensemble1.NOut-1+1);
    
    //
    // Copy
    //
    Ensemble2.NIn := Ensemble1.NIn;
    Ensemble2.NOut := Ensemble1.NOut;
    Ensemble2.WCount := Ensemble1.WCount;
    Ensemble2.EnsembleSize := Ensemble1.EnsembleSize;
    Ensemble2.IsSoftmax := Ensemble1.IsSoftmax;
    Ensemble2.PostProcessing := Ensemble1.PostProcessing;
    Ensemble2.SerializedLen := Ensemble1.SerializedLen;
    I:=0;
    while I<=SSize-1 do
    begin
        Ensemble2.StructInfo[I] := Ensemble1.StructInfo[I];
        Inc(I);
    end;
    APVMove(@Ensemble2.Weights[0], 0, Ensemble1.EnsembleSize*Ensemble1.WCount-1, @Ensemble1.Weights[0], 0, Ensemble1.EnsembleSize*Ensemble1.WCount-1);
    APVMove(@Ensemble2.ColumnMeans[0], 0, Ensemble1.EnsembleSize*CCount-1, @Ensemble1.ColumnMeans[0], 0, Ensemble1.EnsembleSize*CCount-1);
    APVMove(@Ensemble2.ColumnSigmas[0], 0, Ensemble1.EnsembleSize*CCount-1, @Ensemble1.ColumnSigmas[0], 0, Ensemble1.EnsembleSize*CCount-1);
    APVMove(@Ensemble2.SerializedMLP[0], 0, Ensemble1.SerializedLen-1, @Ensemble1.SerializedMLP[0], 0, Ensemble1.SerializedLen-1);
end;


(*************************************************************************
Serialization of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble-   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores ensemble,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPESerialize(var Ensemble : MLPEnsemble;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
var
    I : AlglibInteger;
    SSize : AlglibInteger;
    NTotal : AlglibInteger;
    CCount : AlglibInteger;
    HSize : AlglibInteger;
    Offs : AlglibInteger;
begin
    HSize := 13;
    SSize := Ensemble.StructInfo[0];
    if Ensemble.IsSoftmax then
    begin
        CCount := Ensemble.NIn;
    end
    else
    begin
        CCount := Ensemble.NIn+Ensemble.NOut;
    end;
    NTotal := Ensemble.StructInfo[MLPNTotalOffset];
    RLen := HSize+SSize+Ensemble.EnsembleSize*Ensemble.WCount+2*CCount*Ensemble.EnsembleSize+Ensemble.SerializedLen;
    
    //
    //  RA format:
    //  [0]     RLen
    //  [1]     Version (MLPEVNum)
    //  [2]     EnsembleSize
    //  [3]     NIn
    //  [4]     NOut
    //  [5]     WCount
    //  [6]     IsSoftmax 0/1
    //  [7]     PostProcessing 0/1
    //  [8]     sizeof(StructInfo)
    //  [9]     NTotal (sizeof(Neurons), sizeof(DFDNET))
    //  [10]    CCount (sizeof(ColumnMeans), sizeof(ColumnSigmas))
    //  [11]    data offset
    //  [12]    SerializedLen
    //
    //  [..]    StructInfo
    //  [..]    Weights
    //  [..]    ColumnMeans
    //  [..]    ColumnSigmas
    //
    SetLength(RA, RLen-1+1);
    RA[0] := RLen;
    RA[1] := MLPEVNum;
    RA[2] := Ensemble.EnsembleSize;
    RA[3] := Ensemble.NIn;
    RA[4] := Ensemble.NOut;
    RA[5] := Ensemble.WCount;
    if Ensemble.IsSoftmax then
    begin
        RA[6] := 1;
    end
    else
    begin
        RA[6] := 0;
    end;
    if Ensemble.PostProcessing then
    begin
        RA[7] := 1;
    end
    else
    begin
        RA[7] := 9;
    end;
    RA[8] := SSize;
    RA[9] := NTotal;
    RA[10] := CCount;
    RA[11] := HSize;
    RA[12] := Ensemble.SerializedLen;
    Offs := HSize;
    I:=Offs;
    while I<=Offs+SSize-1 do
    begin
        RA[I] := Ensemble.StructInfo[I-Offs];
        Inc(I);
    end;
    Offs := Offs+SSize;
    APVMove(@RA[0], Offs, Offs+Ensemble.EnsembleSize*Ensemble.WCount-1, @Ensemble.Weights[0], 0, Ensemble.EnsembleSize*Ensemble.WCount-1);
    Offs := Offs+Ensemble.EnsembleSize*Ensemble.WCount;
    APVMove(@RA[0], Offs, Offs+Ensemble.EnsembleSize*CCount-1, @Ensemble.ColumnMeans[0], 0, Ensemble.EnsembleSize*CCount-1);
    Offs := Offs+Ensemble.EnsembleSize*CCount;
    APVMove(@RA[0], Offs, Offs+Ensemble.EnsembleSize*CCount-1, @Ensemble.ColumnSigmas[0], 0, Ensemble.EnsembleSize*CCount-1);
    Offs := Offs+Ensemble.EnsembleSize*CCount;
    APVMove(@RA[0], Offs, Offs+Ensemble.SerializedLen-1, @Ensemble.SerializedMLP[0], 0, Ensemble.SerializedLen-1);
    Offs := Offs+Ensemble.SerializedLen;
end;


(*************************************************************************
Unserialization of MLPEnsemble strucure

INPUT PARAMETERS:
    RA      -   real array which stores ensemble

OUTPUT PARAMETERS:
    Ensemble-   restored structure

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEUnserialize(const RA : TReal1DArray; var Ensemble : MLPEnsemble);
var
    I : AlglibInteger;
    SSize : AlglibInteger;
    NTotal : AlglibInteger;
    CCount : AlglibInteger;
    HSize : AlglibInteger;
    Offs : AlglibInteger;
begin
    Assert(Round(RA[1])=MLPEVNum, 'MLPEUnserialize: incorrect array!');
    
    //
    // load info
    //
    HSize := 13;
    Ensemble.EnsembleSize := Round(RA[2]);
    Ensemble.NIn := Round(RA[3]);
    Ensemble.NOut := Round(RA[4]);
    Ensemble.WCount := Round(RA[5]);
    Ensemble.IsSoftmax := Round(RA[6])=1;
    Ensemble.PostProcessing := Round(RA[7])=1;
    SSize := Round(RA[8]);
    NTotal := Round(RA[9]);
    CCount := Round(RA[10]);
    Offs := Round(RA[11]);
    Ensemble.SerializedLen := Round(RA[12]);
    
    //
    //  Allocate arrays
    //
    SetLength(Ensemble.StructInfo, SSize-1+1);
    SetLength(Ensemble.Weights, Ensemble.EnsembleSize*Ensemble.WCount-1+1);
    SetLength(Ensemble.ColumnMeans, Ensemble.EnsembleSize*CCount-1+1);
    SetLength(Ensemble.ColumnSigmas, Ensemble.EnsembleSize*CCount-1+1);
    SetLength(Ensemble.TmpWeights, Ensemble.WCount-1+1);
    SetLength(Ensemble.TmpMeans, CCount-1+1);
    SetLength(Ensemble.TmpSigmas, CCount-1+1);
    SetLength(Ensemble.Neurons, NTotal-1+1);
    SetLength(Ensemble.DFDNET, NTotal-1+1);
    SetLength(Ensemble.SerializedMLP, Ensemble.SerializedLen-1+1);
    SetLength(Ensemble.Y, Ensemble.NOut-1+1);
    
    //
    // load data
    //
    I:=Offs;
    while I<=Offs+SSize-1 do
    begin
        Ensemble.StructInfo[I-Offs] := Round(RA[I]);
        Inc(I);
    end;
    Offs := Offs+SSize;
    APVMove(@Ensemble.Weights[0], 0, Ensemble.EnsembleSize*Ensemble.WCount-1, @RA[0], Offs, Offs+Ensemble.EnsembleSize*Ensemble.WCount-1);
    Offs := Offs+Ensemble.EnsembleSize*Ensemble.WCount;
    APVMove(@Ensemble.ColumnMeans[0], 0, Ensemble.EnsembleSize*CCount-1, @RA[0], Offs, Offs+Ensemble.EnsembleSize*CCount-1);
    Offs := Offs+Ensemble.EnsembleSize*CCount;
    APVMove(@Ensemble.ColumnSigmas[0], 0, Ensemble.EnsembleSize*CCount-1, @RA[0], Offs, Offs+Ensemble.EnsembleSize*CCount-1);
    Offs := Offs+Ensemble.EnsembleSize*CCount;
    APVMove(@Ensemble.SerializedMLP[0], 0, Ensemble.SerializedLen-1, @RA[0], Offs, Offs+Ensemble.SerializedLen-1);
    Offs := Offs+Ensemble.SerializedLen;
end;


(*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPERandomize(var Ensemble : MLPEnsemble);
var
    I : AlglibInteger;
begin
    I:=0;
    while I<=Ensemble.EnsembleSize*Ensemble.WCount-1 do
    begin
        Ensemble.Weights[I] := RandomReal-Double(0.5);
        Inc(I);
    end;
end;


(*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEProperties(const Ensemble : MLPEnsemble;
     var NIn : AlglibInteger;
     var NOut : AlglibInteger);
begin
    NIn := Ensemble.NIn;
    NOut := Ensemble.NOut;
end;


(*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
function MLPEIsSoftmax(const Ensemble : MLPEnsemble):Boolean;
begin
    Result := Ensemble.IsSoftmax;
end;


(*************************************************************************
Procesing

INPUT PARAMETERS:
    Ensemble-   neural networks ensemble
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NOut-1].

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEProcess(var Ensemble : MLPEnsemble;
     const X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    ES : AlglibInteger;
    WC : AlglibInteger;
    CC : AlglibInteger;
    V : Double;
begin
    ES := Ensemble.EnsembleSize;
    WC := Ensemble.WCount;
    if Ensemble.IsSoftmax then
    begin
        CC := Ensemble.NIn;
    end
    else
    begin
        CC := Ensemble.NIn+Ensemble.NOut;
    end;
    V := AP_Double(1)/ES;
    I:=0;
    while I<=Ensemble.NOut-1 do
    begin
        Y[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=ES-1 do
    begin
        APVMove(@Ensemble.TmpWeights[0], 0, WC-1, @Ensemble.Weights[0], I*WC, (I+1)*WC-1);
        APVMove(@Ensemble.TmpMeans[0], 0, CC-1, @Ensemble.ColumnMeans[0], I*CC, (I+1)*CC-1);
        APVMove(@Ensemble.TmpSigmas[0], 0, CC-1, @Ensemble.ColumnSigmas[0], I*CC, (I+1)*CC-1);
        MLPInternalProcessVector(Ensemble.StructInfo, Ensemble.TmpWeights, Ensemble.TmpMeans, Ensemble.TmpSigmas, Ensemble.Neurons, Ensemble.DFDNET, X, Ensemble.Y);
        APVAdd(@Y[0], 0, Ensemble.NOut-1, @Ensemble.Y[0], 0, Ensemble.NOut-1, V);
        Inc(I);
    end;
end;


(*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Works both for classifier betwork and for regression networks which
are used as classifiers.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
function MLPERelClsError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    MLPEAllErrors(Ensemble, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := RelCls;
end;


(*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if ensemble solves regression task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
function MLPEAvgCE(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    MLPEAllErrors(Ensemble, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := AvgCE;
end;


(*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for classification task
RMS error means error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
function MLPERMSError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    MLPEAllErrors(Ensemble, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := RMS;
end;


(*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
function MLPEAvgError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    MLPEAllErrors(Ensemble, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := Avg;
end;


(*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average relative error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
function MLPEAvgRelError(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    MLPEAllErrors(Ensemble, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := AvgRel;
end;


(*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
Modified Levenberg-Marquardt algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEBaggingLM(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var OOBErrors : MLPCVReport);
begin
    MLPEBaggingInternal(Ensemble, XY, NPoints, Decay, Restarts, Double(0.0), 0, True, Info, Rep, OOBErrors);
end;


(*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
L-BFGS algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEBaggingLBFGS(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var OOBErrors : MLPCVReport);
begin
    MLPEBaggingInternal(Ensemble, XY, NPoints, Decay, Restarts, WStep, MaxIts, False, Info, Rep, OOBErrors);
end;


(*************************************************************************
Training neural networks ensemble using early stopping.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  6, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPETrainES(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
var
    I : AlglibInteger;
    K : AlglibInteger;
    CCount : AlglibInteger;
    PCount : AlglibInteger;
    TrnXY : TReal2DArray;
    ValXY : TReal2DArray;
    TrnSize : AlglibInteger;
    ValSize : AlglibInteger;
    Network : MultiLayerPerceptron;
    TmpInfo : AlglibInteger;
    TmpRep : MLPReport;
begin
    if (NPoints<2) or (Restarts<1) or AP_FP_Less(Decay,0) then
    begin
        Info := -1;
        Exit;
    end;
    if Ensemble.IsSoftmax then
    begin
        I:=0;
        while I<=NPoints-1 do
        begin
            if (Round(XY[I,Ensemble.NIn])<0) or (Round(XY[I,Ensemble.NIn])>=Ensemble.NOut) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
    end;
    Info := 6;
    
    //
    // allocate
    //
    if Ensemble.IsSoftmax then
    begin
        CCount := Ensemble.NIn+1;
        PCount := Ensemble.NIn;
    end
    else
    begin
        CCount := Ensemble.NIn+Ensemble.NOut;
        PCount := Ensemble.NIn+Ensemble.NOut;
    end;
    SetLength(TrnXY, NPoints-1+1, CCount-1+1);
    SetLength(ValXY, NPoints-1+1, CCount-1+1);
    MLPUnserialize(Ensemble.SerializedMLP, Network);
    Rep.NGrad := 0;
    Rep.NHess := 0;
    Rep.NCholesky := 0;
    
    //
    // train networks
    //
    K:=0;
    while K<=Ensemble.EnsembleSize-1 do
    begin
        
        //
        // Split set
        //
        repeat
            TrnSize := 0;
            ValSize := 0;
            I:=0;
            while I<=NPoints-1 do
            begin
                if AP_FP_Less(RandomReal,Double(0.66)) then
                begin
                    
                    //
                    // Assign sample to training set
                    //
                    APVMove(@TrnXY[TrnSize][0], 0, CCount-1, @XY[I][0], 0, CCount-1);
                    TrnSize := TrnSize+1;
                end
                else
                begin
                    
                    //
                    // Assign sample to validation set
                    //
                    APVMove(@ValXY[ValSize][0], 0, CCount-1, @XY[I][0], 0, CCount-1);
                    ValSize := ValSize+1;
                end;
                Inc(I);
            end;
        until (TrnSize<>0) and (ValSize<>0);
        
        //
        // Train
        //
        MLPTrainES(Network, TrnXY, TrnSize, ValXY, ValSize, Decay, Restarts, TmpInfo, TmpRep);
        if TmpInfo<0 then
        begin
            Info := TmpInfo;
            Exit;
        end;
        
        //
        // save results
        //
        APVMove(@Ensemble.Weights[0], K*Ensemble.WCount, (K+1)*Ensemble.WCount-1, @Network.Weights[0], 0, Ensemble.WCount-1);
        APVMove(@Ensemble.ColumnMeans[0], K*PCount, (K+1)*PCount-1, @Network.ColumnMeans[0], 0, PCount-1);
        APVMove(@Ensemble.ColumnSigmas[0], K*PCount, (K+1)*PCount-1, @Network.ColumnSigmas[0], 0, PCount-1);
        Rep.NGrad := Rep.NGrad+TmpRep.NGrad;
        Rep.NHess := Rep.NHess+TmpRep.NHess;
        Rep.NCholesky := Rep.NCholesky+TmpRep.NCholesky;
        Inc(K);
    end;
end;


(*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEAllErrors(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     var RelCls : Double;
     var AvgCE : Double;
     var RMS : Double;
     var Avg : Double;
     var AvgRel : Double);
var
    I : AlglibInteger;
    Buf : TReal1DArray;
    WorkX : TReal1DArray;
    Y : TReal1DArray;
    DY : TReal1DArray;
begin
    SetLength(WorkX, Ensemble.NIn-1+1);
    SetLength(Y, Ensemble.NOut-1+1);
    if Ensemble.IsSoftmax then
    begin
        SetLength(DY, 0+1);
        DSErrAllocate(Ensemble.NOut, Buf);
    end
    else
    begin
        SetLength(DY, Ensemble.NOut-1+1);
        DSErrAllocate(-Ensemble.NOut, Buf);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@WorkX[0], 0, Ensemble.NIn-1, @XY[I][0], 0, Ensemble.NIn-1);
        MLPEProcess(Ensemble, WorkX, Y);
        if Ensemble.IsSoftmax then
        begin
            DY[0] := XY[I,Ensemble.NIn];
        end
        else
        begin
            APVMove(@DY[0], 0, Ensemble.NOut-1, @XY[I][0], Ensemble.NIn, Ensemble.NIn+Ensemble.NOut-1);
        end;
        DSErrAccumulate(Buf, Y, DY);
        Inc(I);
    end;
    DSErrFinish(Buf);
    RelCls := Buf[0];
    AvgCE := Buf[1];
    Rms := Buf[2];
    Avg := Buf[3];
    AvgRel := Buf[4];
end;


(*************************************************************************
Internal bagging subroutine.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPEBaggingInternal(var Ensemble : MLPEnsemble;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     LMAlgorithm : Boolean;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var OOBErrors : MLPCVReport);
var
    XYS : TReal2DArray;
    S : TBoolean1DArray;
    OOBBuf : TReal2DArray;
    OOBCntBuf : TInteger1DArray;
    X : TReal1DArray;
    Y : TReal1DArray;
    DY : TReal1DArray;
    DSBuf : TReal1DArray;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    CCnt : AlglibInteger;
    PCnt : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    TmpRep : MLPReport;
    Network : MultiLayerPerceptron;
begin
    
    //
    // Test for inputs
    //
    if  not LMAlgorithm and AP_FP_Eq(WStep,0) and (MaxIts=0) then
    begin
        Info := -8;
        Exit;
    end;
    if (NPoints<=0) or (Restarts<1) or AP_FP_Less(WStep,0) or (MaxIts<0) then
    begin
        Info := -1;
        Exit;
    end;
    if Ensemble.IsSoftmax then
    begin
        I:=0;
        while I<=NPoints-1 do
        begin
            if (Round(XY[I,Ensemble.NIn])<0) or (Round(XY[I,Ensemble.NIn])>=Ensemble.NOut) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
    end;
    
    //
    // allocate temporaries
    //
    Info := 2;
    Rep.NGrad := 0;
    Rep.NHess := 0;
    Rep.NCholesky := 0;
    OOBErrors.RelClsError := 0;
    OOBErrors.AvgCE := 0;
    OOBErrors.RMSError := 0;
    OOBErrors.AvgError := 0;
    OOBErrors.AvgRelError := 0;
    NIn := Ensemble.NIn;
    NOut := Ensemble.NOut;
    if Ensemble.IsSoftmax then
    begin
        CCnt := NIn+1;
        PCnt := NIn;
    end
    else
    begin
        CCnt := NIn+NOut;
        PCnt := NIn+NOut;
    end;
    SetLength(XYS, NPoints-1+1, CCnt-1+1);
    SetLength(S, NPoints-1+1);
    SetLength(OOBBuf, NPoints-1+1, NOut-1+1);
    SetLength(OOBCntBuf, NPoints-1+1);
    SetLength(X, NIn-1+1);
    SetLength(Y, NOut-1+1);
    if Ensemble.IsSoftmax then
    begin
        SetLength(DY, 0+1);
    end
    else
    begin
        SetLength(DY, NOut-1+1);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        J:=0;
        while J<=NOut-1 do
        begin
            OOBBuf[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        OOBCntBuf[I] := 0;
        Inc(I);
    end;
    MLPUnserialize(Ensemble.SerializedMLP, Network);
    
    //
    // main bagging cycle
    //
    K:=0;
    while K<=Ensemble.EnsembleSize-1 do
    begin
        
        //
        // prepare dataset
        //
        I:=0;
        while I<=NPoints-1 do
        begin
            S[I] := False;
            Inc(I);
        end;
        I:=0;
        while I<=NPoints-1 do
        begin
            J := RandomInteger(NPoints);
            S[J] := True;
            APVMove(@XYS[I][0], 0, CCnt-1, @XY[J][0], 0, CCnt-1);
            Inc(I);
        end;
        
        //
        // train
        //
        if LMAlgorithm then
        begin
            MLPTrainLM(Network, XYS, NPoints, Decay, Restarts, Info, TmpRep);
        end
        else
        begin
            MLPTrainLBFGS(Network, XYS, NPoints, Decay, Restarts, WStep, MaxIts, Info, TmpRep);
        end;
        if Info<0 then
        begin
            Exit;
        end;
        
        //
        // save results
        //
        Rep.NGrad := Rep.NGrad+TmpRep.NGrad;
        Rep.NHess := Rep.NHess+TmpRep.NHess;
        Rep.NCholesky := Rep.NCholesky+TmpRep.NCholesky;
        APVMove(@Ensemble.Weights[0], K*Ensemble.WCount, (K+1)*Ensemble.WCount-1, @Network.Weights[0], 0, Ensemble.WCount-1);
        APVMove(@Ensemble.ColumnMeans[0], K*PCnt, (K+1)*PCnt-1, @Network.ColumnMeans[0], 0, PCnt-1);
        APVMove(@Ensemble.ColumnSigmas[0], K*PCnt, (K+1)*PCnt-1, @Network.ColumnSigmas[0], 0, PCnt-1);
        
        //
        // OOB estimates
        //
        I:=0;
        while I<=NPoints-1 do
        begin
            if  not S[I] then
            begin
                APVMove(@X[0], 0, NIn-1, @XY[I][0], 0, NIn-1);
                MLPProcess(Network, X, Y);
                APVAdd(@OOBBuf[I][0], 0, NOut-1, @Y[0], 0, NOut-1);
                OOBCntBuf[I] := OOBCntBuf[I]+1;
            end;
            Inc(I);
        end;
        Inc(K);
    end;
    
    //
    // OOB estimates
    //
    if Ensemble.IsSoftmax then
    begin
        DSErrAllocate(NOut, DSBuf);
    end
    else
    begin
        DSErrAllocate(-NOut, DSBuf);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        if OOBCntBuf[I]<>0 then
        begin
            V := AP_Double(1)/OOBCntBuf[I];
            APVMove(@Y[0], 0, NOut-1, @OOBBuf[I][0], 0, NOut-1, V);
            if Ensemble.IsSoftmax then
            begin
                DY[0] := XY[I,NIn];
            end
            else
            begin
                APVMove(@DY[0], 0, NOut-1, @XY[I][0], NIn, NIn+NOut-1, V);
            end;
            DSErrAccumulate(DSBuf, Y, DY);
        end;
        Inc(I);
    end;
    DSErrFinish(DSBuf);
    OOBErrors.RelClsError := DSBuf[0];
    OOBErrors.AvgCE := DSBuf[1];
    OOBErrors.RMSError := DSBuf[2];
    OOBErrors.AvgError := DSBuf[3];
    OOBErrors.AvgRelError := DSBuf[4];
end;


end.