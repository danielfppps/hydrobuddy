unit HB_Main;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  ComCtrls, StdCtrls, Menus, ExtCtrls, Buttons, Grids, TAGraph, TASeries,
  hb_load_salts, Dbf, DB, Math, densesolver, hb_commercialnutrient, hb_comparison,
  hb_waterquality, hb_addweight, hb_insprecision, hb_stockanalysis,
  hb_persubstance, hb_datasetname, hb_analysis,
  hb_freedom, dbf_fields, hb_ratios,LCLIntf, Types,IniFiles;

const
  IniFile = 'settings.ini';

type

  { TForm1 }

  TForm1 = class(TForm)
    Button1:    TButton;
    Button10:   TButton;
    Button11:   TButton;
    Button12: TBitBtn;
    Button13: TBitBtn;
    Button14: TBitBtn;
    Button15: TButton;
    Button16: TButton;
    Button17: TButton;
    Button18: TButton;
    Button19: TBitBtn;
    Button2: TBitBtn;
    Button20: TButton;
    Button22: TButton;
    Button23: TButton;
    Button24: TButton;
    Button25: TButton;
    Button26: TBitBtn;
    Button3:    TButton;
    Button4:    TButton;
    Button5:    TButton;
    Button6:    TButton;
    Button7:    TButton;
    Button8:    TButton;
    Button9:    TButton;
    CheckBox1:  TCheckBox;
    CheckBox3: TCheckBox;
    CheckBox5: TCheckBox;
    ComboBox1:  TComboBox;
    ComboBox3: TComboBox;
    Edit1:      TEdit;
    Edit12:     TEdit;
    Edit13:     TEdit;
    Edit14:     TEdit;
    Edit15:     TEdit;
    Edit16:     TEdit;
    Edit2:     TEdit;
    Edit19:     TEdit;
    Edit18:     TEdit;
    Edit11:     TEdit;
    Edit17:     TEdit;
    Edit4:      TEdit;
    Edit3:      TEdit;
    Edit5:      TEdit;
    Edit6:      TEdit;
    Edit7:      TEdit;
    Edit8:      TEdit;
    Edit9:      TEdit;
    Edit10:      TEdit;
    Image1:     TImage;
    Image10:    TImage;
    Image2:     TImage;
    Image3:     TImage;
    Image4:     TImage;
    Image5:     TImage;
    Image6:     TImage;
    Image7:     TImage;
    Image8:     TImage;
    Image9:     TImage;
    Label1:     TLabel;
    Label4:     TLabel;
    Label3:     TLabel;
    Label5:     TLabel;
    Label6:     TLabel;
    Label7:     TLabel;
    Label8:     TLabel;
    Label10:     TLabel;
    Label11:     TLabel;
    Label12:    TLabel;
    Label14:    TLabel;
    Label15:    TLabel;
    Label13:    TLabel;
    Label16:    TLabel;
    Label9:    TLabel;
    Label2:    TLabel;
    Label17:    TLabel;
    Label18:    TLabel;
    Label19:    TLabel;
    Label20:    TLabel;
    Label21:    TLabel;
    Label22:    TLabel;
    Label23:    TLabel;
    Label24:    TLabel;
    Label25:    TLabel;
    Label27:    TLabel;
    Label28:    TLabel;
    Label29:    TLabel;
    Label30:    TLabel;
    Label31:    TLabel;
    Label32:    TLabel;
    Label33:    TLabel;
    Label34:    TLabel;
    Label26:    TLabel;
    LabelVersion: TLabel;
    MenuItem2: TMenuItem;
    PageControl2: TPageControl;
    Panel4:     TPanel;
    Panel5:     TPanel;
    Panel6:     TPanel;
    Panel7:     TPanel;
    Panel8: TPanel;
    Panel9: TPanel;
    PopupMenu1: TPopupMenu;
    RadioButton10: TRadioButton;
    RadioButton11: TRadioButton;
    RadioButton12: TRadioButton;
    RadioButton13: TRadioButton;
    RadioButton14: TRadioButton;
    RadioButton15: TRadioButton;
    RadioButton8: TRadioButton;
    RadioButton9: TRadioButton;
    RLabel1:    TLabel;
    RLabel2:   TLabel;
    RLabel3:    TLabel;
    RLabel4:    TLabel;
    RLabel5:    TLabel;
    RLabel6:    TLabel;
    RLabel7:    TLabel;
    RLabel8:    TLabel;
    RLabel9:    TLabel;
    RLabel10:    TLabel;
    RLabel12:   TLabel;
    RLabel13:   TLabel;
    RLabel14:   TLabel;
    RLabel15:   TLabel;
    RLabel16:   TLabel;
    RLabel11:   TLabel;
    PageControl1: TPageControl;
    Panel1:     TPanel;
    Panel2:     TPanel;
    Panel3:     TPanel;
    RadioButton1: TRadioButton;
    RadioButton2: TRadioButton;
    RadioButton3: TRadioButton;
    RadioButton4: TRadioButton;
    RadioButton5: TRadioButton;
    RadioButton6: TRadioButton;
    RadioButton7: TRadioButton;
    SaveDialog1: TSaveDialog;
    SelectDirectoryDialog1: TSelectDirectoryDialog;
    StringGrid1: TStringGrid;
    StringGrid2: TStringGrid;
    TabSheet1:  TTabSheet;
    TabSheet10: TTabSheet;
    TabSheet11: TTabSheet;
    TabSheet12: TTabSheet;
    TabSheet13: TTabSheet;
    TabSheet2:  TTabSheet;
    TabSheet3:  TTabSheet;
    TabSheet4:  TTabSheet;
    TabSheet5:  TTabSheet;
    TabSheet6:  TTabSheet;
    TabSheet7:  TTabSheet;
    TabSheet8:  TTabSheet;
    TabSheet9:  TTabSheet;
    procedure Button10Click(Sender: TObject);
    procedure Button11Click(Sender: TObject);
    procedure Button12Click(Sender: TObject);
    procedure Button13Click(Sender: TObject);
    procedure Button14Click(Sender: TObject);
    procedure Button15Click(Sender: TObject);
    procedure Button16Click(Sender: TObject);
    procedure Button17Click(Sender: TObject);
    procedure Button18Click(Sender: TObject);
    procedure Button19Click(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure Button20Click(Sender: TObject);
    procedure Button22Click(Sender: TObject);
    procedure Button23Click(Sender: TObject);
    procedure Button24Click(Sender: TObject);
    procedure Button25Click(Sender: TObject);
    procedure Button26Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button4Resize(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
    procedure Button7Click(Sender: TObject);
    procedure Button8Click(Sender: TObject);
    procedure Button9Click(Sender: TObject);
    procedure CheckBox1Change(Sender: TObject);
    procedure CheckBox3Change(Sender: TObject);
    procedure CheckBox4Change(Sender: TObject);
    procedure CheckBox5Change(Sender: TObject);
    procedure ComboBox1Change(Sender: TObject);
    procedure ComboBox1Select(Sender: TObject);
    procedure ComboBox3Change(Sender: TObject);
    procedure FormClose(Sender: TObject; var CloseAction: TCloseAction);
    procedure FormCreate(Sender: TObject);
    procedure FormWindowStateChange(Sender: TObject);
    procedure MenuItem2Click(Sender: TObject);
    procedure PageControl1Exit(Sender: TObject);
    procedure Panel6Click(Sender: TObject);
    procedure RadioButton10Change(Sender: TObject);
    procedure RadioButton11Change(Sender: TObject);
    procedure RadioButton12Change(Sender: TObject);
    procedure RadioButton13Change(Sender: TObject);
    procedure RadioButton1Change(Sender: TObject);
    procedure RadioButton2Change(Sender: TObject);
    procedure RadioButton3Change(Sender: TObject);
    procedure RadioButton4Change(Sender: TObject);
    procedure RadioButton5Change(Sender: TObject);
    procedure RadioButton6Change(Sender: TObject);
    procedure RadioButton7Change(Sender: TObject);
    procedure RadioButton8Change(Sender: TObject);
    procedure RadioButton9Change(Sender: TObject);
    procedure StringGrid2DrawCell(Sender: TObject; aCol, aRow: Integer;
      aRect: TRect; aState: TGridDrawState);
    procedure StringGrid2EditButtonClick(Sender: TObject);
    procedure StringGrid2EditingDone(Sender: TObject);
    procedure TabSheet1ContextPopup(Sender: TObject; MousePos: TPoint;
      var Handled: Boolean);
    procedure ToggleBox1Change(Sender: TObject);
  private
    { private declarations }
    procedure GridShowHint(Sender: TObject; HintInfo: PHintInfo);
  public
    { public declarations }
    const
    NAME_IDX     : integer = 0;
    FORMULA_IDX  : integer = 1;
    AMOUNT_IDX   : integer = 2;
    UNIT_IDX     : integer = 3;
    COST_IDX     : integer = 4;
    SOURCE_IDX   : integer = 5;


    zi : array[0..15] of double =
    (
      1,
      1,
      1,
      2,
      2,
      2,
      2,
      2,
      1,
      2,
      2,
      2,
      1,
      1,
      1,
      1) ;
    var
        water_quality_db,formulations_db,substances_db,substances_used_db: string;
        IsLiquid: array of array of double ;
        procedure weightFineTunning;
        procedure UpdateList;
        procedure UpdateComboBox;
        procedure setinivalues ;
        procedure cleanresults;
        function getratio(first: string; second: string; third: string; items: integer): string ;
        procedure getmolarmasses(var molar_mass: array of double) ;
        procedure getequivalents(var equivalents: array of integer);
        procedure replaceNullWithZeroes();
        procedure LoadValues();
        function round2(const Number: extended; const Places: longint): extended;
  end;

var
  Form1: TForm1;
  datasetname : String ;
  prev_conc :string ;
  selected_cell: string ;

implementation

{ TForm1 }

procedure TForm1.replaceNullWithZeroes();
var
  j: integer;
begin
    for j := 1 to 16 do
    begin
        if (FindComponent('Edit' + IntToStr(j)) as TEdit).Text = '' then (FindComponent('Edit' + IntToStr(j)) as TEdit).Text := '0';
    end;
end;

procedure TForm1.UpdateComboBox;
var
  MyDbf: TDbf;
begin

  HB_Main.Form1.ComboBox1.Items.Clear;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := formulations_db ;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.First;                  // moves to the first data

  while not MyDbf.EOF do
  begin
    HB_Main.Form1.ComboBox1.Items.Add(MyDbf.FieldByName('Name').AsString);
    MyDbf.Next;                                     // use .next here NOT .findnext!
  end;

  MyDbf.Close;

  MyDbf.Free;

end;

function TForm1.round2(const Number: extended; const Places: longint): extended;
var
  t: extended;
begin
  t      := power(10, places);
  round2 := round(Number * t) / t;
end;

procedure TForm1.UpdateList;
var
  MyDbf: TDbf;
  i:     integer;
  j:     integer;
begin

  Form2.ListBox1.Items.Clear;
  Form2.ListBox2.Items.Clear;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := substances_db;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.First;                  // moves to the first data

  while not MyDbf.EOF do
  begin
    hb_load_salts.Form2.ListBox1.Items.Add(MyDbf.FieldByName('Name').AsString);
    MyDbf.Next;                                     // use .next here NOT .findnext!
  end;

  MyDbf.Close;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := substances_used_db;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.First;                  // moves to the first data

  while not MyDbf.EOF do
  begin
    hb_load_salts.Form2.ListBox2.Items.Add(MyDbf.FieldByName('Name').AsString);
    MyDbf.Next;                                     // use .next here NOT .findnext!
  end;

  MyDbf.Close;

  MyDbf.Free;

  for i := 0 to Form2.ListBox2.Items.Count - 1 do

  begin

    j := 0;

    while j <= Form2.ListBox1.Items.Count - 1 do

    begin

      if (Form2.ListBox1.Items[j] = Form2.ListBox2.Items[i]) then

      begin
        Form2.ListBox1.Items.Delete(j);
        j := j + 1;
      end;

      j := j + 1;

    end;

  end;


  // sort listboxes
  hb_load_salts.Form2.ListBox2.Sorted := true ;
  hb_load_salts.Form2.ListBox1.Sorted := true ;

end;

procedure TForm1.Button1Click(Sender: TObject);
begin

  UpdateList;

  Form2.Visible := True;
end;

procedure TForm1.Button20Click(Sender: TObject);
var
  i: integer;
begin
    for i := 1 to 16 do (FindComponent('Edit' + IntToStr(i)) as TEdit).Text := StringGrid1.cells[1,i];
end;

procedure TForm1.Button22Click(Sender: TObject);
begin
  hb_freedom.Form12.Visible := true ;;
end;

procedure TForm1.Button23Click(Sender: TObject);
var
  MyDbf: TDbf;
  FieldDefs: TDbfFieldDefs;
begin
{
MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := 'substances_used.dbf';
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;
  FieldDefs.Add('Density', ftFloat, 0, False);
  FieldDefs.Add('IsLiquid',ftInteger, 0, False);
  MyDbf.RestructureTable(FieldDefs, True);
  MyDbf.Free;
  FieldDefs.Free;
  }

  {
  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := 'formulations.dbf';
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;
  FieldDefs.Add('Units', ftString, 80, False);
  MyDbf.RestructureTable(FieldDefs, True);
  MyDbf.Free;
  FieldDefs.Free;
  }

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := water_quality_db;
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;
  FieldDefs.Add('pH', ftFloat, 0, False);
  FieldDefs.Add('GH', ftFloat, 0, False);
  FieldDefs.Add('KH', ftFloat, 0, False);
  MyDbf.RestructureTable(FieldDefs, True);
  MyDbf.Free;
  FieldDefs.Free;


end;

procedure TForm1.Button24Click(Sender: TObject);
var
  MyDbf: TDbf;
  FieldDefs: TDbfFieldDefs;
begin


  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := substances_db;
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;

  ShowMessage(InttoStr(FieldDefs.Count) + ' '+ MyDbf.TableName) ;

  FieldDefs.Free ;
  MyDbf.Free ;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := substances_used_db;
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;

  ShowMessage(InttoStr(FieldDefs.Count) + ' '+ MyDbf.TableName) ;

  FieldDefs.Free ;
  MyDbf.Free ;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := formulations_db;
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;

  ShowMessage(InttoStr(FieldDefs.Count) + ' '+ MyDbf.TableName) ;

  FieldDefs.Free ;
  MyDbf.Free ;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := water_quality_db;
  MyDbf.Open;
  FieldDefs := TDbfFieldDefs.Create(nil);
  FieldDefs.Assign(MyDbf.DbfFieldDefs);
  MyDbf.Close;

  ShowMessage(InttoStr(FieldDefs.Count) + ' '+ MyDbf.TableName) ;

  FieldDefs.Free ;
  MyDbf.Free ;

  {

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := 'substances.dbf';
MyDbf.Open             ;
MyDbf.Active := true ;

while not MyDbf.EOF do
    begin

        MyDbf.Edit;

        MyDbf.FieldByName('Density').AsFloat := 0 ;

        MyDbf.Post ;

        MyDbf.next;                                   // use .next here NOT .findnext!
    end;


MyDbf.Close ;

MyDbf.Free ;

}

end;

procedure TForm1.Button25Click(Sender: TObject);
begin
  hb_ratios.Form14.Visible := true ;
end;

procedure TForm1.Button26Click(Sender: TObject);
begin
  // This is the paypal link for donations to support the software.
  // Please do not remove or modify.
  OpenURL('https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6YR6X5AAEGBGJ');

end;

procedure TForm1.weightFineTunning;
var
  MyDbf:  TDbf;
  i:      integer;
  j:      integer;
  varcount: integer;
  Volume: double;
  array_size: integer;
  answer: integer;
  arraysize: integer;
  nameSubstance: string;
  weight : double;
  report: DenseSolverLSReport;
  solutions: array of double;
  varnames: array of string;
  temp:   array of string;
  temp2:  array of double;
  vartargetvalue: array of double;
  degree_of_freedom: string;
  Result: array of double;
  test:   double;
  grossError: array[0..15] of double;
  instrumentalError: array[0..15] of double;
  elementInSolutionA: array[0..15] of double;
  elementInSolutionB: array[0..15] of double;
  volumeunit: string;
  waterquality: array of double;
  instrumental_error: array of double;
  all_element_contributions : array of array of double ;
  all_element_names : array of string;
  all_element_targets : array of double ;
  preloaded_weight: array of double ;
  upper:  double;
  lower:  double;
  weight_error: double;
  volume_error: double;
  weight_factor: double;
  ec_contribution: array of double;
  molar_mass: array[0..15] of double ;
  cost:   array of double;
  conc_factor: array of double ;
  equivalents: array[0..15] of integer ;
  predicted_ec: double;
  nameToCompare: string;
  ionic_strength : double;
  all_solids: boolean;
  mixContribution: array[0..15] of double;
  totalWeight: double;
begin
  replaceNullWithZeroes;

   // clear listbox to get rid of old solutions

  hb_ratios.Form14.StringGrid1.Clean;
  hb_ratios.Form14.StringGrid1.RowCount := 1 ;
  hb_persubstance.Form9.StringGrid1.Clean;
  hb_persubstance.Form9.StringGrid1.RowCount := 1 ;
  hb_stockanalysis.Form8.StringGrid1.Clean ;

  //initializing variables for mix label calculation
  all_solids := True;
  totalWeight := 0;

  // EC parameters

  SetLength(ec_contribution, 16);

  if RadioButton15.Checked = True then
  begin
    ec_contribution[0]  := 0.403791;
    ec_contribution[1]  := -0.110791;
    ec_contribution[2]  := 0.142525;
    ec_contribution[3]  := -0.358782;
    ec_contribution[4]  := -0.283549;
    ec_contribution[5]  := 0.276118;
    ec_contribution[6]  := 0;
    ec_contribution[7]  := 0;
    ec_contribution[8]  := 0;
    ec_contribution[9]  := 0;
    ec_contribution[10] := 0;
    ec_contribution[11] := 0;
    ec_contribution[12] := 0;
    ec_contribution[13] := 0;
    ec_contribution[14] := 0;
    ec_contribution[15] := -0.308967;
  end;

  if RadioButton14.Checked = True then
  begin
    ec_contribution[0]  := 71.46 ;
    ec_contribution[1]  := 73 ;
    ec_contribution[2]  := 57 ;
    ec_contribution[3]  := 106 ;
    ec_contribution[4]  := 119 ;
    ec_contribution[5]  := 160;
    ec_contribution[6]  := 108.0;
    ec_contribution[7]  := 0;
    ec_contribution[8]  := 0;
    ec_contribution[9]  := 0;
    ec_contribution[10] := 0;
    ec_contribution[11] := 50.01;
    ec_contribution[12] := 100 ;
    ec_contribution[13] := 76.35 ;
    ec_contribution[14] := 0;
    ec_contribution[15] := 73.5;
  end;

  // define water quality parameters

  SetLength(waterquality, 16);

  waterquality[0]  := StrtoFloat(hb_waterquality.Form6.Edit1.Text);
  waterquality[1]  := StrtoFloat(hb_waterquality.Form6.Edit2.Text);
  waterquality[2]  := StrtoFloat(hb_waterquality.Form6.Edit3.Text);
  waterquality[3]  := StrtoFloat(hb_waterquality.Form6.Edit4.Text);
  waterquality[4]  := StrtoFloat(hb_waterquality.Form6.Edit5.Text);
  waterquality[5]  := StrtoFloat(hb_waterquality.Form6.Edit6.Text);
  waterquality[6]  := StrtoFloat(hb_waterquality.Form6.Edit7.Text);
  waterquality[7]  := StrtoFloat(hb_waterquality.Form6.Edit8.Text);
  waterquality[8]  := StrtoFloat(hb_waterquality.Form6.Edit9.Text);
  waterquality[9]  := StrtoFloat(hb_waterquality.Form6.Edit10.Text);
  waterquality[10] := StrtoFloat(hb_waterquality.Form6.Edit11.Text);
  waterquality[11] := StrtoFloat(hb_waterquality.Form6.Edit12.Text);
  waterquality[12] := StrtoFloat(hb_waterquality.Form6.Edit13.Text);
  waterquality[13] := StrtoFloat(hb_waterquality.Form6.Edit14.Text);
  waterquality[14] := StrtoFloat(hb_waterquality.Form6.Edit15.Text);
  waterquality[15] := StrtoFloat(hb_waterquality.Form6.Edit16.Text);

  // disable input mix analysis button, we will enable it later if all substances are solids
  Button19.Enabled := False ;

  // set instrument precision values

  volume_error := StrtoFloat(hb_insprecision.Form7.Edit1.Text) * 0.001;

  weight_error := StrtoFloat(hb_insprecision.Form7.Edit2.Text);

  // set weight factor (g) or (oz)

  if RadioButton8.Checked then
    weight_factor := 1;

  if RadioButton9.Checked then
    weight_factor := 0.0352739619;

  // update list to get matrix size from used substances

  UpdateList;

  arraysize := Form2.ListBox2.Items.Count;

  // define volume taking into account radio box setting for units
  // the end input volume is always converted to cubic meters
  // which is the unit used within all calculations ppm = g/m3

  // cubic meter
  Volume := StrtoFloat(Edit18.Text);

  // liter conversion
  if Radiobutton1.Checked then
    Volume := Volume / 1000;

  // gallon conversion
  if Radiobutton2.Checked then
    Volume := Volume * 3.78541178 / 1000;

  // correct volume if using concentrated solutions.
  if RadioButton6.Checked then
  Volume := Volume * StrToFloat(Edit17.Text);

  array_size := Form2.ListBox2.Items.Count;

  varcount := 0;

  // set and define molar mass array (contains molar masses for elements)

  getmolarmasses(molar_mass) ;
  getequivalents(equivalents) ;

  // define concentration conversion factor array

  SetLength(conc_factor, 16) ;

  // set size of string grid for solutions

  StringGrid2.RowCount := arraysize + 1 ;


  // set values for the concentration factor according to selected units


  // conversion factor for ppm is 1 for all values (no conversion needed)
  if RadioButton10.Checked then

  begin

  prev_conc := 'ppm' ;

       for i:= 0 to 15 do

       begin

       conc_factor[i] := 1 ;

       end ;

  end ;

  // conversion factor for M (mol/L)
 if RadioButton11.Checked then

  begin

       for i:= 0 to 15 do

       begin

       conc_factor[i] := 1000*molar_mass[i] ;

       end ;

  end ;

 // conversion factor for mM (mmol/L)
if RadioButton12.Checked then

 begin

      for i:= 0 to 15 do

      begin

      conc_factor[i] := molar_mass[i] ;

      end ;

 end;

// conversion factor for mN (meq/L)
if RadioButton13.Checked then

 begin

      for i:= 0 to 15 do

      begin

      conc_factor[i] := molar_mass[i]/equivalents[i] ;

      end ;

 end ;

    varnames := nil ;
    Result := nil ;

    SetLength(varnames, 16);
    SetLength(Result, 16);
    SetLength(all_element_targets, 16);

    for j := 1 to 16 do
    begin

      // load all element names (this time we don't need to discriminate as
      // we simply calculate for everyone
      varnames[j - 1] := (FindComponent('Label' + IntToStr(j)) as TLabel).Caption;
      all_element_targets[j-1] := StrToFloat((FindComponent('Edit' + IntToStr(j)) as TEdit).Text);
      StringGrid1.Cells[2, j] := '0';
      StringGrid1.Cells[3, j] := '0';
      grossError[j-1] := 0;
      instrumentalError[j-1] := 0;
      elementInSolutionA[j-1] := 0;
      elementInSolutionB[j-1] := 0;
      mixContribution[j-1] := 0;
    end;

    // load the database in order to get the weights and find the resulting ppm values

    MyDbf := TDbf.Create(nil);
    MyDbf.FilePathFull := '';
    MyDbf.TableName := substances_used_db;
    MyDbf.Open;
    MyDbf.Active := True;

    for i := 1 to StringGrid2.RowCount - 1 do
    begin

    MyDbf.First;                  // moves to the first data

    while not MyDbf.EOF do
    begin

      nameSubstance := MyDbf.FieldByName('Name').AsString;
      weight := StrToFloat(StringGrid2.Cells[AMOUNT_IDX, i]);
      nameToCompare := StringGrid2.Cells[NAME_IDX, i] ;

      if  MyDbf.FieldByName('isLiquid').AsBoolean = True then all_solids := False;

      If RadioButton6.Checked then
      nameToCompare := Copy(nameToCompare, 5, Length(nameToCompare));

      if nameSubstance = nameToCompare then
      begin

      for j := 0 to 15 do

      begin

        if ((1 / weight_factor) * weight *
          0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
          'Purity').AsFloat / Volume > 0) then

        begin

          upper := ((StrToFloat(StringGrid2.Cells[AMOUNT_IDX, i])+weight_error)) /
            (Volume - volume_error) ;
          lower := ((StrToFloat(StringGrid2.Cells[AMOUNT_IDX, i])-weight_error)) /
            (Volume + volume_error) ;

          if (RadioButton7.Checked = true) and (all_element_targets[j] <> 0) then
          instrumentalError[j] := ((upper - lower) * (100 / all_element_targets[j])) / 2 + instrumentalError[j];

          if (RadioButton6.Checked = true) and (all_element_targets[j] <> 0) then
          instrumentalError[j] := ((upper - lower) * (100 / all_element_targets[j])) /
          (2 * StrtoFloat(Edit17.Text)) + instrumentalError[j];

          hb_persubstance.Form9.StringGrid1.RowCount := hb_persubstance.Form9.StringGrid1.RowCount + 1 ;

          hb_persubstance.Form9.StringGrid1.Cells[0, hb_persubstance.Form9.StringGrid1.RowCount-1] := (nameSubstance);
          hb_persubstance.Form9.StringGrid1.Cells[1, hb_persubstance.Form9.StringGrid1.RowCount-1] :=(varnames[j]);

          if prev_conc = 'ppm' then

          hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1] :=(
            FloattoStr(round2((1 / weight_factor) * weight *
            0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
            'Purity').AsFloat / Volume, 3)));

          if prev_conc <> 'ppm' then

          hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1] :=(
            FloattoStrF((1/conc_factor[j])*(1 / weight_factor) * weight *
            0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
            'Purity').AsFloat / Volume, ffExponent, 4, 2));

          Result[j] := Result[j] + StrToFloat(hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1]);

          if (StringGrid2.Cells[NAME_IDX, i][1] = 'A') and (RadioButton6.Checked) then
          elementInSolutionA[j] := StrToFloat(hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1])+ elementInSolutionA[j] ;

          if (StringGrid2.Cells[NAME_IDX, i][1] = 'B') and (RadioButton6.Checked) then
          elementInSolutionB[j] := StrToFloat(hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1]) + elementInSolutionB[j] ;

          //ShowMessage( nameSubstance + ' ' + varnames[j] + ' ' + FloatToStr(weight * 0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName('Purity').AsFloat));
          mixContribution[j] := mixContribution[j] + weight * 0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName('Purity').AsFloat ;

        end;

      end;

      StringGrid2.Cells[COST_IDX,i] := (FloattoStr(
        round2(weight * MyDbf.FieldByName('Cost').AsFloat * 0.001 *
        (1 / weight_factor), 1)));

      end;

      MyDbf.Next;                                     // use .next here NOT .findnext!
    end;

    end;

    MyDbf.Close;
    MyDbf.Free;

    for i := 1 to StringGrid2.RowCount - 1 do totalWeight := totalWeight+StrToFloat(StringGrid2.Cells[AMOUNT_IDX, i]);

    // save results in main page as well
    for j := 1 to 16 do
    begin

        if all_element_targets[j-1] <> 0 then
        grossError[j-1] :=  (Result[j-1] * 100 / all_element_targets[j-1]) - 100 ;

        (FindComponent('RLabel' + IntToStr(j)) as TLabel).Caption := StringGrid1.Cells[1, i];
        StringGrid1.Cells[3, j] := '+/- ' + FloatToStr(Round(instrumentalError[j-1]*10)/10) + '%' ;
        StringGrid1.Cells[2, j] := FloatToStr(Round(grossError[j-1]*10)/10) + '%' ;

        // save the mix composition here
        hb_analysis.Form11.StringGrid1.Cells[1,j] := FloatToStr(round2(100*mixContribution[j-1]/totalWeight,3));
        if hb_analysis.Form11.StringGrid1.Cells[0,j] = 'K2O' then hb_analysis.Form11.StringGrid1.Cells[1,j] := FloatToStr(round2(1.2047*100*mixContribution[j-1]/totalWeight,3));
        if hb_analysis.Form11.StringGrid1.Cells[0,j] = 'P2O5' then hb_analysis.Form11.StringGrid1.Cells[1,j] := FloatToStr(round2(2.290*100*mixContribution[j-1]/totalWeight,3));

    end;

    if all_solids then Button19.Enabled := True;

    for i := 1 to 16 do
    begin

      for j := 0 to 15 do

      begin

        if (FindComponent('Label' + IntToStr(i)) as TLabel).Caption = varnames[j] then

        begin


          StringGrid1.Cells[0,i] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption;

          if prev_conc = 'ppm' then

          begin


          StringGrid1.Cells[1,i] := (FloatToStr(round2(Result[j] + waterquality[i - 1], 3)));

          (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption :=
            FloatToStr(round2(Result[j] + waterquality[i - 1], 3));

          end ;

          if prev_conc <> 'ppm' then

          begin

          StringGrid1.Cells[1,i] :=(FloatToStrF((1/conc_factor[j])*(Result[j] + waterquality[i - 1]), ffExponent, 4, 2));

          (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption :=
            FloatToStrF((1/conc_factor[j])*(Result[j] + waterquality[i - 1]), ffExponent,4, 2);

          end;

        end;

      end;

    end;

    // calculation of EC by empirical model
    if RadioButton15.Checked = True then
    begin
        predicted_ec := 0;
        for i := 1 to 16 do
        begin
            predicted_ec := conc_factor[i-1]*StrtoFloat(
              (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption) *
              ec_contribution[i - 1] + predicted_ec;
        end;
        predicted_ec := round2(predicted_ec+0.39661671, 3);
    end;

    // Calculation of EC by LMCv2 model
    if RadioButton14.Checked = True then
    begin

        // calculate ionic strength used for conductivity model
        ionic_strength := 0;
        for i := 1 to 16 do ionic_strength := zi[i-1]*zi[i-1]*(StrtoFloat((FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption) /(1000*molar_mass[i-1])) + ionic_strength;

        predicted_ec := 0;
        for i := 1 to 16 do
        begin
            predicted_ec := conc_factor[i-1]
                            * (StrtoFloat((FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption)/(1000*molar_mass[i-1]))
                            * ec_contribution[i - 1]
                            * exp(-0.7025187*sqrt(ionic_strength)*power(zi[i-1],1.5))
                            + predicted_ec;
        end;
        predicted_ec := round2(predicted_ec, 3);
    end;

    Panel6.Caption := 'EC=' + FloattoStr(predicted_ec) + ' mS/cm';

    Button10.Enabled := True;

    //determine volume unit for description label

        if RadioButton2.Checked then
          volumeunit := 'gallons';

        if RadioButton3.Checked then
          volumeunit := 'cubic meters';

        if RadioButton1.Checked then
          volumeunit := 'liters';

          Label20.Caption := 'Values calculated for the preparation of ' +
            Edit18.Text + ' ' + volumeunit;

  // total cost calculation

  test := 0;

  for i := 0 to StringGrid2.RowCount - 2 do

  begin

    test := StrtoFloat(StringGrid2.Cells[COST_IDX,i+1]) + test;

  end;

  Label18.Caption := ('Total Cost is ' + FloattoStr(round2(test, 1)));

  // stock solution analysis

  if (RadioButton6.Checked) then
  begin

  for i := 0 to 15 do
  begin

        //ShowMessage(FloatToStr(elementInSolutionA[i]));
        elementInSolutionA[i] := 100 * elementInSolutionA[i] * StrtoFloat(Edit17.Text) / (1000*1000);

        if (varnames[i] <> 'P') and (varnames[i] <> 'K') then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[0, i+ 1] := (varnames[i]);
          hb_stockanalysis.Form8.StringGrid1.Cells[1, i+ 1] := (FloattoStr(round2(elementInSolutionA[i], 5)));
        end;

        if varnames[i] = 'P' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[0, i+ 1] := ('P2O5');
          hb_stockanalysis.Form8.StringGrid1.Cells[1, i+ 1] := (
            FloattoStr(round2(elementInSolutionA[i] * 2.2915, 3)));
        end;


        if varnames[i] = 'K' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[0, i+ 1] := ('K2O');
          hb_stockanalysis.Form8.StringGrid1.Cells[1, i+ 1] := (
            FloattoStr(round2(elementInSolutionA[i] * 1.2047, 3)));
        end;


  end;

  for i := 0 to 15 do
  begin

        elementInSolutionB[i] := 100 * elementInSolutionB[i] * StrtoFloat(Edit17.Text) / (1000*1000);

        if (varnames[i] <> 'P') and (varnames[i] <> 'K') then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[2, i+ 1] := (FloattoStr(round2(elementInSolutionB[i], 5)));
        end;

        if varnames[i] = 'P' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[2, i+ 1] := (
            FloattoStr(round2(elementInSolutionB[i] * 2.2915, 3)));
        end;


        if varnames[i] = 'K' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[2, i+ 1] := (
            FloattoStr(round2(elementInSolutionB[i] * 1.2047, 3)));
        end;

  end;


  end;

   // post ratios based on results posted on listboxes above

    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N: P: K') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N', 'P', 'K', 3)) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N: P2O5: K2O') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N', 'P2O5', 'K2O', 3)) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N: K') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N', 'K', 'K', 2) ) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N (NO3-): N (NH4+)') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N (NO3-)', 'N (NH4+)', 'K', 2) ) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('Ca: Mg') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('Ca', 'Mg', 'K', 2) ) ;

  end;

procedure TForm1.Button10Click(Sender: TObject);
var
  MyDbf: TDbf;
  i:     integer;
  used_string : string ;
begin

  if RadioButton4.Checked then

  begin

    MyDbf := TDbf.Create(nil);
    MyDbf.FilePathFull := '';
    MyDbf.TableName := substances_used_db;
    MyDbf.Open;
    MyDbf.Active := True;

    for i := 1 to StringGrid2.RowCount - 1 do

    begin

    used_string := StringGrid2.Cells[NAME_IDX,i] ;

    if RadioButton6.checked then
    delete (used_string,1,4);


      MyDbf.Filter := 'Name=' + QuotedStr(used_string);

      MyDbf.Filtered := True;       // This selects the filtered set
      MyDbf.First;

      MyDbf.Edit;

      MyDbf.FieldByName('Weight').AsFloat := StrtoFloat(StringGrid2.Cells[AMOUNT_IDX,i]);

      if RadioButton6.checked then
      MyDbf.FieldByName('Weight').AsFloat := StrtoFloat(StringGrid2.Cells[AMOUNT_IDX,i])/StrToFloat(Edit17.Text);

      MyDbf.Post;

    end;

    MyDbf.Close;

    MyDbf.Free;

    ShowMessage('Weights have been successfully copied to Database');

  end;

  if RadioButton5.Checked then

  begin

    for i := 1 to 16 do

    begin

      (FindComponent('Edit' + IntToStr(i)) as TEdit).Text :=
        (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption;

    end;

  end;

end;

procedure TForm1.Button11Click(Sender: TObject);
begin

  PageControl1.ActivePage := TabSheet1;

end;

procedure TForm1.Button12Click(Sender: TObject);
begin

  hb_stockanalysis.Form8.Visible := True;

end;

procedure TForm1.Button13Click(Sender: TObject);
var
  i: integer;
begin

  with TStringList.Create do

    try

      Add(Label20.Caption);

      Add(' , , , , ');
      Add('Name, Formula, Amount, Units, Cost');

      for i := 0 to StringGrid2.RowCount - 1 do

      begin

        Add(StringGrid2.Cells[NAME_IDX,i]    + ',' +
            StringGrid2.Cells[FORMULA_IDX,i] + ',' +
            StringGrid2.Cells[AMOUNT_IDX,i]  + ',' +
            StringGrid2.Cells[UNIT_IDX,i]    + ',' +
            StringGrid2.Cells[COST_IDX,i]);

      end;

      Add(' , , , , ');
      Add('Element, Results(ppm), GE, IE, Water (ppm)');

      for i := 0 to StringGrid1.RowCount - 1 do

      begin

        Add(StringGrid1.Cells[0,i] + ',' +
            StringGrid1.Cells[1,i] + ',' +
            StringGrid1.Cells[2,i] + ',' +
            StringGrid1.Cells[3,i] + ',' +
            StringGrid1.Cells[4,i]);

      end;

      Add(' , , , ,');
      Add(Panel6.Caption);

      if SaveDialog1.Execute then
        SaveToFile(SaveDialog1.filename);

    finally

      Free;

    end;

end;

procedure TForm1.Button14Click(Sender: TObject);
begin

  hb_persubstance.Form9.Visible := True;

end;

procedure TForm1.Button15Click(Sender: TObject);
begin
  Edit1.Text  := '0';
  Edit3.Text  := '0';
  Edit4.Text  := '0';
  Edit5.Text  := '0';
  Edit6.Text  := '0';
  Edit7.Text  := '0';
  Edit8.Text  := '0';
  Edit10.Text  := '0';
  Edit9.Text  := '0';
  Edit12.Text := '0';
  Edit13.Text := '0';
  Edit14.Text := '0';
  Edit11.Text := '0';
  Edit15.Text := '0';
  Edit16.Text := '0';
  Edit2.Text := '0';
end;

procedure TForm1.Button16Click(Sender: TObject);
var
  i: integer;
  current_weight: double;
begin

  If StringGrid2.RowCount < 2 then exit;

  for i := 1 to StringGrid2.RowCount-1 do
  begin
      current_weight := StrToFloat(StringGrid2.Cells[AMOUNT_IDX, i]);
      StringGrid2.Cells[AMOUNT_IDX, i] := FloatToStr(round2(current_weight*1.05, 4));
  end;

  weightFineTunning;

end;

procedure TForm1.Button17Click(Sender: TObject);
var
  i: integer;
  current_weight: double;
begin

  If StringGrid2.RowCount < 2 then exit;

  for i := 1 to StringGrid2.RowCount-1 do
  begin
      current_weight := StrToFloat(StringGrid2.Cells[AMOUNT_IDX, i]);
      StringGrid2.Cells[AMOUNT_IDX, i] := FloatToStr(round2(current_weight*0.95, 4));
  end;

  weightFineTunning;

end;

procedure TForm1.Button18Click(Sender: TObject);
var
  MyDbf: TDbf;
  i:     integer;
  selected_item: integer;
begin

  Edit19.Text := 'DEFAULT' ;

  if ComboBox1.Items.Count = 0 then
    Exit;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := formulations_db;
  MyDbf.Open;
  MyDbf.Active := True;
  MyDbf.Filter := 'Name=' + QuotedStr('DEFAULT');

  MyDbf.Filtered := True;       // This selects the filtered set
  MyDbf.First;                  // moves the the first filtered data
  ComboBox1.Items.Delete(ComboBox1.ItemIndex);
  MyDbf.Delete;
  MyDbf.Close;
  MyDbf.Free;

  Button4Click(Sender);
end;

procedure TForm1.Button19Click(Sender: TObject);
begin
      hb_analysis.Form11.visible := True ;
end;

procedure TForm1.getmolarmasses(var molar_mass: array of double) ;
begin

  molar_mass[0] :=  14.007 ;
  molar_mass[1] :=  39.098 ;
  molar_mass[2] :=  30.974 ;
  molar_mass[3] :=  24.305 ;
  molar_mass[4] :=  40.078 ;
  molar_mass[5] :=  32.066 ;
  molar_mass[6] :=  55.845 ;
  molar_mass[7] :=  65.409 ;
  molar_mass[8] :=  10.811 ;
  molar_mass[9] :=  63.546 ;
  molar_mass[10] := 95.94  ;
  molar_mass[11] := 22.990  ;
  molar_mass[12] := 28.086  ;
  molar_mass[13] := 35.453  ;
  molar_mass[14] := 54.938  ;
  molar_mass[15] := 14.007  ;

end ;

procedure TForm1.getequivalents(var equivalents: array of integer);
begin

  equivalents[0] := 1 ;
  equivalents[1] := 1 ;
  equivalents[2] := 3 ;
  equivalents[3] := 2 ;
  equivalents[4] := 2 ;
  equivalents[5] := 2 ;
  equivalents[6] := 2 ;
  equivalents[7] := 2 ;
  equivalents[8] := 1 ;
  equivalents[9] := 2 ;
  equivalents[10] :=2  ;
  equivalents[11] :=1  ;
  equivalents[12] := 2  ;
  equivalents[13] := 2  ;
  equivalents[14] := 2  ;
  equivalents[15] := 1  ;

end;

function TForm1.getratio(first: string; second: string; third: string; items: integer): string ;
var
i: integer ;
j: integer ;
ratio: array[0..2] of double ;
temp : array [0..3] of double ;
names: array[0..3] of string ;
special: array[0..2] of integer ;
begin

hb_ratios.Form14.StringGrid1.RowCount := hb_ratios.Form14.StringGrid1.RowCount + 1 ;

temp[0] := 0 ;
temp[1] := 0 ;
temp[2] := 0 ;
temp[3] := 0 ;

names[0] := first ;
names[1] := second ;
names[2] := third ;

special[0] := -1 ;
special[1] := -1 ;
special[2] := -1 ;

if (items > 3) or
   (items < 2) then
   Result := 'no result' ;

for i := 0 to 3 do

         begin

              if names[i] = 'N' then

               begin

               names[i] := 'N (NO3-)' ;
               names[3]  := 'N (NH4+)' ;
               special[0] := i ;

               end ;

              if names[i] = 'P2O5' then

               begin

               names[i] := 'P' ;
               special[1] := i ;

               end ;

              if names[i] = 'K2O' then

               begin

               names[i] := 'K' ;
               special[2] := i ;

               end ;

         end;

for i:= 1 to 16 do

    begin

    for j := 0 to 3 do

        begin

         if StringGrid1.Cells[0,i] = names[j] then
          temp[j] := StrToFloat(StringGrid1.Cells[1,i]) ;

        end;

    end ;

if (temp[0] > 0) and
   (temp[1] > 0) and
   (temp[2] > 0) then

   begin

if (special[0] <> -1) and
   (names[special[0]] = 'N (NO3-)') then

   begin

   temp[special[0]] := temp[3]+temp[special[0]] ;

   end ;

if (special[1] <> -1) and
      (names[special[1]] = 'P') then

   begin

   temp[special[1]] := temp[special[1]]*2.2915 ;

   end ;

   if (special[2] <> -1) and
         (names[special[2]] = 'K') then

      begin

      temp[special[2]] := temp[special[2]]*1.2047 ;

      end ;

ratio[0] := 1 ;
ratio[1] := round2(temp[1]/temp[0],2) ;
ratio[2] := round2(temp[2]/temp[0],2) ;

if Items = 2 then
Result := FloatToStr(ratio[0]) + ': ' + FloatToStr(ratio[1]) ;

if Items = 3 then
Result := FloatToStr(ratio[0]) + ': ' + FloatToStr(ratio[1]) + ': ' + FloatToStr(ratio[2]) ;

   end;


end;

procedure TForm1.Button2Click(Sender: TObject);
var
  MyDbf:  TDbf;
  i:      integer;
  j:      integer;
  k:      integer;
  n:      integer ;
  varcount: integer;
  Volume: double;
  problem_matrix_left: array of array of double;
  problem_matrix_right: array of double;
  array_size: integer;
  answer: integer;
  arraysize: integer;
  name_array: array of array of string;
  report: DenseSolverLSReport;
  solutions: array of double;
  varnames: array of string;
  temp:   array of string;
  temp2:  array of double;
  vartargetvalue: array of double;
  degree_of_freedom: string;
  Result: array of double;
  test:   double;
  gross_error: array of double;
  instrumental_erro: array of double;
  ConcTypeArray: array of string;
  volumeunit: string;
  waterquality: array of double;
  instrumental_error: array of double;
  all_element_contributions : array of array of double ;
  all_element_names : array of string;
  all_element_targets : array of double ;
  preloaded_weight: array of double ;
  upper:  double;
  temp1: double ;
  temp3: double ;
  lower:  double;
  weight_error: double;
  volume_error: double;
  weight_factor: double;
  ec_contribution: array of double;
  molar_mass: array[0..15] of double ;
  cost:   array of double;
  conc_factor: array of double ;
  equivalents: array[0..15] of integer ;
  predicted_ec: double;
  mass_unit: string;
  ionic_strength: double;
  all_solids: boolean;
  mixContribution: array[0..15] of double;
  totalWeight: double;
begin

  //deal with null
  replaceNullWithZeroes();

  // clear listbox to get rid of old solutions

  StringGrid1.Clean;
  StringGrid2.Clean;
  StringGrid2.RowCount:= 1;
  StringGrid2.ColWidths[NAME_IDX] := 223;
  StringGrid2.ColWidths[FORMULA_IDX] := 180;
  StringGrid2.ColWidths[AMOUNT_IDX] := 170;
  StringGrid2.ColWidths[COST_IDX] := 110;
  hb_ratios.Form14.StringGrid1.Clean;
  hb_ratios.Form14.StringGrid1.RowCount := 1 ;
  hb_persubstance.Form9.StringGrid1.Clean;
  hb_persubstance.Form9.StringGrid1.RowCount := 1 ;
  hb_stockanalysis.Form8.StringGrid1.Clean ;

  //initializing variables for mix label calculation
  all_solids := True;
  totalWeight := 0;
  Button19.Enabled := False;
  for j:=0 to 15 do mixContribution[j] := 0;

  // EC parameters

  SetLength(ec_contribution, 16);

  if RadioButton15.Checked = True then
  begin
    ec_contribution[0]  := 0.403791;
    ec_contribution[1] := -0.308967;
    ec_contribution[2]  := 0.142525;
    ec_contribution[3]  := -0.110791;
    ec_contribution[4]  := -0.358782;
    ec_contribution[5]  := -0.283549;
    ec_contribution[6]  := 0.276118;
    ec_contribution[7]  := 0;
    ec_contribution[8]  := 0;
    ec_contribution[9]  := 0;
    ec_contribution[10]  := 0;
    ec_contribution[11] := 0;
    ec_contribution[12] := 0;
    ec_contribution[13] := 0;
    ec_contribution[14] := 0;
    ec_contribution[15] := 0;

  end;

  if RadioButton14.Checked = True then
  begin
    ec_contribution[0]  := 71.46 ;
    ec_contribution[1] := 73.5;
    ec_contribution[2]  := 57 ;
    ec_contribution[3]  := 73 ;
    ec_contribution[4]  := 106 ;
    ec_contribution[5]  := 119 ;
    ec_contribution[6]  := 160 ;
    ec_contribution[7]  := 108.0;
    ec_contribution[8]  := 0;
    ec_contribution[9]  := 0;
    ec_contribution[10]  := 0;
    ec_contribution[11] := 0;
    ec_contribution[12] := 0; // at the pH used in hydroponics, silicon does not conduct
    ec_contribution[13] := 76.35 ;
    ec_contribution[14] := 50.01 ;
    ec_contribution[15] := 0;
  end;


  // define water quality parameters
  SetLength(waterquality, 16);

  waterquality[0]  := StrtoFloat(hb_waterquality.Form6.Edit1.Text);
  waterquality[1]  := StrtoFloat(hb_waterquality.Form6.Edit2.Text);
  waterquality[2]  := StrtoFloat(hb_waterquality.Form6.Edit3.Text);
  waterquality[3]  := StrtoFloat(hb_waterquality.Form6.Edit4.Text);
  waterquality[4]  := StrtoFloat(hb_waterquality.Form6.Edit5.Text);
  waterquality[5]  := StrtoFloat(hb_waterquality.Form6.Edit6.Text);
  waterquality[6]  := StrtoFloat(hb_waterquality.Form6.Edit7.Text);
  waterquality[7]  := StrtoFloat(hb_waterquality.Form6.Edit8.Text);
  waterquality[8]  := StrtoFloat(hb_waterquality.Form6.Edit9.Text);
  waterquality[9]  := StrtoFloat(hb_waterquality.Form6.Edit10.Text);
  waterquality[10] := StrtoFloat(hb_waterquality.Form6.Edit11.Text);
  waterquality[11] := StrtoFloat(hb_waterquality.Form6.Edit12.Text);
  waterquality[12] := StrtoFloat(hb_waterquality.Form6.Edit13.Text);
  waterquality[13] := StrtoFloat(hb_waterquality.Form6.Edit14.Text);
  waterquality[14] := StrtoFloat(hb_waterquality.Form6.Edit15.Text);
  waterquality[15] := StrtoFloat(hb_waterquality.Form6.Edit16.Text);

  // set instrument precision values

  volume_error := StrtoFloat(hb_insprecision.Form7.Edit1.Text) * 0.001;

  weight_error := StrtoFloat(hb_insprecision.Form7.Edit2.Text);

  // set weight factor (g) or (oz)

  if RadioButton8.Checked then
  begin
    weight_factor := 1;
    mass_unit := 'g';
  end;

  if RadioButton9.Checked then
  begin
    weight_factor := 0.0352739619;
    mass_unit := 'oz';
  end;

  // update list to get matrix size from used substances

  UpdateList;

  arraysize := Form2.ListBox2.Items.Count;

  // define volume taking into account radio box setting for units
  // the end input volume is always converted to cubic meters
  // which is the unit used within all calculations ppm = g/m3

  // cubic meter
  Volume := StrtoFloat(Edit18.Text);

  // liter conversion
  if Radiobutton1.Checked then
    Volume := Volume / 1000;

  // gallon conversion
  if Radiobutton2.Checked then
    Volume := Volume * 3.78541178 / 1000;

  array_size := Form2.ListBox2.Items.Count;

  varcount := 0;

  // set and define molar mass array (contains molar masses for elements)

  getmolarmasses(molar_mass) ;
  getequivalents(equivalents) ;

  // define concentration conversion factor array

  SetLength(conc_factor, 16) ;

  // set size of string grid for solutions

  StringGrid2.RowCount := arraysize + 1 ;


  // set values for the concentration factor according to selected units


  // conversion factor for ppm is 1 for all values (no conversion needed)
  if RadioButton10.Checked then

  begin

  prev_conc := 'ppm' ;

       for i:= 0 to 15 do

       begin

       conc_factor[i] := 1 ;

       end ;

  end ;

  // conversion factor for M (mol/L)
 if RadioButton11.Checked then

  begin

       for i:= 0 to 15 do

       begin

       conc_factor[i] := 1000*molar_mass[i] ;

       end ;

  end ;

 // conversion factor for mM (mmol/L)
if RadioButton12.Checked then

 begin

      for i:= 0 to 15 do

      begin

      conc_factor[i] := molar_mass[i] ;

      end ;

 end;

// conversion factor for mN (meq/L)
if RadioButton13.Checked then

 begin

      for i:= 0 to 15 do

      begin

      conc_factor[i] := molar_mass[i]/equivalents[i] ;

      end ;

 end ;

  // This IF statement controls the type of calculation run
  // according to the radio buttons which control it

  // VERY IMPORTANT IF STATEMENT //
  if RadioButton4.Checked then
    // VERY IMPORTANT IF STATEMENT //
  begin

  // define arrays used to solve the linear equation system

    problem_matrix_left := nil ;

    problem_matrix_right := nil ;

    solutions := nil ;

    name_array := nil ;

    ConcTypeArray := nil ;

    cost := nil ;

    all_element_contributions := nil ;

    all_element_names := nil ;

    all_element_targets := nil ;

    SetLength(preloaded_weight, arraysize) ;

    SetLength(solutions, arraysize);

    SetLength(name_array, arraysize, 3);

    SetLength(ConcTypeArray, arraysize);

    SetLength(cost, arraysize);

    SetLength(IsLiquid, 2, arraysize) ;

    SetLength(all_element_contributions, 16, arraysize);

    SetLength(all_element_names, 16);

    SetLength(all_element_targets, 16);

     // set array to 0

    for i := 0 to arraysize - 1 do
    begin
      for j := 0 to 15 do
      begin
        all_element_contributions[j][i] := 0;

      end;
    end;

    // choose element to use as degree of freedom

    degree_of_freedom := hb_freedom.Form12.ComboBox1.Items[hb_freedom.Form12.ComboBox1.ItemIndex] ;

     // assign all element names to all_element_names array

    for i := 1 to 16 do

    begin

    all_element_names[i - 1] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption   ;

    if (all_element_names[i - 1] <> 'Si') then
    all_element_targets[i - 1] := (conc_factor[i-1])*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text ) ;

    if (all_element_names[i - 1] = 'Si') and (ComboBox3.ItemIndex = 0) then
    all_element_targets[i - 1] := (conc_factor[i-1])*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text ) ;

    if (all_element_names[i - 1] = 'Si') and (ComboBox3.ItemIndex = 1) then
    all_element_targets[i - 1] := 0.4684*(conc_factor[i-1])*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text ) ;
   ;

    end ;

    // load the database in order to get the weights and find the resulting ppm values

    MyDbf := TDbf.Create(nil);
    MyDbf.FilePathFull := '';
    MyDbf.TableName := substances_used_db ;
    MyDbf.Open;
    MyDbf.Active := True;

    MyDbf.First; // moves to the first data

    i := 0 ;

    while not MyDbf.EOF do
    begin

      IsLiquid[0][i] := MyDbf.FieldByName('IsLiquid').AsFloat ;
      IsLiquid[1][i] := MyDbf.FieldByName('Density').AsFloat ;

      for j := 0 to 15 do

      begin
        if IsLiquid[0][i] = 0 then preloaded_weight[i] :=  (1 / weight_factor) * MyDbf.FieldByName('Weight').AsFloat ;
        if IsLiquid[0][i] = 1 then preloaded_weight[i] :=   MyDbf.FieldByName('Weight').AsFloat ;
        if IsLiquid[0][i] = 1 then all_solids := False;

        all_element_contributions[j][i] :=
        0.01 * MyDbf.FieldByName(all_element_names[j]).AsFloat * MyDbf.FieldByName('Purity').AsFloat / Volume;

       // if all_element_contributions[j][i] <> 0 then
       // ShowMessage(all_element_names[j] + ' ' +FloatToStr(all_element_contributions[j][i]));

      end;



       if (preloaded_weight[i] > 0) and (CheckBox3.Checked = false) then
       begin

            ShowMessage('You have selected a fixed weight for ' + MyDbf.FieldByName('Name').AsString + '. The program will find a solution with this constraint but it will add more of this substance if it leads to a better solution');

            if RadioButton6.Checked then
            preloaded_weight[i] := preloaded_weight[i]/StrToInt(Edit17.Text);

       end;

       i := i + 1;
      MyDbf.Next;// use .next here NOT .findnext!
    end;

    MyDbf.Close;
    MyDbf.Free;

    // define variable count, names and intended values

    j := 0;

    varnames := nil ;
    vartargetvalue  := nil ;
    temp := nil ;
    temp2  := nil ;


    for i := 1 to 16 do
    begin

      if (StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text) > 0) and
        ((FindComponent('Label' + IntToStr(i)) as TLabel).Caption <> degree_of_freedom) then
      begin

        varcount := varcount + 1;

        if j > 0 then
        begin

          for k := 0 to j - 1 do

          begin

            temp[k]  := varnames[k];
            temp2[k] := vartargetvalue[k];

          end;
        end;


        //increase size of arrays to fit new data
        SetLength(varnames, j + 1);
        SetLength(vartargetvalue, j + 1);



        if j > 0 then

        begin

          for k := 0 to j - 1 do

          begin

            varnames[k] := temp[k];
            vartargetvalue[k] := temp2[k];

          end;
        end;

        temp3 := 0 ;

        for n := 0 to arraysize - 1 do

            begin

            temp3 :=  preloaded_weight[n]*all_element_contributions[i- 1][n] + temp3 ;

             end ;

        varnames[j] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption;

        if (varnames[j] <> 'Si') or (ComboBox3.ItemIndex = 0) then
        vartargetvalue[j] := conc_factor[i-1]*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text) - waterquality[i - 1] - temp3;

        if (varnames[j] = 'Si') and (ComboBox3.ItemIndex = 1) then
        vartargetvalue[j] := 0.4684*conc_factor[i-1]*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text) - waterquality[i - 1] - temp3;

        SetLength(temp, j + 1);
        SetLength(temp2, j + 1);

        j := j + 1;

      end;

    end;



    // assign element target values

    for i := 1 to 16 do

    begin

    all_element_targets[i - 1] := conc_factor[i-1]*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text )  ;

    if (i=13) and (ComboBox3.ItemIndex = 1) then
    all_element_targets[i - 1] := 0.4684*conc_factor[i-1]*StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text )  ;


    end ;

    // definition of arrays which require varcount

       SetLength(problem_matrix_left, varcount, arraysize);

       SetLength(problem_matrix_right, varcount);


    // set array to 0

    for i := 0 to arraysize - 1 do
    begin


      for j := 0 to varcount - 1 do
      begin
        problem_matrix_left[j][i] := 0;

      end;
    end;

    // first get names from substances_used database

    i := 0;

    MyDbf := TDbf.Create(nil);
    MyDbf.FilePathFull := '';
    MyDbf.TableName := substances_used_db ;
    MyDbf.Open;
    MyDbf.Active := True;

    MyDbf.First;                  // moves to the first data

    while not MyDbf.EOF do
    begin
      name_array[i][0] := MyDbf.FieldByName('Name').AsString;
      name_array[i][1] := MyDbf.FieldByName('Formula').AsString;
      name_array[i][2] := MyDbf.FieldByName('Source').AsString;

      // if conditional for when A+B solutions are needed
      if Radiobutton6.Checked then

      begin

        ConcTypeArray[i] := MyDbf.FieldByName('ConcType').AsString;

        if ConcTypeArray[i] = '0' then

        begin

          ShowMessage('Substance ' + name_array[i][0] +
            ' is incompatible with concentrated solutions');
          Exit;

        end;

      end;

      for j := 0 to varcount - 1 do
      begin

        problem_matrix_left[j][i] :=
          0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName('Purity').AsFloat / Volume;

      end;


      cost[i] := MyDbf.FieldByName('Cost').AsFloat;

      i := i + 1;
      MyDbf.Next;                                     // use .next here NOT .findnext!
    end;

    MyDbf.Close;
    MyDbf.Free;

    // right hand sided of the problem matrix (intended concentrations)

    for j := 0 to varcount - 1 do

    begin

      problem_matrix_right[j] := vartargetvalue[j];

    end;



    // check that there is at least one value for every element being evaluated

    for i := 0 to varcount - 1 do

    begin

      test := 0;

      for j := 0 to arraysize - 1 do

      begin

        test := problem_matrix_left[i][j] + test;

      end;

      if test = 0 then

      begin

        ShowMessage(
          'Error: There must be at least ONE substance containing each element for which a non-zero concentration is desired');
        ShowMessage('No substance is providing ' + varnames[i]);
        Exit;

      end;

    end;
    // solve equations using the MatrixSoleLS function

    RMatrixSolveLS(problem_matrix_left, varcount, arraysize, problem_matrix_right,
      0.0, answer,
      report, solutions);

    // calculate resulting ppm from solution


    SetLength(Result, 16);
    SetLength(gross_error, 16);
    SetLength(instrumental_error, 16);

    for i := 0 to 15 do

    begin

      test  := 0;
      upper := 0;
      lower := 0;

      for j := 0 to arraysize - 1 do

      begin

        if (solutions[j] > 0)  then

        begin
          test  := solutions[j] * all_element_contributions[i][j] +  preloaded_weight[j]*all_element_contributions[i][j] + test;
          upper := (preloaded_weight[j] + solutions[j] + weight_error) * all_element_contributions[i][j] * Volume /
            (Volume - volume_error) + upper;
          lower := (preloaded_weight[j] + solutions[j] - weight_error) * all_element_contributions[i][j] * Volume /
            (Volume + volume_error) + lower;
        end;

        if (solutions[j] <= 0) and (preloaded_weight[j] > 0)  then

        begin
          test  :=  preloaded_weight[j]*all_element_contributions[i][j] + test;
          upper := (preloaded_weight[j]  + weight_error) * all_element_contributions[i][j] * Volume /
            (Volume - volume_error) + upper;
          lower := (preloaded_weight[j]  - weight_error) * all_element_contributions[i][j] * Volume /
            (Volume + volume_error) + lower;
        end;

      end;

      Result[i] := test;

      if (RadioButton7.Checked = true) and (all_element_targets[i] <> 0) then
        instrumental_error[i] := ((upper - lower) * (100 / all_element_targets[i])) / 2;

      if (RadioButton6.Checked = true) and (all_element_targets[i] <> 0) then
        instrumental_error[i] := ((upper - lower) * (100 / all_element_targets[i])) /
          (2 * StrtoFloat(Edit17.Text));

      if all_element_targets[i] <> 0 then
      gross_error[i] := ((test+waterquality[i]) * 100 / (all_element_targets[i])) - 100;


      if   all_element_targets[i] = 0 then
      begin

      instrumental_error[i] := 0 ;
      gross_error[i] := 0 ;

      end ;


    end;

    // warn about values which are not used

    for i := 0 to arraysize - 1 do

    begin

     if (solutions[i] <= 0) and (preloaded_weight[i] <= 0) then

           begin

           ShowMessage(name_array[i][0] + ' is not used within the final result. Please remove it to obtain a better solution');

           end;

    end;


    // post results on listbox


    for i := 0 to arraysize - 1 do

    begin

      if (solutions[i] > 0) then

      begin

        StringGrid2.Cells[FORMULA_IDX,i+1] := (name_array[i][1]);

        //determine volume unit for description label

        if RadioButton2.Checked then
          volumeunit := 'gallons';

        if RadioButton3.Checked then
          volumeunit := 'cubic meters';

        if RadioButton1.Checked then
          volumeunit := 'liters';

        // list additions for "direct addition" calculation type

        if RadioButton6.Checked = False then
        begin

          StringGrid2.Cells[NAME_IDX,i+1] := (name_array[i][0]);
          StringGrid2.Cells[SOURCE_IDX,i+1] := (name_array[i][2]);
          if IsLiquid[0][i] = 0 then  StringGrid2.Cells[AMOUNT_IDX,i+1] := (FloatToStr(round2(solutions[i] * weight_factor+ preloaded_weight[i] * weight_factor, 3)) );
          if IsLiquid[0][i] = 1 then  StringGrid2.Cells[AMOUNT_IDX,i+1] := (FloatToStr(round2(solutions[i] + preloaded_weight[i], 3)) );
          StringGrid2.Cells[COST_IDX,i+1] := (FloatToStr(round2(0.001 * solutions[i] * cost[i]+ 0.001*preloaded_weight[i]* cost[i], 1)));

          Label20.Caption := 'Values calculated for the preparation of ' +
            Edit18.Text + ' ' + volumeunit;

        end;



        // List additions for A+B calculations

        if RadioButton6.Checked = True then

        begin

          StringGrid2.Cells[NAME_IDX,i+1] := (ConcTypeArray[i] + ' - ' + name_array[i][0]);
          StringGrid2.Cells[SOURCE_IDX,i+1] := (name_array[i][2]);

          if IsLiquid[0][i] = 0 then StringGrid2.Cells[AMOUNT_IDX,i+1] :=(FloatToStr(round2(solutions[i] * StrtoFloat(Edit17.Text) * weight_factor +preloaded_weight[i]* weight_factor*StrtoFloat(Edit17.Text), 3)));
          if IsLiquid[0][i] = 1 then StringGrid2.Cells[AMOUNT_IDX,i+1] :=(FloatToStr(round2(solutions[i] * StrtoFloat(Edit17.Text) +preloaded_weight[i]*StrtoFloat(Edit17.Text), 3)));

          StringGrid2.Cells[COST_IDX,i+1] :=(FloatToStr(round2(0.001 * solutions[i] * cost[i] * StrtoFloat(Edit17.Text)+0.001*preloaded_weight[i] * cost[i]*StrtoFloat(Edit17.Text), 1)));

          Label20.Caption := 'Values calculated for the preparation of ' +
            Edit18.Text + ' ' + volumeunit + ' of A and ' + Edit18.Text + ' ' +
            volumeunit + ' of B solution. Please use ' +
            FloatToStr(round2(1000 / StrtoFloat(Edit17.Text), 3)) +
            'mL of A and B within every Liter of final solution';

        end;

        end;

        if (solutions[i] <= 0) and (preloaded_weight[i] > 0) then

      begin

        StringGrid2.Cells[FORMULA_IDX,i+1] := (name_array[i][1]);

        //determine volume unit for description label

        if RadioButton2.Checked then
          volumeunit := 'gallons';

        if RadioButton3.Checked then
          volumeunit := 'cubic meters';

        if RadioButton1.Checked then
          volumeunit := 'liters';

        // list additions for "direct addition" calculation type

        if RadioButton6.Checked = False then
        begin

          StringGrid2.Cells[NAME_IDX,i+1] := (name_array[i][0]);
          StringGrid2.Cells[SOURCE_IDX,i+1] := (name_array[i][2]);
          StringGrid2.Cells[AMOUNT_IDX,i+1] := (FloatToStr(round2( preloaded_weight[i], 3)));
          StringGrid2.Cells[COST_IDX,i+1] :=(FloatToStr(round2(0.001*preloaded_weight[i]* cost[i], 1)));

          Label20.Caption := 'Values calculated for the preparation of ' +
            Edit18.Text + ' ' + volumeunit;

        end;

        // List additions for A+B calculations

        if RadioButton6.Checked = True then

        begin

          StringGrid2.Cells[NAME_IDX,i+1] :=(ConcTypeArray[i] + ' - ' + name_array[i][0]);
          StringGrid2.Cells[SOURCE_IDX,i+1] := (name_array[i][2]);

          StringGrid2.Cells[AMOUNT_IDX,i+1] :=(FloatToStr(
            round2(preloaded_weight[i]*StrtoFloat(Edit17.Text), 3)));

          StringGrid2.Cells[COST_IDX,i+1] :=(FloatToStr(
            round2(0.001*preloaded_weight[i] * cost[i]*StrtoFloat(Edit17.Text), 1)));

          Label20.Caption := 'Values calculated for the preparation of ' +
            Edit18.Text + ' ' + volumeunit + ' of A and ' + Edit18.Text + ' ' +
            volumeunit + ' of B solution. Please use ' +
            FloatToStr(round2(1000 / StrtoFloat(Edit17.Text), 3)) +
            'mL of A and B within every Liter of final solution';

        end;

      end;

        //add units to use
        if IsLiquid[0][i] = 0 then StringGrid2.Cells[UNIT_IDX,i+1] := mass_unit;
        if IsLiquid[0][i] = 1 then StringGrid2.Cells[UNIT_IDX,i+1] := 'mL';

    end;


    for i := 0 to 15 do

    begin

      test := 0;

      StringGrid1.Cells[0,i+1] := all_element_names[i];

      test := waterquality[i];;

      if prev_conc = 'ppm' then
      StringGrid1.Cells[1,i+1] := (FloatToStr(round2((1/conc_factor[i])*(Result[i] + test), 3)));

      if prev_conc <> 'ppm' then
      StringGrid1.Cells[1,i+1] := (FloatToStrF((1/conc_factor[i])*(Result[i] + test), ffExponent, 4, 2));

      StringGrid1.Cells[2,i+1] := (FloatToStr(round2(gross_error[i], 1)) + '%');

      StringGrid1.Cells[3,i+1] :=('+/- ' + FloatToStr(round2(instrumental_error[i], 1)) + '%');

      if prev_conc = 'ppm' then
      (FindComponent('RLabel' + IntToStr(i+1)) as TLabel).Caption := FloattoStr(round2((1/conc_factor[i])*Result[i] + test, 3)) ;

      if prev_conc <> 'ppm' then
      (FindComponent('RLabel' + IntToStr(i+1)) as TLabel).Caption := FloattoStrF((1/conc_factor[i])*Result[i] + test, ffExponent, 4, 2) ;

      if (prev_conc = 'ppm') and (i=12) and (ComboBox3.ItemIndex = 1) then
      (FindComponent('RLabel' + IntToStr(i+1)) as TLabel).Caption := FloattoStr(round2((1/0.4684)*(1/conc_factor[i])*Result[i] + test, 3)) ;

      if (prev_conc <> 'ppm') and (i=12) and (ComboBox3.ItemIndex = 1) then
      (FindComponent('RLabel' + IntToStr(i+1)) as TLabel).Caption := FloattoStrF((1/0.4684)*(1/conc_factor[i])*Result[i] + test, ffExponent, 4, 2) ;

    end;





    // check magnitude of instrumental and gross error

       // check for popup disabling
    if CheckBox3.Checked = false then

    begin

    for i := 0 to 15 do

    begin

      if (instrumental_error[i] > 100) then
        ShowMessage('Instrumental error is too high on ' + all_element_names[i] +
          ', prepare more solution volume or a more concentrated solution to reduce this error');

      if (gross_error[i] > 100) then
        ShowMessage('Gross errors on ' + all_element_names[i] +
          ' are too high, probably an adequate solution was NOT found. Please try a different substance combination');

    end;

    end ;

    Button10.Enabled := True;

    // calculation of EC by empirical model
    if RadioButton15.Checked = True then
    begin
        predicted_ec := 0;
        for i := 1 to 16 do
        begin
            predicted_ec := conc_factor[i-1]*StrtoFloat(
              (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption) *
              ec_contribution[i - 1] + predicted_ec;
        end;
        predicted_ec := round2(predicted_ec+0.39661671, 3);
    end;

    // Calculation of EC by LMCv2 model
    if RadioButton14.Checked = True then
    begin

        // calculate ionic strength used for conductivity model
        ionic_strength := 0;
        for i := 1 to 16 do ionic_strength := zi[i-1]*zi[i-1]*(StrtoFloat((FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption) /(1000*molar_mass[i-1])) + ionic_strength;

        predicted_ec := 0;
        for i := 1 to 16 do
        begin
            predicted_ec := conc_factor[i-1]
                            * (StrtoFloat((FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption)/(1000*molar_mass[i-1]))
                            * ec_contribution[i - 1]
                            * exp(-0.7025187*sqrt(ionic_strength)*power(zi[i-1],1.5))
                            + predicted_ec;
        end;
        predicted_ec := round2(predicted_ec, 3);
    end;

    Panel6.Caption := 'EC=' + FloattoStr(predicted_ec) + ' mS/cm';

    // analysis of stock solutions

    if RadioButton6.Checked then

    begin

      //first analyze A solution      23.0  158.5   30.3  284.6  146.2   58.0  28.3

      for i := 0 to 15 do

      begin

        test := 0;

        for j := 0 to arraysize - 1 do

        begin

          if ((solutions[j] > 0) or (preloaded_weight[j] > 0)) and (ConcTypeArray[j] = 'A') then

          begin

            test := preloaded_weight[j]*all_element_contributions[i][j]*Volume + solutions[j] * all_element_contributions[i][j] * Volume + test;

          end;

        end;

        // the 100 is to convert to percentage and the 100,000 to convert cubic meter to mL
        // in order to do the %W/V calculation (g to mL).

        test := test * StrtoFloat(Edit17.Text) * 100 / (Volume * 1000000);

        if (all_element_names[i] <> 'P') and (all_element_names[i] <> 'K') then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[0, i+ 1] := (all_element_names[i]);
          hb_stockanalysis.Form8.StringGrid1.Cells[1, i+ 1] := (FloattoStr(round2(test, 5)));
        end;

        if all_element_names[i] = 'P' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[0, i+ 1] := ('P2O5');
          hb_stockanalysis.Form8.StringGrid1.Cells[1, i+ 1] := (
            FloattoStr(round2(test * 2.2915, 3)));
        end;


        if all_element_names[i] = 'K' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[0, i+ 1] := ('K2O');
          hb_stockanalysis.Form8.StringGrid1.Cells[1, i+ 1] := (
            FloattoStr(round2(test * 1.2047, 3)));
        end;

      end;

      //now analyze B

      for i := 0 to 15 do

      begin

        test := 0;

        for j := 0 to arraysize - 1 do

        begin

          if ((solutions[j] > 0) or (preloaded_weight[j] > 0)) and (ConcTypeArray[j] = 'B') then

          begin

            test :=  preloaded_weight[j]*all_element_contributions[i][j]*Volume + solutions[j] * all_element_contributions[i][j] * Volume + test;

          end;

        end;

        // the 100 is to convert to percentage and the 100,000 to convert cubic meter to mL
        // in order to do the %W/V calculation (g to mL).

        test := test * StrtoFloat(Edit17.Text) * 100 / (Volume * 1000000);

        if (all_element_names[i] <> 'P') and (all_element_names[i] <> 'K') then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[2, i+ 1] := (FloattoStr(round2(test, 5)));
        end;

        if all_element_names[i] = 'P' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[2, i+ 1] :=(
            FloattoStr(round2(test * 2.2915, 3)));
        end;


        if all_element_names[i] = 'K' then

        begin
          hb_stockanalysis.Form8.StringGrid1.Cells[2, i+ 1] :=(
            FloattoStr(round2(test * 1.2047, 3)));
        end;

      end;

    end;

    // carry out detailed analysis for each substance

    for i := 0 to arraysize - 1 do

    begin

      for  j := 0 to 15 do

      begin

        if (all_element_contributions[j][i] > 0) and (solutions[i] > 0) then

        begin

          hb_persubstance.Form9.StringGrid1.RowCount := hb_persubstance.Form9.StringGrid1.RowCount + 1 ;

          hb_persubstance.Form9.StringGrid1.Cells[0, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(name_array[i][0]);
          hb_persubstance.Form9.StringGrid1.Cells[1, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(all_element_names[j]);

          temp3 := solutions[i] * all_element_contributions[j][i]+preloaded_weight[i]*all_element_contributions[j][j] ;

           if prev_conc = 'ppm' then
          hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(
            FloattoStr(round2((1/conc_factor[j])*temp3, 3)));

            if prev_conc <> 'ppm' then
            hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(
            FloattoStrF((1/conc_factor[j])*temp3, ffExponent, 4, 2));

        end;

        if (all_element_contributions[j][i] > 0) and (solutions[i] <= 0) and (preloaded_weight[i] > 0) then

        begin

          hb_persubstance.Form9.StringGrid1.RowCount := hb_persubstance.Form9.StringGrid1.RowCount + 1 ;

          hb_persubstance.Form9.StringGrid1.Cells[0, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(name_array[i][0]);
          hb_persubstance.Form9.StringGrid1.Cells[1, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(all_element_names[j]);

          temp3 := preloaded_weight[i]*all_element_contributions[j][i] ;

          if prev_conc = 'ppm' then
          hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(
            FloattoStr(round2((1/conc_factor[j])*temp3, 3)));

            if prev_conc <> 'ppm' then
            hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount - 1] :=(
            FloattoStrF((1/conc_factor[j])*temp3, ffExponent, 4, 2));

        end;

      end;

    end;


    // check for basic compatibility problems of stock solutions

    if  RadioButton6.Checked then

    begin

        if (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 6]) > 0) and
           (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 5]) > 0)  then
           ShowMessage('Your stock solutions have not been properly designed. Currently there is calcium and sulfate within the same solution (a big problem). Please do NOT carry out this preparation');

        if (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[2, 6]) > 0) and
           (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[2, 5]) > 0)  then
           ShowMessage('Your stock solutions have not been properly designed. Currently there is calcium and sulfate within the same solution (a big problem). Please do NOT carry out this preparation');

         if (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 7]) > 0) and
           (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 6]) > 0)  then
           ShowMessage('Your stock solutions have not been properly designed. To avoid precipitation problems please keep sulfate and iron appart');

          if (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 6]) > 0) and
           (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 7]) > 0)  then
           ShowMessage('Your stock solutions have not been properly designed. To avoid precipitation problems please keep sulfate and iron appart');

          if (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 7]) > 0) and
           (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 3]) > 0)  then
           ShowMessage('Your stock solutions have not been properly designed. To avoid precipitation problems please keep phosphate and iron appart');

          if (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 7]) > 0) and
           (StrtoFloat(hb_stockanalysis.Form8.StringGrid1.Cells[1, 3]) > 0)  then
           ShowMessage('Your stock solutions have not been properly designed. To avoid precipitation problems please keep phosphate and iron appart');

          temp1 :=  0 ;

          for i := 1 to hb_stockanalysis.Form8.StringGrid1.RowCount - 1 do

          begin

          temp1 := StrtoFloat( hb_stockanalysis.Form8.StringGrid1.Cells[1, i]) + temp1 ;

          end ;

          if temp1 > 30 then
          ShowMessage('Your A solution might be too concentrated, please make sure you run solubility tests before preparing any large batches') ;


          temp1 :=  0 ;

          for i := 1 to hb_stockanalysis.Form8.StringGrid1.RowCount - 1 do

          begin

          temp1 := StrtoFloat( hb_stockanalysis.Form8.StringGrid1.Cells[1, i]) + temp1 ;

          end ;

          if temp1 > 30 then
          ShowMessage('Your B solution might be too concentrated, please make sure you run solubility tests before preparing any large batches') ;



    end ;



    // END OF IMPORTANT IF STATEMENT
  end;
  // END OF IMPORTANT IF STATEMENT

      {

      From here we will be calulating concentrations from the custom masses input for
      each of the salts added to the substances_used database (the ones on the ListBox2 of
      the salt selection form).

      }


  // SECOND IMPORTANT IF STATEMENT
  if RadioButton5.Checked then
    // SECOND IMPORTANT IF STATEMENT

  begin

    varnames := nil ;
    Result := nil ;
    name_array := nil ;

    SetLength(varnames, 16);
    SetLength(Result, 16);
    SetLength(name_array, arraysize, 2);
    SetLength(IsLiquid, 2, arraysize) ;
    SetLength(all_element_contributions, 16, arraysize);

    for j := 1 to 16 do
    begin

      // load all element names (this time we don't need to discriminate as
      // we simply calculate for everyone
      varnames[j - 1] := (FindComponent('Label' + IntToStr(j)) as TLabel).Caption;

    end;

    // initialize the counter
    i := 0;

    // load the database in order to get the weights and find the resulting ppm values

    MyDbf := TDbf.Create(nil);
    MyDbf.FilePathFull := '';
    MyDbf.TableName := substances_used_db ;
    MyDbf.Open;
    MyDbf.Active := True;

    MyDbf.First;                  // moves to the first data

    while not MyDbf.EOF do
    begin
      name_array[i][0] := MyDbf.FieldByName('Name').AsString;
      name_array[i][1] := MyDbf.FieldByName('Formula').AsString;

      IsLiquid[0][i] := MyDbf.FieldByName('IsLiquid').AsFloat ;
      IsLiquid[1][i] := MyDbf.FieldByName('Density').AsFloat ;

      for j := 0 to 15 do
      begin
        if IsLiquid[0][i] = 1 then all_solids := False;
        all_element_contributions[j][i] :=
        0.01 * MyDbf.FieldByName((FindComponent('Label' + IntToStr(j+1)) as TLabel).Caption).AsFloat * MyDbf.FieldByName('Purity').AsFloat / Volume;
      end;

      if IsLiquid[0][i] = 0 then StringGrid2.Cells[UNIT_IDX,i+1] := mass_unit;
      if IsLiquid[0][i] = 1 then StringGrid2.Cells[UNIT_IDX,i+1] := 'mL';

      for j := 0 to 15 do

      begin

        // ppm values are very easily calculated using all the information
        // within the DB
        Result[j] := (1 / weight_factor) * MyDbf.FieldByName('Weight').AsFloat *
          0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
          'Purity').AsFloat / Volume + Result[j];

        if ((1 / weight_factor) * MyDbf.FieldByName('Weight').AsFloat *
          0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
          'Purity').AsFloat / Volume > 0) then

        begin

          hb_persubstance.Form9.StringGrid1.RowCount := hb_persubstance.Form9.StringGrid1.RowCount + 1 ;

          hb_persubstance.Form9.StringGrid1.Cells[0, hb_persubstance.Form9.StringGrid1.RowCount-1] := (name_array[i][1]);
          hb_persubstance.Form9.StringGrid1.Cells[1, hb_persubstance.Form9.StringGrid1.RowCount-1] :=(varnames[j]);

          if prev_conc = 'ppm' then

          hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1] :=(
            FloattoStr(round2((1 / weight_factor) * MyDbf.FieldByName('Weight').AsFloat *
            0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
            'Purity').AsFloat / Volume, 3)));

          if prev_conc <> 'ppm' then

          hb_persubstance.Form9.StringGrid1.Cells[2, hb_persubstance.Form9.StringGrid1.RowCount-1] :=(
            FloattoStrF((1/conc_factor[j])*(1 / weight_factor) * MyDbf.FieldByName('Weight').AsFloat *
            0.01 * MyDbf.FieldByName(varnames[j]).AsFloat * MyDbf.FieldByName(
            'Purity').AsFloat / Volume, ffExponent, 4, 2));



        end;

      end;

      StringGrid2.RowCount := i+2 ;

      StringGrid2.Cells[NAME_IDX,i+1] := (name_array[i][0]);
      StringGrid2.Cells[FORMULA_IDX,i+1] := (name_array[i][1]);
      StringGrid2.Cells[COST_IDX,i+1] := (FloattoStr(
        round2(MyDbf.FieldByName('Weight').AsFloat * MyDbf.FieldByName('Cost').AsFloat * 0.001 *
        (1 / weight_factor), 1)));
      StringGrid2.Cells[AMOUNT_IDX,i+1] := (MyDbf.FieldByName('Weight').AsString);


      i := i + 1;
      MyDbf.Next;                                     // use .next here NOT .findnext!
    end;

    MyDbf.Close;
    MyDbf.Free;


    for i := 1 to 16 do
    begin

      for j := 0 to 15 do

      begin

        if (FindComponent('Label' + IntToStr(i)) as TLabel).Caption = varnames[j] then

        begin


          StringGrid1.Cells[0,i] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption;

          if prev_conc = 'ppm' then

          begin


          StringGrid1.Cells[1,i] := (FloatToStr(round2(Result[j] + waterquality[i - 1], 3)));

          (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption :=
            FloatToStr(round2(Result[j] + waterquality[i - 1], 3));

          end ;

          if prev_conc <> 'ppm' then

          begin

          StringGrid1.Cells[1,i] :=(FloatToStrF((1/conc_factor[j])*(Result[j] + waterquality[i - 1]), ffExponent, 4, 2));

          (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption :=
            FloatToStrF((1/conc_factor[j])*(Result[j] + waterquality[i - 1]), ffExponent,4, 2);

          end;

        end;

      end;

    end;



    // calculation of EC by empirical model
    if RadioButton15.Checked = True then
    begin
        predicted_ec := 0;
        for i := 1 to 16 do
        begin
            predicted_ec := conc_factor[i-1]*StrtoFloat(
              (FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption) *
              ec_contribution[i - 1] + predicted_ec;
        end;
        predicted_ec := round2(predicted_ec+0.39661671, 3);
    end;

    // Calculation of EC by LMCv2 model
    if RadioButton14.Checked = True then
    begin

        // calculate ionic strength used for conductivity model
        ionic_strength := 0;
        for i := 1 to 16 do ionic_strength := zi[i-1]*zi[i-1]*(StrtoFloat((FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption) /(1000*molar_mass[i-1])) + ionic_strength;

        predicted_ec := 0;
        for i := 1 to 16 do
        begin
            predicted_ec := conc_factor[i-1]
                            * (StrtoFloat((FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption)/(1000*molar_mass[i-1]))
                            * ec_contribution[i - 1]
                            * exp(-0.7025187*sqrt(ionic_strength)*power(zi[i-1],1.5))
                            + predicted_ec;
        end;
        predicted_ec := round2(predicted_ec, 3);
    end;

    Panel6.Caption := 'EC=' + FloattoStr(predicted_ec) + ' mS/cm';

    Button10.Enabled := True;

    //determine volume unit for description label

        if RadioButton2.Checked then
          volumeunit := 'gallons';

        if RadioButton3.Checked then
          volumeunit := 'cubic meters';

        if RadioButton1.Checked then
          volumeunit := 'liters';

    Label20.Caption := 'Values calculated for the preparation of ' +
            Edit18.Text + ' ' + volumeunit;


    // END OF SECOND IMPORTANT IF STATEMENT
  end;
  // END OF SECOND IMPORTANT IF STATEMENT




  //check and assign any empty elements in StringGrid

  for i := 0 to StringGrid2.RowCount - 2 do

  begin

    if (StringGrid2.Cells[NAME_IDX,i+1]) = '' then
    begin
    StringGrid2.Cells[NAME_IDX,i+1] := name_array[i][0] ;
    StringGrid2.Cells[FORMULA_IDX,i+1] := name_array[i][1] ;
    StringGrid2.Cells[AMOUNT_IDX,i+1] := '0' ;
    StringGrid2.Cells[COST_IDX,i+1] := '0' ;
    end;

  end;

  // total cost and mix calculation

   test := 0;

   for i := 1 to StringGrid2.RowCount - 1 do
   begin
      test := StrtoFloat(StringGrid2.Cells[COST_IDX,i]) + test;
      totalWeight := totalWeight + StrtoFloat(StringGrid2.Cells[AMOUNT_IDX,i]);
      for j:= 0 to 15 do mixContribution[j] := mixContribution[j] + StrtoFloat(StringGrid2.Cells[AMOUNT_IDX,i])*all_element_contributions[j][i-1]*Volume;
   end;

   for j := 1 to 16 do
   begin
        hb_analysis.Form11.StringGrid1.Cells[1,j] := FloatToStr(round2(100*mixContribution[j-1]/totalWeight,3));
        if hb_analysis.Form11.StringGrid1.Cells[0,j] = 'K2O' then hb_analysis.Form11.StringGrid1.Cells[1,j] := FloatToStr(round2(1.2047*100*mixContribution[j-1]/totalWeight,3));
        if hb_analysis.Form11.StringGrid1.Cells[0,j] = 'P2O5' then hb_analysis.Form11.StringGrid1.Cells[1,j] := FloatToStr(round2(2.290*100*mixContribution[j-1]/totalWeight,3));
   end;

   if all_solids then Button19.Enabled := True;

  Label18.Caption := ('Total Cost is ' + FloattoStr(round2(test, 1)));


   // post ratios based on results posted on listboxes above

    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N: P: K') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N', 'P', 'K', 3)) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N: P2O5: K2O') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N', 'P2O5', 'K2O', 3)) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N: K') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N', 'K', 'K', 2) ) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('N (NO3-): N (NH4+)') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('N (NO3-)', 'N (NH4+)', 'K', 2) ) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('Ca: Mg') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('Ca', 'Mg', 'K', 2) ) ;
    hb_ratios.Form14.StringGrid1.Cells[0, hb_ratios.Form14.StringGrid1.RowCount - 1] :=('K: Ca') ;
    hb_ratios.Form14.StringGrid1.Cells[1, hb_ratios.Form14.StringGrid1.RowCount - 2] :=(getratio('K', 'Ca', 'Ca', 2) ) ;

   // enable or disable stock solution analysis button
  if RadioButton6.Checked then Button12.Enabled := True else  Button12.Enabled := False;

  // set water quality values
  for j := 1 to 16 do StringGrid1.Cells[4,j] := FloatToStr(waterquality[j-1]);

  if CheckBox3.Checked = false then
  ShowMessage('Calculation carried out successfully :o)');

end;

procedure TForm1.Button3Click(Sender: TObject);
var

  MyDbf: TDbf;
  oldvalues1: array of array of string;
  oldvalues: array of array of double;
  i: integer;
  j: integer;

begin

  MyDbf := TDbf.Create(nil);
  try
    { use relative path to "data" directory }
    MyDbf.FilePath   := '';
    { we want to use Visual dBase VII compatible tables }
    MyDbf.TableLevel := 7;
    MyDbf.Exclusive  := True;
    MyDbf.TableName  := substances_db ;
    with MyDbf.FieldDefs do
    begin
      Add('Name', ftString, 80, False);
      Add('Formula', ftString, 80, False);
      Add('Source', ftString, 80, False);
      Add('Purity', ftFloat, 0, False);
      Add('N (NO3-)', ftFloat, 0, False);
      Add('N (NH4+)', ftFloat, 0, False);
      Add('P', ftFloat, 0, False);
      Add('K', ftFloat, 0, False);
      Add('Mg', ftFloat, 0, False);
      Add('Ca', ftFloat, 0, False);
      Add('S', ftFloat, 0, False);
      Add('B', ftFloat, 0, False);
      Add('Fe', ftFloat, 0, False);
      Add('Zn', ftFloat, 0, False);
      Add('Mn', ftFloat, 0, False);
      Add('Cu', ftFloat, 0, False);
      Add('Mo', ftFloat, 0, False);
      Add('Na', ftFloat, 0, False);
      Add('Si', ftFloat, 0, False);
      Add('Cl', ftFloat, 0, False);
      Add('isLiquid', ftFloat, 0, False);
      Add('Density', ftFloat, 0, False);
      Add('Cost', ftFloat, 0, False);
      Add('ConcType', ftString, 80, False);
    end;
    MyDbf.CreateTable;
    MyDbf.FieldDefs.Clear;
    MyDbf.Open;
    MyDbf.AddIndex('name', 'Name', [ixCaseInsensitive]);
    MyDbf.AddIndex('formula', 'Formula', [ixCaseInsensitive]);
    MyDbf.Close;
  finally
    MyDbf.Free;
  end;

  SetLength(oldvalues1, 50, 3);
  SetLength(oldvalues, 50, 17);

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := 'substances_win_backup.dbf';
  MyDbf.Open;

  MyDbf.First;

  i := 0;

  while not MyDbf.EOF do
  begin

    oldvalues1[i][0] := MyDbf.FieldByName('Name').AsString;
    oldvalues1[i][1] := MyDbf.FieldByName('Formula').AsString;
    oldvalues1[i][2] := MyDbf.FieldByName('ConcType').AsString;
    oldvalues[i][0]  := MyDbf.FieldByName('Purity').AsFloat;
    ;
    oldvalues[i][1]  := MyDbf.FieldByName('N (NO3-)').AsFloat;
    oldvalues[i][2]  := MyDbf.FieldByName('P').AsFloat;
    oldvalues[i][3]  := MyDbf.FieldByName('K').AsFloat;
    oldvalues[i][4]  := MyDbf.FieldByName('Mg').AsFloat;
    oldvalues[i][5]  := MyDbf.FieldByName('Ca').AsFloat;
    oldvalues[i][6]  := MyDbf.FieldByName('S').AsFloat;
    oldvalues[i][7]  := MyDbf.FieldByName('B').AsFloat;
    oldvalues[i][8]  := MyDbf.FieldByName('Fe').AsFloat;
    oldvalues[i][9]  := MyDbf.FieldByName('Zn').AsFloat;
    oldvalues[i][10] := MyDbf.FieldByName('Cu').AsFloat;
    oldvalues[i][11] := MyDbf.FieldByName('Mo').AsFloat;
    oldvalues[i][12] := MyDbf.FieldByName('Na').AsFloat;
    oldvalues[i][13] := MyDbf.FieldByName('Si').AsFloat;
    oldvalues[i][14] := MyDbf.FieldByName('Cl').AsFloat;
    oldvalues[i][15] := MyDbf.FieldByName('Mn').AsFloat;
    oldvalues[i][16] := MyDbf.FieldByName('N (NH4+)').AsFloat;


    i := i + 1;

    MyDbf.Next;                                     // use .next here NOT .findnext!

  end;

  MyDbf.Close;
  MyDbf.Free;


  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := substances_db;
  MyDbf.Open;

  MyDbf.Active := True;

  for j := 0 to i do

  begin

    MyDbf.Insert;

    MyDbf.FieldByName('Name').AsString := oldvalues1[j][0];
    MyDbf.FieldByName('Formula').AsString := oldvalues1[j][1];
    MyDbf.FieldByName('ConcType').AsString := oldvalues1[j][2];
    MyDbf.FieldByName('Purity').AsFloat := oldvalues[j][0];
    MyDbf.FieldByName('N (NO3-)').AsFloat := oldvalues[j][1];
    MyDbf.FieldByName('P').AsFloat    := oldvalues[j][2];
    MyDbf.FieldByName('K').AsFloat    := oldvalues[j][3];
    MyDbf.FieldByName('Mg').AsFloat   := oldvalues[j][4];
    MyDbf.FieldByName('Ca').AsFloat   := oldvalues[j][5];
    MyDbf.FieldByName('S').AsFloat    := oldvalues[j][6];
    MyDbf.FieldByName('B').AsFloat    := oldvalues[j][7];
    MyDbf.FieldByName('Fe').AsFloat   := oldvalues[j][8];
    MyDbf.FieldByName('Zn').AsFloat   := oldvalues[j][9];
    MyDbf.FieldByName('Cu').AsFloat   := oldvalues[j][10];
    MyDbf.FieldByName('Mo').AsFloat   := oldvalues[j][11];
    MyDbf.FieldByName('Na').AsFloat   := oldvalues[j][12];
    MyDbf.FieldByName('Si').AsFloat   := oldvalues[j][13];
    MyDbf.FieldByName('Cl').AsFloat   := oldvalues[j][14];
    MyDbf.FieldByName('Mn').AsFloat   := oldvalues[j][15];
    MyDbf.FieldByName('N (NH4+)').AsFloat := oldvalues[j][16];
    MyDbf.FieldByName('Cost').AsFloat := 0;
    ;

    MyDbf.Post;

  end;

  MyDbf.Close;
  MyDbf.Free;

end;

procedure TForm1.Button4Click(Sender: TObject);
var
  MyDbf: TDbf;
begin

  if Edit19.Text = 'Input Formulation Name Here' then
  begin
    ShowMessage('Please input a name within the name input field') ;
    Exit;
  end;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := formulations_db;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.Insert;

  MyDbf.FieldByName('Name').AsString := Edit19.Text;
  MyDbf.FieldByName('N (NO3-)').AsFloat := StrtoFloat(Edit1.Text);
  MyDbf.FieldByName('N (NH4+)').AsFloat := StrtoFloat(Edit2.Text);
  MyDbf.FieldByName('P').AsFloat     := StrtoFloat(Edit3.Text);
  MyDbf.FieldByName('K').AsFloat     := StrtoFloat(Edit4.Text);
  MyDbf.FieldByName('Mg').AsFloat    := StrtoFloat(Edit5.Text);
  MyDbf.FieldByName('Ca').AsFloat    := StrtoFloat(Edit6.Text);
  MyDbf.FieldByName('S').AsFloat     := StrtoFloat(Edit7.Text);
  MyDbf.FieldByName('Fe').AsFloat    := StrtoFloat(Edit8.Text);
  MyDbf.FieldByName('Mn').AsFloat    := StrtoFloat(Edit9.Text);
  MyDbf.FieldByName('Zn').AsFloat    := StrtoFloat(Edit10.Text);
  MyDbf.FieldByName('B').AsFloat     := StrtoFloat(Edit11.Text);
  MyDbf.FieldByName('Cu').AsFloat    := StrtoFloat(Edit12.Text);
  MyDbf.FieldByName('Si').AsFloat    := StrtoFloat(Edit13.Text);
  MyDbf.FieldByName('Mo').AsFloat    := StrtoFloat(Edit14.Text);
  MyDbf.FieldByName('Na').AsFloat    := StrtoFloat(Edit15.Text);
  MyDbf.FieldByName('Cl').AsFloat    := StrtoFloat(Edit16.Text);


  if RadioButton10.Checked then
  MyDbf.FieldByName('Units').AsString := 'ppm' ;

  if RadioButton11.Checked then
  MyDbf.FieldByName('Units').AsString := 'M' ;

  if RadioButton12.Checked then
  MyDbf.FieldByName('Units').AsString := 'mM' ;

  if RadioButton13.Checked then
  MyDbf.FieldByName('Units').AsString := 'mN' ;

  MyDbf.Post;

  MyDbf.Close;

  MyDbf.Free;

  if Edit19.Text <> 'DEFAULT' then
  ShowMessage('Formulation named ' + Edit19.Text + ' has been saved to the Database');

  Form1.UpdateComboBox;

  Button6.Enabled := True;

end;

procedure TForm1.Button4Resize(Sender: TObject);
begin

end;


procedure TForm1.Button5Click(Sender: TObject);

var
  MyDbf: TDbf;

begin

  MyDbf := TDbf.Create(nil);
  try
    { use relative path to "data" directory }
    MyDbf.FilePath   := '';
    { we want to use Visual dBase VII compatible tables }
    MyDbf.TableLevel := 7;
    MyDbf.Exclusive  := True;
    MyDbf.TableName  := 'substances_used_win.dbf';
    with MyDbf.FieldDefs do
    begin
   {   Add('Name', ftString, 80, False);
      Add('Time', ftDateTime, 0, False);
      Add('pH1', ftFloat, 0, False);
      Add('pH2', ftFloat, 0, False);
      Add('EC1', ftFloat, 0, False);
      Add('EC2', ftFloat, 0, False);
      Add('Default', ftString, 80, False);
      }
      Add('Name', ftString, 80, False);
      Add('Formula', ftString, 80, False);
      Add('Source', ftString, 80, False);
      Add('Purity', ftFloat, 0, False);
      Add('Cost', ftFloat, 0, False);
      Add('Weight', ftFloat, 0, False);
      Add('N (NO3-)', ftFloat, 0, False);
      Add('P', ftFloat, 0, False);
      Add('K', ftFloat, 0, False);
      Add('Mg', ftFloat, 0, False);
      Add('Ca', ftFloat, 0, False);
      Add('S', ftFloat, 0, False);
      Add('B', ftFloat, 0, False);
      Add('Fe', ftFloat, 0, False);
      Add('Zn', ftFloat, 0, False);
      Add('Cu', ftFloat, 0, False);
      Add('Mo', ftFloat, 0, False);
      Add('Mn', ftFloat, 0, False);
      Add('Na', ftFloat, 0, False);
      Add('Si', ftFloat, 0, False);
      Add('Cl', ftFloat, 0, False);
      Add('N (NH4+)', ftFloat, 0, False);
   //  Add('Default', ftInteger, 0, False);
      Add('isLiquid', ftFloat, 0, False);
      Add('Density', ftFloat, 0, False);
      Add('ConcType', ftString, 80, False);
    end;
    MyDbf.CreateTable;
    MyDbf.FieldDefs.Clear;
    MyDbf.Open;
    //MyDbf.AddIndex('name', 'Name', [ixCaseInsensitive]);
    //MyDbf.AddIndex('default', 'Default', [ixCaseInsensitive]);
    // MyDbf.AddIndex('formula', 'Formula', [ixCaseInsensitive]);
    MyDbf.Close;
  finally
    MyDbf.Free;
  end;

end;

procedure TForm1.Button6Click(Sender: TObject);
var
  MyDbf: TDbf;
  i:     integer;
  selected_item: integer;
begin

  if ComboBox1.Text = 'Select formulation from DB' then
  begin
    ShowMessage('Please select a formulation to delete from the drop down menu') ;
    Exit;
  end;

  if ComboBox1.Items.Count = 0 then
    Exit;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := formulations_db;
  MyDbf.Open;
  MyDbf.Active := True;


  MyDbf.Filter := 'Name=' + QuotedStr(ComboBox1.Items[ComboBox1.ItemIndex]);

  MyDbf.Filtered := True;       // This selects the filtered set
  MyDbf.First;                  // moves the the first filtered data
  ComboBox1.Items.Delete(ComboBox1.ItemIndex);
  MyDbf.Delete;

  MyDbf.Close;

  MyDbf.Free;

  if ComboBox1.Items.Count = 0 then

  begin

    ComboBox1.Text  := 'Select formulation from DB';
    Button6.Enabled := False;

  end;

end;

procedure TForm1.Button7Click(Sender: TObject);
var
  MyDbf: TDbf;
  i:     integer;
  j:     integer;
begin

  Form5.ComboBox4.Items.Clear;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := substances_db;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.First;                  // moves to the first data

  while not MyDbf.EOF do
  begin
    Form5.ComboBox4.Items.Add(MyDbf.FieldByName('Name').AsString);
    MyDbf.Next;                                     // use .next here NOT .findnext!
  end;

  MyDbf.Close;

  Form5.ComboBox4.Sorted := true ;
  hb_commercialnutrient.Form5.Visible := True;
end;

procedure TForm1.Button8Click(Sender: TObject);
begin

  hb_waterquality.Form6.Visible := True;

  hb_waterquality.Form6.UpdateComboBox;

  if hb_waterquality.Form6.ComboBox1.Items.Count = 0 then

  begin

    hb_waterquality.Form6.ComboBox1.Text  := 'Select Water Quality Data From DB';
    hb_waterquality.Form6.Button2.Enabled := False;

  end;

end;

procedure TForm1.Button9Click(Sender: TObject);
begin

  hb_insprecision.Form7.Visible := True;

end;

procedure TForm1.CheckBox1Change(Sender: TObject);
begin

  if CheckBox1.Checked then
  begin

    Button3.Enabled := True;
    Button5.Enabled := True;

  end;

  if CheckBox1.Checked = False then
  begin

    Button3.Enabled := False;
    Button5.Enabled := False;

  end;

end;

procedure TForm1.setinivalues ;
var
  update_ini : TStringList ;
  answer : array of string ;
begin

  SetLength(answer, 3) ;

  if CheckBox3.Checked then
  answer[1] := 'POPUP = 1'
  else
  answer[1] := 'POPUP = 0' ;

  if CheckBox5.Checked then
  answer[2] := 'SMALLWINDOW = 1'
  else
  answer[2] := 'SMALLWINDOW = 0' ;

  update_ini :=  TStringList.Create;

  update_ini.Add(answer[0]) ;

  update_ini.Add(answer[1]) ;

  update_ini.Add(answer[2]) ;

  update_ini.SaveToFile('update.ini');

  update_ini.Free ;

end;

procedure TForm1.CheckBox3Change(Sender: TObject);
begin

  setinivalues;


end;


procedure TForm1.CheckBox4Change(Sender: TObject);
begin

  setinivalues ;

end;

procedure TForm1.CheckBox5Change(Sender: TObject);
begin

  setinivalues ;

end;

procedure TForm1.ComboBox1Change(Sender: TObject);
begin

end;


procedure TForm1.ComboBox1Select(Sender: TObject);
var
  i:     integer;
  selected_item: integer;
  MyDbf: TDbf;
  Units : string ;
begin

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := formulations_db;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.Filter := 'Name=' + QuotedStr(ComboBox1.Items[ComboBox1.ItemIndex]);

  MyDbf.Filtered := True;       // This selects the filtered set
  MyDbf.First;                  // moves the the first filtered data

  Units := MyDbf.FieldByName('Units').AsString ;

  if Units = 'ppm' then
  RadioButton10.Checked := true ;

  if Units = 'M' then
  RadioButton11.Checked := true ;

  if Units = 'mM' then
  RadioButton12.Checked := true ;

  if Units = 'mN' then
  RadioButton13.Checked := true ;

  Edit19.Text := MyDbf.FieldByName('Name').AsString;
  Edit1.Text  := MyDbf.FieldByName('N (NO3-)').AsString;
  Edit2.Text := MyDbf.FieldByName('N (NH4+)').AsString;
  Edit3.Text  := MyDbf.FieldByName('P').AsString;
  Edit4.Text  := MyDbf.FieldByName('K').AsString;
  Edit5.Text  := MyDbf.FieldByName('Mg').AsString;
  Edit6.Text  := MyDbf.FieldByName('Ca').AsString;
  Edit7.Text  := MyDbf.FieldByName('S').AsString;
  Edit8.Text  := MyDbf.FieldByName('Fe').AsString;
  Edit9.Text := MyDbf.FieldByName('Mn').AsString;
  Edit10.Text  := MyDbf.FieldByName('Zn').AsString;
  Edit11.Text  := MyDbf.FieldByName('B').AsString;
  Edit12.Text := MyDbf.FieldByName('Cu').AsString;
  Edit13.Text := MyDbf.FieldByName('Si').AsString;
  Edit14.Text := MyDbf.FieldByName('Mo').AsString;
  Edit15.Text := MyDbf.FieldByName('Na').AsString;
  Edit16.Text := MyDbf.FieldByName('Cl').AsString;

  MyDbf.Close;

  MyDbf.Free;

  Button6.Enabled := True;

end;

procedure TForm1.ComboBox3Change(Sender: TObject);
begin

end;

procedure TForm1.FormClose(Sender: TObject; var CloseAction: TCloseAction);
var
    Sett : TIniFile;
    j: integer;
begin
    //save program variable states on exit
    Sett := TIniFile.Create(IniFile);
    for j := 1 to 19 do Sett.WriteString('Main', 'Form1.Edit' + IntToStr(j), (FindComponent('Edit' + IntToStr(j)) as TEdit).Text);
    for j := 1 to 16 do Sett.WriteString('Main', 'Form1.RLabel' + IntToStr(j), (FindComponent('RLabel' + IntToStr(j)) as TLabel).Caption);
    for j := 1 to 34 do Sett.WriteString('Main', 'Form1.Label' + IntToStr(j), (FindComponent('Label' + IntToStr(j)) as TLabel).Caption);
    for j := 1 to 15 do Sett.WriteBool('Main', 'Form1.RadioButton' + IntToStr(j), (FindComponent('RadioButton' + IntToStr(j)) as TRadioButton).Checked);

    Sett.WriteString('Main', 'prev_conc', prev_conc);
    Sett.WriteString('Main', 'Form1.Label20', Label20.Caption);
    Sett.WriteString('Main', 'Form1.Label18', Label18.Caption);
    Sett.WriteString('Main', 'Form1.Panel6', Panel6.Caption);
    Sett.WriteBool('Main', 'Form1.Checkbox3', Checkbox3.Checked);
    Sett.WriteBool('Main', 'Form1.Checkbox5', Checkbox5.Checked);
    Sett.WriteInteger('Main', 'Form1.ComboBox3', ComboBox3.ItemIndex);
    Sett.WriteBool('Main', 'Form2.CheckBox1', hb_load_salts.Form2.Checkbox1.checked);

    hb_comparison.Form15.StringGrid1.SavetoCSVFile('hb_comparison.csv');
    hb_stockanalysis.Form8.StringGrid1.SavetoCSVFile('hb_stockanalysis.csv');
    hb_persubstance.Form9.StringGrid1.SavetoCSVFile('hb_persubstance.csv');
    hb_ratios.Form14.StringGrid1.SavetoCSVFile('hb_ratios.csv');
    StringGrid1.SavetoCSVFile('hb_ppm_results.csv');
    StringGrid2.SavetoCSVFile('hb_results.csv');

    if hb_comparison.Form15.StringGrid1.ColCount = 1 then  DeleteFile('hb_comparison.csv');
    Sett.Free;
end;

procedure TForm1.LoadValues;
var
    Sett : TIniFile;
    j: integer;
begin
    //load program variables
    Sett := TIniFile.Create(IniFile);
    for j := 1 to 19 do (FindComponent('Edit' + IntToStr(j)) as TEdit).Text := Sett.ReadString('Main', 'Form1.Edit' + IntToStr(j), (FindComponent('Edit' + IntToStr(j)) as TEdit).Text);
    for j := 1 to 16 do (FindComponent('RLabel' + IntToStr(j)) as TLabel).Caption := Sett.ReadString('Main', 'Form1.RLabel' + IntToStr(j), '0');
    for j := 1 to 34 do (FindComponent('Label' + IntToStr(j)) as TLabel).Caption := Sett.ReadString('Main', 'Form1.Label' + IntToStr(j), (FindComponent('Label' + IntToStr(j)) as TLabel).Caption);

    prev_conc := Sett.ReadString('Main', 'prev_conc', prev_conc);
    for j := 1 to 15 do (FindComponent('RadioButton' + IntToStr(j)) as TRadioButton).Checked := Sett.ReadBool('Main', 'Form1.RadioButton' + IntToStr(j), (FindComponent('RadioButton' + IntToStr(j)) as TRadioButton).Checked);

    Label20.Caption := Sett.ReadString('Main', 'Form1.Label20', Label20.Caption);
    Label18.Caption := Sett.ReadString('Main', 'Form1.Label18', Label18.Caption);
    Panel6.Caption := Sett.ReadString('Main', 'Form1.Panel6', Panel6.Caption);
    Checkbox3.Checked := Sett.ReadBool('Main', 'Form1.Checkbox3', Checkbox3.Checked);
    Checkbox5.Checked := Sett.ReadBool('Main', 'Form1.Checkbox5', Checkbox3.Checked);
    ComboBox3.ItemIndex := Sett.ReadInteger('Main', 'Form1.ComboBox3', ComboBox3.ItemIndex);
    hb_load_salts.Form2.Checkbox1.checked := Sett.ReadBool('Main', 'Form2.CheckBox1', hb_load_salts.Form2.Checkbox1.checked);

    if FileExists('hb_comparison.csv') then hb_comparison.Form15.StringGrid1.LoadFromCSVFile('hb_comparison.csv');
    if FileExists('hb_stockanalysis.csv') then hb_stockanalysis.Form8.StringGrid1.LoadFromCSVFile('hb_stockanalysis.csv');
    if FileExists('hb_persubstance.csv') then hb_persubstance.Form9.StringGrid1.LoadFromCSVFile('hb_persubstance.csv');
    if FileExists('hb_ratios.csv') then hb_ratios.Form14.StringGrid1.LoadFromCSVFile('hb_ratios.csv');
    if FileExists('hb_ppm_results.csv') then StringGrid1.LoadFromCSVFile('hb_ppm_results.csv');
    if FileExists('hb_results.csv') then StringGrid2.LoadFromCSVFile('hb_results.csv');
    Sett.Free;
end;

procedure TForm1.FormCreate(Sender: TObject);
begin
  StringGrid1.ShowHint:=True;
  StringGrid1.OnShowHint:=@GridShowHint;
end;

procedure TForm1.GridShowHint(Sender: TObject; HintInfo: PHintInfo);
var
  col: integer = -1;
  row: integer = -1;
  grid: TStringGrid absolute Sender;
begin
  grid.MouseToCell(HintInfo^.CursorPos.X, HintInfo^.CursorPos.Y, col, row);
  if col = 1 then HintInfo^.HintStr:='Final ppm of solution';
  if col = 2 then HintInfo^.HintStr:='Gross Error';
  if col = 3 then HintInfo^.HintStr:='Instrumental Error';
  if col = 4 then HintInfo^.HintStr:='Source water ppm contribution';
  HintInfo^.HideTimeout:=5000; // long-lasting hint
end;


procedure TForm1.FormWindowStateChange(Sender: TObject);
begin

  If WindowState = wsMaximized then
  begin
  WindowState := wsNormal ;
  ShowMessage('Maximizing is currently not supported due to limitations of the window resizing Lazarus implementations') ;
  end;

end;

procedure TForm1.MenuItem2Click(Sender: TObject);
begin

  StringGrid2.CopyToClipboard(true);
  StringGrid2.CopyToClipboard ;

end;

procedure TForm1.PageControl1Exit(Sender: TObject);
begin

end;

procedure TForm1.Panel6Click(Sender: TObject);
begin

end;

procedure TForm1.cleanresults;
var
  i: integer ;
begin

StringGrid1.Clean ;
StringGrid2.Clean ;
{
hb_ratios.Form14.StringGrid1.Clean;
hb_persubstance.Form9.StringGrid1.Clean;
hb_stockanalysis.Form8.StringGrid1.Clean;
}

for i := 1 to 16 do
(FindComponent('RLabel' + IntToStr(i)) as TLabel).Caption := '0' ;

end;

procedure TForm1.RadioButton10Change(Sender: TObject);
var
  i: integer ;
  molar_masses : array[0..15] of double ;
  equivalents :  array[0..15] of integer ;
  prev_value : double ;
begin

  Label29.Caption := 'Conc. (ppm)' ;
  Label32.Caption := '(ppm)' ;
  StringGrid1.Cells[1,0] := 'Result (ppm)' ;
  hb_persubstance.Form9.StringGrid1.Cells[2, 0] := 'Contribution (ppm)' ;
  hb_ratios.Form14.StringGrid1.Cells[1, 0] := 'Ratio (ppm: ppm)' ;

  getmolarmasses(molar_masses) ;
  getequivalents(equivalents) ;

  ///////

  if prev_conc = 'mol/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*molar_masses[i-1]*1000, ffGeneral, 4, 2)

  end ;

  end;

  /////

  if prev_conc = 'mmol/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*molar_masses[i-1], ffGeneral, 4, 2)

  end ;

  end;

   ///////

  if prev_conc = 'meq/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*molar_masses[i-1]/equivalents[i-1], ffGeneral, 4, 2)

  end ;


  end;

  prev_conc := 'ppm' ;
  cleanresults ;

end;

procedure TForm1.RadioButton11Change(Sender: TObject);
var
 i: integer ;
  molar_masses : array[0..15] of double ;
  equivalents :  array[0..15] of integer ;
  prev_value : double ;
begin

  Label29.Caption := 'Conc. (mol/L)' ;
  Label32.Caption := '(mol/L)' ;
  StringGrid1.Cells[1,0] := 'Result (mol/L)' ;
  hb_persubstance.Form9.StringGrid1.Cells[2, 0] := 'Contribution (mol/L)' ;
  hb_ratios.Form14.StringGrid1.Cells[1, 0] := 'Ratio (mol/L: mol/L)' ;

  getmolarmasses(molar_masses) ;
  getequivalents(equivalents)  ;

   ///////

  if prev_conc = 'ppm' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF((prev_value/molar_masses[i-1])/1000, ffExponent, 4, 2)

  end ;

  end;

  /////

  if prev_conc = 'mmol/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value/1000, ffExponent, 4, 2)

  end ;

  end;

   ///////

  if prev_conc = 'meq/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(0.001*prev_value/equivalents[i-1], ffExponent, 4, 2)

  end ;


  end;

  prev_conc := 'mol/L' ;
  cleanresults ;

end;

procedure TForm1.RadioButton12Change(Sender: TObject);
var
   i: integer ;
  molar_masses : array[0..15] of double ;
  equivalents :  array[0..15] of integer ;
  prev_value : double ;
begin

  Label29.Caption := 'Conc. (mmol/L)' ;
  Label32.Caption := '(mmol/L)' ;
  StringGrid1.Cells[1,0] := 'Result (mmol/L)' ;
  hb_persubstance.Form9.StringGrid1.Cells[2, 0] := 'Contribution (mmol/L)' ;
  hb_ratios.Form14.StringGrid1.Cells[1, 0]  := 'Ratio (mmol/L: mmol/L)' ;

  getmolarmasses(molar_masses) ;
  getequivalents(equivalents) ;

   ///////

  if prev_conc = 'ppm' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value/molar_masses[i-1], ffExponent, 4, 2)

  end ;

  end;

  /////

  if prev_conc = 'mol/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*1000, ffExponent, 4, 2)

  end ;

  end;

   ///////

  if prev_conc = 'meq/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value/equivalents[i-1], ffExponent, 4, 2)

  end ;


  end;

  prev_conc := 'mmol/L' ;
  cleanresults ;

end;

procedure TForm1.RadioButton13Change(Sender: TObject);
var
   i: integer ;
  molar_masses : array[0..15] of double ;
  equivalents :  array[0..15] of integer ;
  prev_value : double ;
begin

  Label29.Caption := 'Conc. (meq/L)' ;
  Label32.Caption := '(meq/L)' ;
  StringGrid1.Cells[1,0] := 'Result (meq/L)' ;
  hb_persubstance.Form9.StringGrid1.Cells[2, 0] := 'Contribution (meq/L)' ;
  hb_ratios.Form14.StringGrid1.Cells[1, 0] := 'Ratio (meq/L: meq/L)' ;

  getmolarmasses(molar_masses) ;
  getequivalents(equivalents)  ;

   ///////

  if prev_conc = 'ppm' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*equivalents[i-1]/molar_masses[i-1], ffExponent, 4, 2)

  end ;

  end;

  /////

  if prev_conc = 'mmol/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*equivalents[i-1], ffExponent, 4, 2)

  end ;

  end;

   ///////

  if prev_conc = 'mol/L' then

  begin

  for i := 1 to 16 do

  begin

  prev_value := StrtoFloat((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)  ;

  (FindComponent('Edit' + IntToStr(i)) as TEdit).Text  := FloatToStrF(prev_value*equivalents[i-1]*1000, ffExponent, 4, 2)

  end ;


  end;

  prev_conc := 'meq/L' ;
  cleanresults ;

end;

procedure TForm1.RadioButton1Change(Sender: TObject);
begin
  cleanresults;
end;

procedure TForm1.RadioButton2Change(Sender: TObject);
begin
  cleanresults;
end;

procedure TForm1.RadioButton3Change(Sender: TObject);
begin
  cleanresults;
end;

procedure TForm1.RadioButton4Change(Sender: TObject);
var
  i: integer;
begin

  if RadioButton4.Checked then

  begin

    Button10.Caption := 'Copy Weight Results to DB';
    Button10.Enabled   := False;
    Button12.Enabled   := False;

    for i := 1 to 16 do

    begin

      (FindComponent('Edit' + IntToStr(i)) as TEdit).Enabled := True;
      Panel2.Enabled := True;

    end;

  end;

  cleanresults;

end;

procedure TForm1.RadioButton5Change(Sender: TObject);
var
  i: integer;
begin

  if RadioButton5.Checked then

  begin

    Button10.Caption := 'Copy ppm results to formulation';
    Button10.Enabled := False;
    Button12.Enabled := False;

    for i := 1 to 16 do

    begin

      (FindComponent('Edit' + IntToStr(i)) as TEdit).Enabled := False;
      RadioButton7.Checked := True;
      Panel2.Enabled := False;

    end;

  end;

  cleanresults;

end;

procedure TForm1.RadioButton6Change(Sender: TObject);
begin

  Label17.Caption := 'Stock solution volume';
  Edit17.Enabled  := True;
  cleanresults;

end;

procedure TForm1.RadioButton7Change(Sender: TObject);
begin

  Label17.Caption := 'Volume';
  Edit17.Enabled  := False;
  cleanresults;

end;

procedure TForm1.RadioButton8Change(Sender: TObject);
begin

  StringGrid2.Cells[AMOUNT_IDX,0]  := 'Mass (g)';
  hb_addweight.Form4.Label1.Caption := 'Mass of Substance Used (g)';
  cleanresults;

end;

procedure TForm1.RadioButton9Change(Sender: TObject);
begin

  StringGrid2.Cells[AMOUNT_IDX,0]  := 'Mass (oz)';
  hb_addweight.Form4.Label1.Caption := 'Mass of Substance Used (oz)';
  cleanresults;

end;

procedure TForm1.StringGrid2DrawCell(Sender: TObject; aCol, aRow: Integer;
  aRect: TRect; aState: TGridDrawState);
begin

end;

procedure TForm1.StringGrid2EditButtonClick(Sender: TObject);
var
  i, j: integer;
begin

     i := StringGrid2.Row;
     j := StringGrid2.Col;

     if StringGrid2.Cells[SOURCE_IDX,i] <> '' then
     begin
         OpenURL(StringGrid2.Cells[SOURCE_IDX,i]);
     end;

end;

procedure TForm1.StringGrid2EditingDone(Sender: TObject);
begin
  weightFineTunning;
end;

procedure TForm1.TabSheet1ContextPopup(Sender: TObject; MousePos: TPoint;
  var Handled: Boolean);
begin

end;

procedure TForm1.ToggleBox1Change(Sender: TObject);
begin
  hb_datasetname.Form10.Visible := true;
end;





initialization
  {$I hb_main.lrs}

end.

