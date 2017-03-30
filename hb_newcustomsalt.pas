unit hb_newcustomsalt;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Buttons, Dbf, db, Dbf_Common ;

type

  { TForm3 }

  TForm3 = class(TForm)
    Button1: TBitBtn;
    Button2: TBitBtn;
    CheckBox2: TCheckBox;
    ComboBox1: TComboBox;
    ComboBox2: TComboBox;
    ComboBox3: TComboBox;
    Edit1: TEdit;
    Edit10: TEdit;
    Edit11: TEdit;
    Edit12: TEdit;
    Edit13: TEdit;
    Edit14: TEdit;
    Edit15: TEdit;
    Edit16: TEdit;
    Edit17: TEdit;
    Edit18: TEdit;
    Edit19: TEdit;
    Edit2: TEdit;
    Edit20: TEdit;
    Edit21: TEdit;
    Edit22: TEdit;
    Edit3: TEdit;
    Edit4: TEdit;
    Edit5: TEdit;
    Edit6: TEdit;
    Edit7: TEdit;
    Edit8: TEdit;
    Edit9: TEdit;
    Label1: TLabel;
    Label10: TLabel;
    Label11: TLabel;
    Label12: TLabel;
    Label13: TLabel;
    Label14: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    Formula: TLabel;
    Label17: TLabel;
    Label18: TLabel;
    Label19: TLabel;
    Label2: TLabel;
    Label20: TLabel;
    Label21: TLabel;
    Label22: TLabel;
    Label23: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure CheckBox2Change(Sender: TObject);
    procedure ComboBox1Change(Sender: TObject);
    procedure ComboBox2Change(Sender: TObject);
    procedure ComboBox3Change(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
    end;

var
  Form3: TForm3; 

implementation

uses HB_Main ;

{ TForm3 }

procedure TForm3.Button1Click(Sender: TObject);
var
MyDbf: TDbf;
currentValP: integer;
currentValK: Integer;
currentValSi: Integer;
begin

currentValP := ComboBox1.ItemIndex;
currentValK := ComboBox2.ItemIndex;
currentValSi := ComboBox3.ItemIndex;

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := 'substances.dbf';
MyDbf.Open             ;
MyDbf.Active := true ;

MyDbf.Insert ;

MyDbf.FieldByName('Name').AsString:= Edit15.Text ;
MyDbf.FieldByName('Formula').AsString:= Edit17.Text;
MyDbf.FieldByName('Purity').AsFloat:=StrtoFloat(Edit16.Text)/100 ;

if currentValP = 0 then
MyDbf.FieldByName('P').AsFloat:=StrtoFloat(Edit3.Text);

if currentValK = 0 then
MyDbf.FieldByName('K').AsFloat:=StrtoFloat(Edit2.Text);

if currentValSi = 0 then
MyDbf.FieldByName('Si').AsFloat:=StrtoFloat(Edit13.Text);

if currentValP = 1 then
MyDbf.FieldByName('P').AsFloat:=(StrtoFloat(Edit3.Text)*0.4364);

if currentValK = 1 then
MyDbf.FieldByName('K').AsFloat:=(StrtoFloat(Edit2.Text)*0.8301);

if currentValSi = 1 then
MyDbf.FieldByName('Si').AsFloat:=(StrtoFloat(Edit13.Text)*0.4684);


if CheckBox2.Checked = false then
   MyDbf.FieldByName('IsLiquid').AsInteger:=0;


if CheckBox2.Checked  then
   MyDbf.FieldByName('IsLiquid').AsInteger:=1;

MyDbf.FieldByName('N (NO3-)').AsFloat:=StrtoFloat(Edit1.Text);
MyDbf.FieldByName('N (NH4+)').AsFloat:=StrtoFloat(Edit19.Text);
MyDbf.FieldByName('Mg').AsFloat:=StrtoFloat(Edit4.Text);
MyDbf.FieldByName('Ca').AsFloat:=StrtoFloat(Edit5.Text);
MyDbf.FieldByName('S').AsFloat:=StrtoFloat(Edit6.Text);
MyDbf.FieldByName('B').AsFloat:=StrtoFloat(Edit9.Text);
MyDbf.FieldByName('Fe').AsFloat:=StrtoFloat(Edit7.Text);
MyDbf.FieldByName('Zn').AsFloat:=StrtoFloat(Edit8.Text);
MyDbf.FieldByName('Cu').AsFloat:=StrtoFloat(Edit10.Text);
MyDbf.FieldByName('Mo').AsFloat:=StrtoFloat(Edit11.Text);
MyDbf.FieldByName('Na').AsFloat:=StrtoFloat(Edit12.Text);
MyDbf.FieldByName('Cl').AsFloat:=StrtoFloat(Edit14.Text);
MyDbf.FieldByName('Mn').AsFloat:=StrtoFloat(Edit18.Text);
MyDbf.FieldByName('Cost').AsFloat:=StrtoFloat(Edit21.Text);
MyDbf.FieldByName('ConcType').AsString:=Edit20.Text;
MyDbf.FieldByName('Density').AsString:=Edit22.Text;

MyDbf.Post ;

MyDbf.Close ;

MyDbf.Free ;

Form1.UpdateList ;

if currentValP = 1  then
ShowMessage('P will be converted and saved as P%, to see P2O5 again in the future simply select it from the dropbox for automatic conversion');

if currentValK = 1  then
ShowMessage('K will be converted and saved as K%, to see K2O again in the future simply select it from the dropbox for automatic conversion');

if currentValSi = 1  then
ShowMessage('Si willconverted and be saved as Si%, to see SiO2 again in the future simply select it from the dropbox for automatic conversion');

ComboBox1.ItemIndex := 0;

ComboBox2.ItemIndex := 0 ;

ComboBox3.ItemIndex := 0 ;

Form3.Visible := False ;

end;

procedure TForm3.Button2Click(Sender: TObject);
var
MyDbf: TDbf;
currentValP: integer;
currentValK: Integer;
currentValSi: Integer;
begin

currentValP := ComboBox1.ItemIndex;
currentValK := ComboBox2.ItemIndex;
currentValSi := ComboBox3.ItemIndex;

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := 'substances.dbf';
MyDbf.Open             ;
MyDbf.Active := true ;

MyDbf.Filter := 'Name=' + QuotedStr(Label23.Caption) ;

    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;

    MyDbf.Edit;

               MyDbf.FieldByName('Name').AsString:= Edit15.Text ;
               MyDbf.FieldByName('Formula').AsString:= Edit17.Text;
               MyDbf.FieldByName('Purity').AsFloat:=StrtoFloat(Edit16.Text)/100 ;
               MyDbf.FieldByName('N (NO3-)').AsFloat:=StrtoFloat(Edit1.Text);


    if currentValP = 0 then
    MyDbf.FieldByName('P').AsFloat:=StrtoFloat(Edit3.Text);

    if currentValK = 0 then
    MyDbf.FieldByName('K').AsFloat:=StrtoFloat(Edit2.Text);

    if currentValSi = 0 then
    MyDbf.FieldByName('Si').AsFloat:=StrtoFloat(Edit13.Text);

    if currentValP = 1 then
    MyDbf.FieldByName('P').AsFloat:=(StrtoFloat(Edit3.Text)*0.4364);

    if currentValK = 1 then
    MyDbf.FieldByName('K').AsFloat:=(StrtoFloat(Edit2.Text)*0.8301);

    if currentValSi = 1 then
    MyDbf.FieldByName('Si').AsFloat:=(StrtoFloat(Edit13.Text)*0.4684);

    if CheckBox2.Checked = false then
    MyDbf.FieldByName('IsLiquid').AsInteger:=0;


   if CheckBox2.Checked  then
   MyDbf.FieldByName('IsLiquid').AsInteger:=1;

               MyDbf.FieldByName('Mg').AsFloat:=StrtoFloat(Edit4.Text);
               MyDbf.FieldByName('Ca').AsFloat:=StrtoFloat(Edit5.Text);
               MyDbf.FieldByName('S').AsFloat:=StrtoFloat(Edit6.Text);
               MyDbf.FieldByName('B').AsFloat:=StrtoFloat(Edit9.Text);
               MyDbf.FieldByName('Fe').AsFloat:=StrtoFloat(Edit7.Text);
               MyDbf.FieldByName('Zn').AsFloat:=StrtoFloat(Edit8.Text);
               MyDbf.FieldByName('Cu').AsFloat:=StrtoFloat(Edit10.Text);
               MyDbf.FieldByName('Mo').AsFloat:=StrtoFloat(Edit11.Text);
               MyDbf.FieldByName('Na').AsFloat:=StrtoFloat(Edit12.Text);
               MyDbf.FieldByName('Cl').AsFloat:=StrtoFloat(Edit14.Text);
               MyDbf.FieldByName('Mn').AsFloat:=StrtoFloat(Edit18.Text);
               MyDbf.FieldByName('Cost').AsFloat:=StrtoFloat(Edit21.Text);
               MyDbf.FieldByName('N (NH4+)').AsFloat:=StrtoFloat(Edit19.Text);
               MyDbf.FieldByName('ConcType').AsString:=Edit20.Text;
               MyDbf.FieldByName('Density').AsString:=Edit22.Text;

    MyDbf.Post ;

MyDbf.Close ;

MyDbf.Free ;

ComboBox1.ItemIndex := 0;

ComboBox2.ItemIndex := 0 ;

ComboBox3.ItemIndex := 0 ;

HB_Main.Form1.UpdateList ;

Form3.Visible := False ;



end;


procedure TForm3.CheckBox2Change(Sender: TObject);
begin

if Checkbox2.Checked then

   begin

   Edit22.Visible := true ;
   Label22.Visible := true ;

   end;

    if Checkbox2.Checked = false then

   begin

   Edit22.Visible := false ;
   Label22.Visible := false ;

   end;

end;

procedure TForm3.ComboBox1Change(Sender: TObject);
var
currentVal: integer;
begin

currentVal := ComboBox1.ItemIndex;


    if (currentVal = 1) and (Button2.Enabled) then

       begin

       Edit3.Text := FloattoStr(Form1.round2(StrtoFloat(Edit3.Text)*2.2915, 3))   ;

       end ;

   if (currentVal = 0) and (Button2.Enabled) then

      begin

      Edit3.Text := FloattoStr(Form1.round2(StrtoFloat(Edit3.Text)*(1/2.2915), 3))   ;

      end ;

end;

procedure TForm3.ComboBox2Change(Sender: TObject);
var
currentVal: integer;
begin

currentVal := ComboBox2.ItemIndex;

   if (currentVal = 1) and (Button2.Enabled) then

      begin

      Edit2.Text := FloattoStr(Form1.round2(StrtoFloat(Edit2.Text)*1.2047, 3))   ;

      end ;

   if (currentVal = 0) and (Button2.Enabled) then

      begin

      Edit2.Text := FloattoStr(Form1.round2(StrtoFloat(Edit2.Text)*(1/1.2047), 3))   ;

      end ;

end;

procedure TForm3.ComboBox3Change(Sender: TObject);
var
currentVal: integer;
begin

currentVal := ComboBox3.ItemIndex;

   if (currentVal = 1) and (Button2.Enabled) then

      begin

      Edit13.Text := FloattoStr(Form1.round2(StrtoFloat(Edit13.Text)*2.1348, 3))   ;

      end ;

   if (currentVal = 0) and (Button2.Enabled) then

      begin

      Edit13.Text := FloattoStr(Form1.round2(StrtoFloat(Edit13.Text)*(1/2.1348), 3))   ;

      end ;

end;

initialization
  {$I hb_newcustomsalt.lrs}

end.

