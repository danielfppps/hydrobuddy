unit hb_commercialnutrient;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  ComCtrls, StdCtrls, Menus, ExtCtrls, Buttons, Dbf, hb_comparison;

type

  { TForm5 }

  TForm5 = class(TForm)
    Button1: TBitBtn;
    Button2: TButton;
    Button3: TButton;
    ComboBox4: TComboBox;
    Edit1: TEdit;
    Edit10: TEdit;
    Edit11: TEdit;
    Edit12: TEdit;
    Edit13: TEdit;
    Edit14: TEdit;
    Edit17: TEdit;
    Edit15: TEdit;
    Edit16: TEdit;
    Edit2: TEdit;
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
    Label17: TLabel;
    Label29: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    Label2: TLabel;
    Label20: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Panel1: TPanel ;
    CheckBox1: TCheckBox;
    RadioButton3: TRadioButton;
    RadioButton4: TRadioButton;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure ComboBox4Change(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
    is_liquid: boolean;
  end; 

var
  Form5: TForm5; 

implementation

uses HB_Main ;

{ TForm5 }


procedure TForm5.Button1Click(Sender: TObject);
var
i : integer ;
j : integer ;
varnames : array of string ;
result : array of double ;
test : double ;
Volume : double ;
begin

SetLength (varnames, 16) ;
SetLength (result, 16) ;

if  RadioButton3.Checked then
Volume := 1      ;

if RadioButton4.Checked then
Volume := 3.78541178 ;

for i := 1 to 16 do
 begin

 // load all element names (this time we don't need to discriminate as
 // we simply calculate for everyone
 varnames[ i - 1 ] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption ;

 end;


    for i := 1 to 16 do

    begin

      test := StrtoFloat(((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)) ;
      result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit17.Text) * 1000 )/Volume, 2) ;

    end;

    HB_Main.Form1.cleanresults ;
    HB_Main.Form1.RadioButton10.Checked := true ;


    // finally copy values to edit boxes



           if CheckBox1.Checked = false then

           begin

          HB_Main.Form1.Edit1.Text := FloattoStr(result[0]) ;
          HB_Main.Form1.Edit3.Text := FloattoStr(result[2]) ;
          HB_Main.Form1.Edit2.Text := FloattoStr(result[1]) ;
          HB_Main.Form1.Edit13.Text := FloattoStr(result[13]) ;
          HB_Main.Form1.Edit4.Text := FloattoStr(result[3]) ;
          HB_Main.Form1.Edit5.Text := FloattoStr(result[4]) ;
          HB_Main.Form1.Edit6.Text := FloattoStr(result[5]) ;
          HB_Main.Form1.Edit7.Text := FloattoStr(result[6]) ;
          HB_Main.Form1.Edit8.Text := FloattoStr(result[7]) ;
          HB_Main.Form1.Edit9.Text := FloattoStr(result[8]) ;
          HB_Main.Form1.Edit10.Text := FloattoStr(result[9]) ;
          HB_Main.Form1.Edit11.Text := FloattoStr(result[10]) ;
          HB_Main.Form1.Edit12.Text := FloattoStr(result[11]) ;
          HB_Main.Form1.Edit13.Text := FloattoStr(result[12]) ;
          HB_Main.Form1.Edit15.Text := FloattoStr(result[14]) ;
          HB_Main.Form1.Edit16.Text := FloattoStr(result[15]) ;

          end ;




          if CheckBox1.Checked then

          begin

          HB_Main.Form1.Edit1.Text := FloattoStr(result[0]+ StrtoFloat(HB_Main.Form1.Edit1.Text)) ;

          HB_Main.Form1.Edit3.Text := FloattoStr(result[2]+ StrtoFloat(HB_Main.Form1.Edit3.Text)) ;

          HB_Main.Form1.Edit2.Text := FloattoStr(result[1]+ StrtoFloat(HB_Main.Form1.Edit2.Text)) ;  ;

          HB_Main.Form1.Edit13.Text := FloattoStr(result[13]+ StrtoFloat(HB_Main.Form1.Edit13.Text)) ;




          HB_Main.Form1.Edit4.Text := FloattoStr(result[3]+ StrtoFloat(HB_Main.Form1.Edit4.Text)) ;
          HB_Main.Form1.Edit5.Text := FloattoStr(result[4]+ StrtoFloat(HB_Main.Form1.Edit5.Text)) ;
          HB_Main.Form1.Edit6.Text := FloattoStr(result[5]+ StrtoFloat(HB_Main.Form1.Edit6.Text)) ;
          HB_Main.Form1.Edit7.Text := FloattoStr(result[6]+ StrtoFloat(HB_Main.Form1.Edit7.Text)) ;
          HB_Main.Form1.Edit8.Text := FloattoStr(result[7]+ StrtoFloat(HB_Main.Form1.Edit8.Text)) ;
          HB_Main.Form1.Edit9.Text := FloattoStr(result[8]+ StrtoFloat(HB_Main.Form1.Edit9.Text)) ;
          HB_Main.Form1.Edit10.Text := FloattoStr(result[9]+ StrtoFloat(HB_Main.Form1.Edit10.Text)) ;
          HB_Main.Form1.Edit11.Text := FloattoStr(result[10]+ StrtoFloat(HB_Main.Form1.Edit11.Text));
          HB_Main.Form1.Edit12.Text := FloattoStr(result[11]+ StrtoFloat(HB_Main.Form1.Edit12.Text));
          HB_Main.Form1.Edit13.Text := FloattoStr(result[12]+ StrtoFloat(HB_Main.Form1.Edit13.Text));
          HB_Main.Form1.Edit14.Text := FloattoStr(result[13]+ StrtoFloat(HB_Main.Form1.Edit14.Text));
          HB_Main.Form1.Edit15.Text := FloattoStr(result[14]+ StrtoFloat(HB_Main.Form1.Edit15.Text));
          HB_Main.Form1.Edit16.Text :=FloattoStr(result[15]+ StrtoFloat(HB_Main.Form1.Edit16.Text));

          end ;


          Form5.Visible := false ;


          end ;

procedure TForm5.Button2Click(Sender: TObject);
var
i : integer ;
j : integer ;
varnames : array of string ;
result : array of double ;
test : double ;
Volume : double ;
colCount: integer;
addition_units: string;
s:TTextStyle;
begin


SetLength (varnames, 16) ;
SetLength (result, 16) ;

if (is_liquid = True) and (RadioButton4.Checked) then addition_units := 'mL/gal' ;
if (is_liquid = False) and (RadioButton4.Checked) then addition_units := 'g/gal' ;
if (is_liquid = True) and (RadioButton3.Checked) then addition_units := 'mL/L'   ;
if (is_liquid = False) and (RadioButton3.Checked) then addition_units := 'g/L'   ;

if  RadioButton3.Checked then
Volume := 1      ;

if RadioButton4.Checked then
Volume := 3.78541178 ;

for i := 1 to 16 do
 begin

 // load all element names (this time we don't need to discriminate as
 // we simply calculate for everyone
 varnames[ i - 1 ] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption ;

 end;


    for i := 1 to 16 do

    begin

    test := StrtoFloat(((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)) ;
    result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit17.Text) * 1000 )/Volume, 2) ;

    end;


    // finally copy values for comparison

    hb_comparison.Form15.StringGrid1.ColCount:= hb_comparison.Form15.StringGrid1.ColCount + 1 ;

    colCount := hb_comparison.Form15.StringGrid1.ColCount ;

    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 0] := ComboBox4.Items[ComboBox4.ItemIndex] ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 1] := FloattoStr(result[0]) ;

    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 2] := FloattoStr(result[2]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 3] := FloattoStr(result[1]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 14] := FloattoStr(result[13]);
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 4] := FloattoStr(result[3]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 5] := FloattoStr(result[4]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 6] := FloattoStr(result[5]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 7] := FloattoStr(result[6]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 8] := FloattoStr(result[7]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 9] := FloattoStr(result[8]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 10] := FloattoStr(result[9]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 11] := FloattoStr(result[10]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 12] := FloattoStr(result[11]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 13] := FloattoStr(result[12]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 15] := FloattoStr(result[14]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 16] := FloattoStr(result[15]) ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 17] := Edit17.Text + addition_units ;

    ShowMessage('Product final ppm values added to comparison chart') ;
    hb_comparison.Form15.StringGrid1.AutoSizeColumn(colCount-1);
    s := hb_comparison.Form15.StringGrid1.DefaultTextStyle;
    s.Alignment:=taCenter;
    hb_comparison.Form15.StringGrid1.DefaultTextStyle := s;
end;

procedure TForm5.Button3Click(Sender: TObject);
begin
  hb_comparison.Form15.Visible := true ;
end;

procedure TForm5.ComboBox4Change(Sender: TObject);
  var
i : integer ;
selected_item : integer ;
MyDbf: TDbf;
begin

   selected_item := ComboBox4.ItemIndex;

   MyDbf := TDbf.Create(nil) ;
   MyDbf.FilePathFull := '';
   MyDbf.TableName := Form1.substances_db;
   MyDbf.Open             ;
   MyDbf.Active := true ;


         MyDbf.Filter := 'Name=' + QuotedStr(ComboBox4.Items[selected_item]) ;

    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;                  // moves the the first filtered data

    Edit1.text := MyDbf.FieldByName('N (NO3-)').AsString ;
    Edit3.text := MyDbf.FieldByName('P').AsString ;
    Edit2.text := MyDbf.FieldByName('K').AsString ;
    Edit4.text := MyDbf.FieldByName('Mg').AsString ;
    Edit5.text := MyDbf.FieldByName('Ca').AsString ;
    Edit6.text := MyDbf.FieldByName('S').AsString ;
    Edit7.text := MyDbf.FieldByName('Fe').AsString ;
    Edit9.text := MyDbf.FieldByName('B').AsString ;
    Edit8.text := MyDbf.FieldByName('Zn').AsString ;
    Edit10.text := MyDbf.FieldByName('Cu').AsString ;
    Edit11.text := MyDbf.FieldByName('Mo').AsString ;
    Edit12.text := MyDbf.FieldByName('Na').AsString ;
    Edit15.text := MyDbf.FieldByName('Mn').AsString ;
    Edit13.text := MyDbf.FieldByName('Si').AsString ;
    Edit14.text := MyDbf.FieldByName('Cl').AsString ;
    Edit16.text := MyDbf.FieldByName('N (NH4+)').AsString ;

    if MyDbf.FieldByName('IsLiquid').AsInteger = 0 then
    begin
         is_liquid := False ;
         Label20.Caption := 'Mass of addition (g)'
    end;

    if MyDbf.FieldByName('IsLiquid').AsInteger = 1 then
    begin
         is_liquid := True ;
         Label20.Caption := 'Volume of addition (mL)'
    end;

    MyDbf.Close ;

    MyDbf.Free ;


end;


initialization
  {$I hb_commercialnutrient.lrs}

end.

