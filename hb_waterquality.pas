unit hb_waterquality;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Buttons, Dbf, db, Dbf_Common, hb_ph;

type

  { TForm6 }

  TForm6 = class(TForm)
    Button1: TBitBtn;
    Button2: TBitBtn;
    Button3: TBitBtn;
    Button4: TButton;
    Button5: TButton;
    ComboBox1: TComboBox;
    Edit1: TEdit;
    Edit10: TEdit;
    Edit11: TEdit;
    Edit12: TEdit;
    Edit13: TEdit;
    Edit14: TEdit;
    Edit25: TEdit;
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
    Label26: TLabel;
    Label25: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure ComboBox1Select(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
    procedure UpdateComboBox ;
  end; 

var
  Form6: TForm6; 

implementation

uses HB_Main;

procedure TForm6.Button1Click(Sender: TObject);
 var
MyDbf: TDbf;
begin

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := Form1.water_quality_db;
MyDbf.Open             ;
MyDbf.Active := true ;

MyDbf.Insert ;

MyDbf.FieldByName('Name').AsString:= Edit25.Text ;
MyDbf.FieldByName('P').AsFloat:=StrtoFloat(Edit3.Text);
MyDbf.FieldByName('K').AsFloat:=StrtoFloat(Edit2.Text);
MyDbf.FieldByName('N (NO3-)').AsFloat:=StrtoFloat(Edit1.Text);
MyDbf.FieldByName('N (NH4+)').AsFloat:=StrtoFloat(Edit16.Text);
MyDbf.FieldByName('Mg').AsFloat:=StrtoFloat(Edit4.Text);
MyDbf.FieldByName('Ca').AsFloat:=StrtoFloat(Edit5.Text);
MyDbf.FieldByName('S').AsFloat:=StrtoFloat(Edit6.Text);
MyDbf.FieldByName('B').AsFloat:=StrtoFloat(Edit9.Text);
MyDbf.FieldByName('Fe').AsFloat:=StrtoFloat(Edit7.Text);
MyDbf.FieldByName('Zn').AsFloat:=StrtoFloat(Edit8.Text);
MyDbf.FieldByName('Cu').AsFloat:=StrtoFloat(Edit10.Text);
MyDbf.FieldByName('Mo').AsFloat:=StrtoFloat(Edit11.Text);
MyDbf.FieldByName('Na').AsFloat:=StrtoFloat(Edit12.Text);
MyDbf.FieldByName('Si').AsFloat:=StrtoFloat(Edit13.Text);
MyDbf.FieldByName('Cl').AsFloat:=StrtoFloat(Edit14.Text);
MyDbf.FieldByName('Mn').AsFloat:=StrtoFloat(Edit15.Text);


MyDbf.FieldByName('Default').AsInteger := 0;

MyDbf.Post ;

MyDbf.Close ;

MyDbf.Free ;

ShowMessage('Water Quality data named ' + Edit25.Text + ' has been saved to the Database');

Form6.UpdateComboBox ;

Button3.Enabled := True ;

Button2.Enabled := true ;

end;

procedure TForm6.Button2Click(Sender: TObject);
  var
MyDbf: TDbf;
i : integer ;
selected_item : integer ;
begin

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := Form1.water_quality_db;
MyDbf.Open             ;
MyDbf.Active := true ;


MyDbf.Filter := 'Name=' + QuotedStr(ComboBox1.Items[ComboBox1.ItemIndex]) ;

    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;                  // moves the the first filtered data
    ComboBox1.Items.Delete(ComboBox1.ItemIndex) ;
    MyDbf.Delete ;

MyDbf.Close ;

MyDbf.Free ;

if ComboBox1.Items.Count = 0 then

   begin

   ComboBox1.Text := 'Select Water Quality Data From DB' ;
   Button2.Enabled := false ;

   end;

end;

procedure TForm6.Button3Click(Sender: TObject);
  var
MyDbf: TDbf;
begin

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := Form1.water_quality_db;
MyDbf.Open             ;
MyDbf.Active := true ;

while not MyDbf.EOF do
    begin

        MyDbf.Edit;

        MyDbf.FieldByName('Default').AsInteger := 0 ;

        MyDbf.Post ;

        MyDbf.next;                                   // use .next here NOT .findnext!
    end;

MyDbf.Filter := 'Name=' + QuotedStr(Edit25.Text) ;

    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;

    MyDbf.Edit;

              MyDbf.FieldByName('Default').AsInteger := 1 ;

    MyDbf.Post ;

MyDbf.Close ;

MyDbf.Free ;

ShowMessage(Edit25.Text + ' set as default water quality set') ;



end;

procedure TForm6.Button4Click(Sender: TObject);
begin

Form6.Visible := false ;

end;

procedure TForm6.Button5Click(Sender: TObject);
begin
  hb_ph.Form13.Visible := true ;
end;


procedure TForm6.ComboBox1Select(Sender: TObject);
var
i : integer ;
selected_item : integer ;
MyDbf: TDbf;
begin

   MyDbf := TDbf.Create(nil) ;
   MyDbf.FilePathFull := '';
   MyDbf.TableName := Form1.water_quality_db;
   MyDbf.Open             ;
   MyDbf.Active := true ;

         MyDbf.Filter := 'Name=' + QuotedStr(ComboBox1.Items[ComboBox1.ItemIndex]) ;

    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;                  // moves the the first filtered data

    Edit25.text := MyDbf.FieldByName('Name').AsString;
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


    MyDbf.Close ;

    MyDbf.Free ;

  Button3.Enabled := True ;




end;


procedure TForm6.UpdateComboBox ;
var
MyDbf: TDbf;
i : integer ;
j : integer ;
begin

ComboBox1.Items.Clear ;

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := Form1.water_quality_db;
MyDbf.Open             ;
MyDbf.Active := true ;

    MyDbf.First;                  // moves to the first data

    while not MyDbf.EOF do
    begin
        ComboBox1.Items.Add(MyDbf.FieldByName('Name').AsString);
        MyDbf.next;                                     // use .next here NOT .findnext!
    end;

MyDbf.Close ;

MyDbf.Free ;



end ;

initialization
  {$I hb_waterquality.lrs}

end.

