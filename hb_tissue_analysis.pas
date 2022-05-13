unit hb_tissue_analysis;

{$mode objfpc}

interface

uses
  Classes, SysUtils, LResources, Forms, Controls, Graphics, Dialogs, StdCtrls,
  Buttons, Grids, Dbf, db, Dbf_Common;

type

  { TForm16 }

  TForm16 = class(TForm)
    Button1: TBitBtn;
    Button2: TBitBtn;
    Button3: TBitBtn;
    Button4: TBitBtn;
    Button5: TBitBtn;
    Edit1: TEdit;
    Edit10: TEdit;
    Edit11: TEdit;
    Edit12: TEdit;
    Edit13: TEdit;
    Edit14: TEdit;
    Edit15: TEdit;
    Edit16: TEdit;
    Edit17: TEdit;
    Edit2: TEdit;
    Edit25: TEdit;
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
    Label17: TLabel;
    Label26: TLabel;
    Label29: TLabel;
    Label3: TLabel;
    Label30: TLabel;
    Label31: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    ListBox1: TListBox;
    StringGrid1: TStringGrid;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ListBox1Click(Sender: TObject);
    procedure ListBox1SelectionChange(Sender: TObject; User: boolean);
  private

  public
      procedure UpdateTissueList;
  end;

var
  Form16: TForm16;

implementation

uses HB_Main ;

{ TForm16 }

procedure  TForm16.UpdateTissueList;
var
  MyDbf: TDbf;
begin

  ListBox1.Items.Clear;

  MyDbf := TDbf.Create(nil);
  MyDbf.FilePathFull := '';
  MyDbf.TableName := Form1.tissue_analysis_db ;
  MyDbf.Open;
  MyDbf.Active := True;

  MyDbf.First;                  // moves to the first data

  while not MyDbf.EOF do
  begin
    ListBox1.Items.Add(MyDbf.FieldByName('Name').AsString);
    MyDbf.Next;                                     // use .next here NOT .findnext!
  end;

  MyDbf.Close;

  MyDbf.Free;

end;

procedure TForm16.FormCreate(Sender: TObject);
begin

end;

procedure TForm16.ListBox1Click(Sender: TObject);
begin

end;

procedure TForm16.ListBox1SelectionChange(Sender: TObject; User: boolean);
var
  i,selected_idx : integer ;
  item_selected : boolean ;
  MyDbf: TDbf;
  wue: double;
begin

item_selected := false ;
wue := StrToFloat(Edit17.Text)  ;

for i := 0 to ListBox1.Items.Count - 1 do

    begin

    if (ListBox1.Selected [i]) then
    begin
         selected_idx := i;
         item_selected := true ;
    end;

    end ;

if item_selected then

begin

    MyDbf := TDbf.Create(nil) ;
    MyDbf.FilePathFull := '';
    MyDbf.TableName := Form1.tissue_analysis_db;
    MyDbf.Open             ;
    MyDbf.Active := true ;
    MyDbf.Filter := 'Name=' + QuotedStr(ListBox1.Items[selected_idx]) ;
    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;                  // moves the the first filtered data


    Edit25.Text :=  ListBox1.Items[selected_idx];

    Edit1.Text :=  FloattoStr(MyDbf.FieldByName('N').AsFloat);
    Edit3.Text :=  FloattoStr(MyDbf.FieldByName('P').AsFloat);
    Edit4.Text :=  FloattoStr(MyDbf.FieldByName('K').AsFloat);
    Edit5.Text :=  FloattoStr(MyDbf.FieldByName('Mg').AsFloat);
    Edit6.Text :=  FloattoStr(MyDbf.FieldByName('Ca').AsFloat);
    Edit7.Text :=  FloattoStr(MyDbf.FieldByName('S').AsFloat);
    Edit13.Text :=  FloattoStr(MyDbf.FieldByName('Si').AsFloat);
    Edit8.Text :=  FloattoStr(MyDbf.FieldByName('Fe').AsFloat);
    Edit9.Text :=  FloattoStr(MyDbf.FieldByName('Mn').AsFloat);
    Edit10.Text :=  FloattoStr(MyDbf.FieldByName('Zn').AsFloat);
    Edit11.Text :=  FloattoStr(MyDbf.FieldByName('B').AsFloat);
    Edit12.Text :=  FloattoStr(MyDbf.FieldByName('Cu').AsFloat);
    Edit14.Text :=  FloattoStr(MyDbf.FieldByName('Mo').AsFloat);
    Edit15.Text :=  FloattoStr(MyDbf.FieldByName('Na').AsFloat);
    Edit16.Text :=  FloattoStr(MyDbf.FieldByName('Cl').AsFloat);

    StringGrid1.Cells[1,1] :=  FloattoStr((MyDbf.FieldByName('N').AsFloat/100.0)*(wue*1000)) ;
    StringGrid1.Cells[1,2] :=  FloattoStr((MyDbf.FieldByName('P').AsFloat/100.0)*(wue*1000)) ;
    StringGrid1.Cells[1,3] :=  FloattoStr((MyDbf.FieldByName('K').AsFloat/100.0)*(wue*1000)) ;
    StringGrid1.Cells[1,4] :=  FloattoStr((MyDbf.FieldByName('Mg').AsFloat/100.0)*(wue*1000));
    StringGrid1.Cells[1,5] :=  FloattoStr((MyDbf.FieldByName('Ca').AsFloat/100.0)*(wue*1000))   ;
    StringGrid1.Cells[1,6] :=  FloattoStr((MyDbf.FieldByName('S').AsFloat/100.0)*(wue*1000));
    StringGrid1.Cells[1,7] :=  FloattoStr((MyDbf.FieldByName('Fe').AsFloat/(10000*100.0))*(wue*1000))   ;
    StringGrid1.Cells[1,8] :=  FloattoStr((MyDbf.FieldByName('Mn').AsFloat/(10000*100.0))*(wue*1000)) ;
    StringGrid1.Cells[1,9] :=  FloattoStr((MyDbf.FieldByName('Zn').AsFloat/(10000*100.0))*(wue*1000)) ;
    StringGrid1.Cells[1,10] :=  FloattoStr((MyDbf.FieldByName('B').AsFloat/(10000*100.0))*(wue*1000));
    StringGrid1.Cells[1,11] :=  FloattoStr((MyDbf.FieldByName('Cu').AsFloat/(10000*100.0))*(wue*1000)) ;
    StringGrid1.Cells[1,12] :=  FloattoStr((MyDbf.FieldByName('Si').AsFloat/100.0)*(wue*1000)) ;
    StringGrid1.Cells[1,13] :=  FloattoStr((MyDbf.FieldByName('Mo').AsFloat/(10000*100.0))*(wue*1000)) ;
    StringGrid1.Cells[1,14] :=  FloattoStr((MyDbf.FieldByName('Na').AsFloat/(10000*100.0))*(wue*1000));
    StringGrid1.Cells[1,15] :=  FloattoStr((MyDbf.FieldByName('Cl').AsFloat/(10000*100.0))*(wue*1000)) ;

    MyDbf.Close ;
    MyDbf.Free ;

    Button1.Enabled := False ;
    Button2.Enabled := True ;
    Button4.Enabled := True ;
    Button5.Enabled := True ;
    Button3.Enabled := True ;

end ;

if item_selected = False then
begin
    Button1.Enabled := True ;
    Button2.Enabled := False ;
    Button4.Enabled := True ;
    Button5.Enabled := False ;
    Button3.Enabled := False ;

    for i := 1 to StringGrid1.RowCount - 1 do StringGrid1.Cells[1,i] := '';

    Edit1.Text := '0';
    Edit3.Text := '0';
    Edit4.Text := '0';
    Edit5.Text := '0';
    Edit6.Text := '0';
    Edit7.Text := '0';
    Edit8.Text := '0';
    Edit9.Text := '0';
    Edit10.Text := '0';
    Edit11.Text := '0';
    Edit12.Text := '0';
    Edit13.Text := '0';
    Edit14.Text := '0';
    Edit15.Text := '0';
    Edit16.Text := '0';
    Edit25.Text := 'Name for DB';
end;


end;

procedure TForm16.Button1Click(Sender: TObject);
   var
MyDbf: TDbf;
begin

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := Form1.tissue_analysis_db;
MyDbf.Open             ;
MyDbf.Active := true ;

MyDbf.Insert ;

MyDbf.FieldByName('Name').AsString:= Edit25.Text ;

MyDbf.FieldByName('N').AsFloat:=StrtoFloat(Edit1.Text);
MyDbf.FieldByName('P').AsFloat:=StrtoFloat(Edit3.Text);
MyDbf.FieldByName('K').AsFloat:=StrtoFloat(Edit4.Text);
MyDbf.FieldByName('Mg').AsFloat:=StrtoFloat(Edit5.Text);
MyDbf.FieldByName('Ca').AsFloat:=StrtoFloat(Edit6.Text);
MyDbf.FieldByName('S').AsFloat:=StrtoFloat(Edit7.Text);
MyDbf.FieldByName('Fe').AsFloat:=StrtoFloat(Edit8.Text);
MyDbf.FieldByName('Mn').AsFloat:=StrtoFloat(Edit9.Text);
MyDbf.FieldByName('Zn').AsFloat:=StrtoFloat(Edit10.Text);
MyDbf.FieldByName('B').AsFloat:=StrtoFloat(Edit11.Text);
MyDbf.FieldByName('Cu').AsFloat:=StrtoFloat(Edit12.Text);
MyDbf.FieldByName('Si').AsFloat:=StrtoFloat(Edit13.Text);
MyDbf.FieldByName('Mo').AsFloat:=StrtoFloat(Edit14.Text);
MyDbf.FieldByName('Na').AsFloat:=StrtoFloat(Edit15.Text);
MyDbf.FieldByName('Cl').AsFloat:=StrtoFloat(Edit16.Text);

MyDbf.Post ;

MyDbf.Close ;

MyDbf.Free ;

Form16.UpdateTissueList ;
end;

procedure TForm16.Button2Click(Sender: TObject);
   var
   MyDbf: TDbf;
   i : integer ;
   selected_item : integer ;
   begin

   MyDbf := TDbf.Create(nil) ;
   MyDbf.FilePathFull := '';
   MyDbf.TableName := Form1.tissue_analysis_db;
   MyDbf.Open             ;
   MyDbf.Active := true ;

   if ListBox1.SelCount = 0 then // No ítems selected
        Exit;

   For i := 0 to ListBox1.Items.Count - 1 do

            begin
                                  if ListBox1.Selected [i] then
                                  begin
                                   selected_item := i ;
                                  end;
            end;

   MyDbf.Filter := 'Name=' + QuotedStr(ListBox1.Items[selected_item]) ;

       MyDbf.Filtered := true;       // This selects the filtered set
       MyDbf.First;                  // moves the the first filtered data
       ShowMessage('Deleting ' + MyDbf.FieldByName('Name').AsString + ' from database');
       MyDbf.Delete ;

   MyDbf.Close ;

   MyDbf.Free ;

   ListBox1.Items.Delete(selected_item);

   end;

procedure TForm16.Button3Click(Sender: TObject);
   var
     i,selected_idx : integer ;
     item_selected : boolean ;
     MyDbf: TDbf;
     wue: double;
   begin

   item_selected := false ;
   wue := StrToFloat(Edit17.Text)  ;

   for i := 0 to ListBox1.Items.Count - 1 do

       begin

       if (ListBox1.Selected [i]) then
       begin
            selected_idx := i;
            item_selected := true ;
       end;

       end ;

   if item_selected then

   begin

       MyDbf := TDbf.Create(nil) ;
       MyDbf.FilePathFull := '';
       MyDbf.TableName := Form1.tissue_analysis_db;
       MyDbf.Open             ;
       MyDbf.Active := true ;
       MyDbf.Filter := 'Name=' + QuotedStr(ListBox1.Items[selected_idx]) ;
       MyDbf.Filtered := true;       // This selects the filtered set
       MyDbf.First;                  // moves the the first filtered data

       HB_Main.Form1.Edit1.Text :=  FloattoStr((MyDbf.FieldByName('N').AsFloat/100.0)*(wue*1000)) ;
       HB_Main.Form1.Edit3.Text :=  FloattoStr((MyDbf.FieldByName('P').AsFloat/100.0)*(wue*1000)) ;
       HB_Main.Form1.Edit4.Text :=  FloattoStr((MyDbf.FieldByName('K').AsFloat/100.0)*(wue*1000)) ;
       HB_Main.Form1.Edit5.Text :=  FloattoStr((MyDbf.FieldByName('Mg').AsFloat/100.0)*(wue*1000));
       HB_Main.Form1.Edit6.Text :=  FloattoStr((MyDbf.FieldByName('Ca').AsFloat/100.0)*(wue*1000))   ;
       HB_Main.Form1.Edit7.Text :=  FloattoStr((MyDbf.FieldByName('S').AsFloat/100.0)*(wue*1000));
       HB_Main.Form1.Edit13.Text :=  FloattoStr((MyDbf.FieldByName('Si').AsFloat/100.0)*(wue*1000)) ;
       HB_Main.Form1.Edit8.Text :=  FloattoStr((MyDbf.FieldByName('Fe').AsFloat/(10000*100.0))*(wue*1000))   ;
       HB_Main.Form1.Edit9.Text :=  FloattoStr((MyDbf.FieldByName('Mn').AsFloat/(10000*100.0))*(wue*1000)) ;
       HB_Main.Form1.Edit10.Text :=  FloattoStr((MyDbf.FieldByName('Zn').AsFloat/(10000*100.0))*(wue*1000)) ;
       HB_Main.Form1.Edit11.Text :=  FloattoStr((MyDbf.FieldByName('B').AsFloat/(10000*100.0))*(wue*1000));
       HB_Main.Form1.Edit12.Text :=  FloattoStr((MyDbf.FieldByName('Cu').AsFloat/(10000*100.0))*(wue*1000)) ;
       HB_Main.Form1.Edit14.Text :=  FloattoStr((MyDbf.FieldByName('Mo').AsFloat/(10000*100.0))*(wue*1000)) ;
       HB_Main.Form1.Edit15.Text :=  FloattoStr((MyDbf.FieldByName('Na').AsFloat/(10000*100.0))*(wue*1000));
       HB_Main.Form1.Edit16.Text :=  FloattoStr((MyDbf.FieldByName('Cl').AsFloat/(10000*100.0))*(wue*1000)) ;

       MyDbf.Close ;
       MyDbf.Free ;

   end ;

   if item_selected = False then
   begin
       Button1.Enabled := True ;
       Button2.Enabled := False ;
       Button4.Enabled := True ;
       Button5.Enabled := False ;
       Button3.Enabled := False ;
       for i := 1 to StringGrid1.RowCount - 1 do StringGrid1.Cells[1,i] := '';
       Edit1.Text := '0';
       Edit3.Text := '0';
        Edit4.Text := '0';
        Edit5.Text := '0';
        Edit6.Text := '0';
        Edit7.Text := '0';
        Edit8.Text := '0';
        Edit9.Text := '0';
        Edit10.Text := '0';
        Edit11.Text := '0';
        Edit12.Text := '0';
        Edit13.Text := '0';
        Edit14.Text := '0';
        Edit15.Text := '0';
        Edit16.Text := '0';
        Edit25.Text := 'Name for DB';
   end;

end;

procedure TForm16.Button4Click(Sender: TObject);
var
  i: integer;
begin
    ListBox1.Clear;
    UpdateTissueList;
    Button1.Enabled := True ;
    Button2.Enabled := False ;
    Button4.Enabled := True ;
    Button5.Enabled := False ;
    Button3.Enabled := False ;
    for i := 1 to StringGrid1.RowCount - 1 do StringGrid1.Cells[1,i] := '';
    Edit1.Text := '0';
    Edit3.Text := '0';
    Edit4.Text := '0';
    Edit5.Text := '0';
    Edit6.Text := '0';
    Edit7.Text := '0';
    Edit8.Text := '0';
    Edit9.Text := '0';
    Edit10.Text := '0';
    Edit11.Text := '0';
    Edit12.Text := '0';
    Edit13.Text := '0';
    Edit14.Text := '0';
    Edit15.Text := '0';
    Edit16.Text := '0';
    Edit25.Text := 'Name for DB';
end;

procedure TForm16.Button5Click(Sender: TObject);
   var
     MyDbf: TDbf;
     selected_item, i : integer ;
   begin

   MyDbf := TDbf.Create(nil) ;
   MyDbf.FilePathFull := '';
   MyDbf.TableName := Form1.tissue_analysis_db;
   MyDbf.Open             ;
   MyDbf.Active := true ;

   if ListBox1.SelCount = 0 then // No ítems selected
        Exit;

   For i := 0 to ListBox1.Items.Count - 1 do

            begin
                                  if ListBox1.Selected [i] then
                                  begin
                                   selected_item := i ;
                                  end;
            end;

   MyDbf.Filter := 'Name=' + QuotedStr(ListBox1.Items[selected_item]) ;

       MyDbf.Filtered := true;       // This selects the filtered set
       MyDbf.First;

       MyDbf.Edit;

       MyDbf.FieldByName('Name').AsString:= Edit25.Text ;
       MyDbf.FieldByName('N').AsFloat:=StrtoFloat(Edit1.Text);
       MyDbf.FieldByName('P').AsFloat:=StrtoFloat(Edit3.Text);
       MyDbf.FieldByName('K').AsFloat:=StrtoFloat(Edit4.Text);
       MyDbf.FieldByName('Mg').AsFloat:=StrtoFloat(Edit5.Text);
       MyDbf.FieldByName('Ca').AsFloat:=StrtoFloat(Edit6.Text);
       MyDbf.FieldByName('S').AsFloat:=StrtoFloat(Edit7.Text);
       MyDbf.FieldByName('Fe').AsFloat:=StrtoFloat(Edit8.Text);
       MyDbf.FieldByName('Mn').AsFloat:=StrtoFloat(Edit9.Text);
       MyDbf.FieldByName('Zn').AsFloat:=StrtoFloat(Edit10.Text);
       MyDbf.FieldByName('B').AsFloat:=StrtoFloat(Edit11.Text);
       MyDbf.FieldByName('Cu').AsFloat:=StrtoFloat(Edit12.Text);
       MyDbf.FieldByName('Si').AsFloat:=StrtoFloat(Edit13.Text);
       MyDbf.FieldByName('Mo').AsFloat:=StrtoFloat(Edit14.Text);
       MyDbf.FieldByName('Na').AsFloat:=StrtoFloat(Edit15.Text);
       MyDbf.FieldByName('Cl').AsFloat:=StrtoFloat(Edit16.Text);

       MyDbf.Post ;

   MyDbf.Close ;

   MyDbf.Free ;

   UpdateTissueList ;




   end;

initialization
  {$I hb_tissue_analysis.lrs}

end.

