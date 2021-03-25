unit hb_addweight;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Dbf, db, Dbf_Common ;

type

  { TForm4 }

  TForm4 = class(TForm)
    Button1: TButton;
    Edit1: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
    is_liquid: integer;
  end; 

var
  Form4: TForm4; 

implementation

{ TForm4 }

uses HB_Main;

procedure TForm4.Button1Click(Sender: TObject);
var
MyDbf: TDbf;
begin

MyDbf := TDbf.Create(nil) ;
MyDbf.FilePathFull := '';
MyDbf.TableName := substances_used_db;
MyDbf.Open             ;
MyDbf.Active := true ;

MyDbf.Filter := 'Name=' + QuotedStr(Label2.Caption) ;

    MyDbf.Filtered := true;       // This selects the filtered set
    MyDbf.First;

    MyDbf.Edit;

               MyDbf.FieldByName('Weight').AsFloat:= StrtoFloat(Edit1.Text) ;

    MyDbf.Post ;

MyDbf.Close ;

MyDbf.Free ;

Form4.Visible := False ;



end;

initialization
  {$I hb_addweight.lrs}

end.

