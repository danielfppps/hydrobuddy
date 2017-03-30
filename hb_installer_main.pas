unit hb_installer_main;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, ComCtrls,
  StdCtrls, ExtCtrls, Buttons, Zipper, Windows,LCLIntf;

type

  { TForm1 }

  TForm1 = class(TForm)
    BitBtn1: TBitBtn;
    BitBtn2: TBitBtn;
    BitBtn3: TBitBtn;
    BitBtn4: TBitBtn;
    BitBtn5: TBitBtn;
    BitBtn6: TBitBtn;
    Button1: TBitBtn;
    Button2: TButton;
    Edit1: TEdit;
    Image1: TImage;
    Image2: TImage;
    Image3: TImage;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    PageControl1: TPageControl;
    ProgressBar1: TProgressBar;
    SelectDirectoryDialog1: TSelectDirectoryDialog;
    TabSheet1: TTabSheet;
    TabSheet2: TTabSheet;
    TabSheet3: TTabSheet;
    TabSheet4: TTabSheet;
    procedure BitBtn1Click(Sender: TObject);
    procedure BitBtn2Click(Sender: TObject);
    procedure BitBtn3Click(Sender: TObject);
    procedure BitBtn4Click(Sender: TObject);
    procedure BitBtn5Click(Sender: TObject);
    procedure BitBtn6Click(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end;

var
  Form1: TForm1;

implementation

{$R *.lfm}
{$R hb_install.rc}

{ TForm1 }

procedure TForm1.Button1click(Sender: TObject);
begin
    PageControl1.ActivePage := TabSheet1;
end;

procedure TForm1.Button2Click(Sender: TObject);
begin

    if SelectDirectoryDialog1.Execute then
    Edit1.Text := SelectDirectoryDialog1.FileName + '\' + 'HydroBuddy';
end;

procedure TForm1.BitBtn2Click(Sender: TObject);
begin
  PageControl1.ActivePage := TabSheet4;
end;

procedure TForm1.BitBtn3Click(Sender: TObject);
  var
  S: TResourceStream;
  F: TFileStream;
  UnZipper: TUnZipper;
begin

  ProgressBar1.Position := 1;

  if DirectoryExists(Edit1.Text) = false then
  CreateDir(Edit1.Text);

  ProgressBar1.Position := 2;
  // create a resource stream which points to our resource
  S := TResourceStream.Create(HInstance, 'HB_DATA', RT_RCDATA);
  // Please be aware of writing an apostrophes in resource type - source will not be axtracted!!!
  try
    // create a file mydata.dat in the application directory
    F := TFileStream.Create(Edit1.Text + '/hb_install.zip', fmCreate);
    try
      F.CopyFrom(S, S.Size); // copy data from the resource stream to file stream
    finally
      F.Free; // destroy the file stream
    end;
  finally
    S.Free; // destroy the resource stream
  end;

  ProgressBar1.Position := 3;

  SetCurrentDir(Edit1.Text);

  UnZipper := TUnZipper.Create;
  try
    UnZipper.FileName := Edit1.Text + '/hb_install.zip';
    UnZipper.OutputPath := Edit1.Text;
    UnZipper.Examine;
    UnZipper.UnZipAllFiles;
  finally
    UnZipper.Free;
  end;

  ProgressBar1.Position := 4;

  PageControl1.ActivePage := TabSheet3;


end;

procedure TForm1.BitBtn4Click(Sender: TObject);
begin
 PageControl1.ActivePage := TabSheet1;
end;

procedure TForm1.BitBtn5Click(Sender: TObject);
begin
  OpenURL('https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6YR6X5AAEGBGJ');

end;

procedure TForm1.BitBtn6Click(Sender: TObject);
begin
  Application.Terminate;
end;

procedure TForm1.BitBtn1Click(Sender: TObject);
begin
  if Edit1.Text <> 'Input your installation path here' then
  PageControl1.ActivePage := TabSheet2 else
    ShowMessage('Please select an installation folder before proceeding.');
end;

end.

