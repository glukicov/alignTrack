from ROOT import TStyle, gROOT, gStyle, TColor 

def SetMyStyle():
  print("\n ~/rootlogon.C loaded with !!4 sig.fig.!! for custom Opt Fit and Stat!\n")
  MyStyle = myStyle()
  gROOT.SetStyle("MyStyle")
  gROOT.ForceStyle()


def myStyle():
  myStyle  = TStyle("MyStyle", "My Root Styles")
  #Canvas
  myStyle.SetCanvasBorderMode(0)  # Transparent
  myStyle.SetCanvasColor(0) # Transparent 
  #Paper, Pad, Palette, Frame
  myStyle.SetPadBorderMode(0) # Transparent 
  myStyle.SetPadColor(0) # Transparent 
  myStyle.SetPalette(1) # Default 
  myStyle.SetFrameBorderMode(1) # Border
   # Axis 
  myStyle.SetLabelSize(0.04, "xyz") # size of axis values
  myStyle.SetTitleSize(0.04, "xyz")
  myStyle.SetPadTickX(1)
  myStyle.SetPadTickY(1)
  # Title 
  myStyle.SetTitleColor(1) # Black 
  myStyle.SetTitleStyle(0) # Transparent 
  myStyle.SetTitleBorderSize(0) # Transparent
  myStyle.SetTitleY(0.97) # Set y-position (fraction of pad size)
  myStyle.SetTitleX(0.4) # Set x-position (fraction of pad size)
  # #Stat box dimensions, position and style 
  myStyle.SetStatY(0.89) # Set y-position (fraction of pad size)
  myStyle.SetStatX(0.89) # Set x-position (fraction of pad size)
  myStyle.SetStatW(0.36) # Set width of stat-box (fraction of pad size)
  myStyle.SetStatH(0.12) # Set height of stat-box (fraction of pad size)
  myStyle.SetStatStyle(0) # Transparent 
  myStyle.SetStatColor(0)  # Transparent
  myStyle.SetStatBorderSize(1) # Transparent
  # Histo Filling (visual)
  myStyle.SetHistFillColor(3)
  myStyle.SetHistLineColor(3)
  myStyle.SetHistFillStyle(1001)           
  # Stats display options 
  #myStyle.SetOptStat("ourRmMe") #over/under -flows, Rms and Means with errors, number of entries
  myStyle.SetOptStat("neouRM") #over/under -flows, Rms and Means with errors, number of entries
  myStyle.SetOptFit(1111)  #probability, Chi2, errors, name/values of parameters
  myStyle.SetStatFormat("11.4f")  # 4 sig.fig, f=float

  return myStyle