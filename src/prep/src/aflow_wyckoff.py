# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOWpi software.
#
#  AFLOWpi is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************



class space_group_sym():
  def __init__(self,sg_num):
    self.spacegroup="P1";
    self.spacegrouplabel="#1";
    self.spacegroupnumber=1;
    self.spacegroupoption="";
    self.sg_num=sg_num
    self.sym_ops=[]
    self.option=1
    self.origin = np.array([0.0,0.0,0.0])
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 1  P1 #1
    if(self.sg_num==1) :
      self.spacegroup="P1"
      self.spacegrouplabel="#1"
      self.spacegroupnumber=1
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 2 P-1 #2
    if(self.sg_num==2) :
      self.spacegroup="P-1"
      self.spacegrouplabel="#2"
      self.spacegroupnumber=2
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 3  P2 #3
    if(self.sg_num==3 and self.option==1) :
      self.spacegroup="P2"
      self.spacegrouplabel="#3"
      self.spacegroupnumber=3
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y,-z")
    
    if(self.sg_num==3 and self.option==2) :
      self.spacegroup="P2"
      self.spacegrouplabel="#3"
      self.spacegroupnumber=3
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 4  P2_:1 #4
    if(self.sg_num==4 and self.option==1) :
      self.spacegroup="P2_:1"
      self.spacegrouplabel="#4"
      self.spacegroupnumber=4
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y+1/2,-z")
    
    if(self.sg_num==4 and self.option==2) :
      self.spacegroup="P2_:1"
      self.spacegrouplabel="#4"
      self.spacegroupnumber=4
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 5  C2 #5
    if(self.sg_num==5 and self.option==1) :
      self.spacegroup="C2"
      self.spacegrouplabel="#5"
      self.spacegroupnumber=5
      self.spacegroupoption="unique axis b"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y,-z")
      
    
    if(self.sg_num==5 and self.option==2) :
      self.spacegroup="C2"
      self.spacegrouplabel="#5"
      self.spacegroupnumber=5
      self.spacegroupoption="unique axis c"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 6  Pm #6
    if(self.sg_num==6 and self.option==1) :
      self.spacegroup="Pm"
      self.spacegrouplabel="#6"
      self.spacegroupnumber=6
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x,-y,z")
    
    if(self.sg_num==6 and self.option==2) :
      self.spacegroup="Pm"
      self.spacegrouplabel="#6"
      self.spacegroupnumber=6
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x,y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 7  Pc #7
    if(self.sg_num==7 and self.option==1) :
      self.spacegroup="Pc"
      self.spacegrouplabel="#7"
      self.spacegroupnumber=7
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x,-y,z+1/2")
    
    if(self.sg_num==7 and self.option==2) :
      self.spacegroup="Pc"
      self.spacegrouplabel="#7"
      self.spacegroupnumber=7
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x+1/2,y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 8  Cm #8
    if(self.sg_num==8 and self.option==1) :
      self.spacegroup="Cm"
      self.spacegrouplabel="#8"
      self.spacegroupnumber=8
      self.spacegroupoption="unique axis b"

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x,-y,z")
      
    
    if(self.sg_num==8 and self.option==2) :
      self.spacegroup="Cm"
      self.spacegrouplabel="#8"
      self.spacegroupnumber=8
      self.spacegroupoption="unique axis c"


      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x,y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 9  Cc #9
    if(self.sg_num==9 and self.option==1) :
      self.spacegroup="Cc"
      self.spacegrouplabel="#9"
      self.spacegroupnumber=9
      self.spacegroupoption="unique axis b"

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x,-y,z+1/2")
      
    
    if(self.sg_num==9 and self.option==2) :
      self.spacegroup="Cc"
      self.spacegrouplabel="#9"
      self.spacegroupnumber=9
      self.spacegroupoption="unique axis c"

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("x+1/2,y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 10  P2/m #10
    if(self.sg_num==10 and self.option==1) :
      self.spacegroup="P2/m"
      self.spacegrouplabel="#10"
      self.spacegroupnumber=10
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,-y,z")
    
    if(self.sg_num==10 and self.option==2) :
      self.spacegroup="P2/m"
      self.spacegrouplabel="#10"
      self.spacegroupnumber=10
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 11  P2_:1/m #11
    if(self.sg_num==11 and self.option==1) :
      self.spacegroup="P2_:1/m"
      self.spacegrouplabel="#11"
      self.spacegroupnumber=11
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,-y+1/2,z")
    
    if(self.sg_num==11 and self.option==2) :
      self.spacegroup="P2_:1/m"
      self.spacegrouplabel="#11"
      self.spacegroupnumber=11
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 12  C2/m #12
    if(self.sg_num==12 and self.option==1) :
      self.spacegroup="C2/m"
      self.spacegrouplabel="#12"
      self.spacegroupnumber=12
      self.spacegroupoption="unique axis b"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,-y,z")
      
    
    if(self.sg_num==12 and self.option==2) :
      self.spacegroup="C2/m"
      self.spacegrouplabel="#12"
      self.spacegroupnumber=12
      self.spacegroupoption="unique axis c"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 13  P2/c #13
    if(self.sg_num==13 and self.option==1) :
      self.spacegroup="P2/c"
      self.spacegrouplabel="#13"
      self.spacegroupnumber=13
      self.spacegroupoption="unique axis b"
      
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,-y,z+1/2")
    
    if(self.sg_num==13 and self.option==2) :
      self.spacegroup="P2/c"
      self.spacegrouplabel="#13"
      self.spacegroupnumber=13
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 14  P2_:1/c #14
    if(self.sg_num==14 and self.option==1) :
      self.spacegroup="P2_:1/c"
      self.spacegrouplabel="#14"
      self.spacegroupnumber=14
      self.spacegroupoption="unique axis b"
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,-y+1/2,z+1/2")
    
    if(self.sg_num==14 and self.option==2) :
      self.spacegroup="P2_:1/c"
      self.spacegrouplabel="#14"
      self.spacegroupnumber=14
      self.spacegroupoption="unique axis c"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 15  C2/c #15
    if(self.sg_num==15 and self.option==1) :
      self.spacegroup="C2/c"
      self.spacegrouplabel="#15"
      self.spacegroupnumber=15
      self.spacegroupoption="unique axis b"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,-y,z+1/2")
      
    
    if(self.sg_num==15 and self.option==2) :
      self.spacegroup="C2/c"
      self.spacegrouplabel="#15"
      self.spacegroupnumber=15
      self.spacegroupoption="unique axis c"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 16  P222 #16
    if(self.sg_num==16) :
      self.spacegroup="P222"
      self.spacegrouplabel="#16"
      self.spacegroupnumber=16
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 17  P222_:1 #17
    if(self.sg_num==17) :
      self.spacegroup="P222_:1"
      self.spacegrouplabel="#17"
      self.spacegroupnumber=17
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 18  P2_:12_:12 #18
    if(self.sg_num==18) :
      self.spacegroup="P2_:12_:12"
      self.spacegrouplabel="#18"
      self.spacegroupnumber=18
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 19  P2_:12_:12_:1 #19
    if(self.sg_num==19) :
      self.spacegroup="P2_:12_:12_:1"
      self.spacegrouplabel="#19"
      self.spacegroupnumber=19
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 20  C222_:1 #20
    if(self.sg_num==20) :
      self.spacegroup="C222_:1"
      self.spacegrouplabel="#20"
      self.spacegroupnumber=20
      self.spacegroupoption=""	
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 21  C222 #21
    if(self.sg_num==21) :
      self.spacegroup="C222"
      self.spacegrouplabel="#21"
      self.spacegroupnumber=21
      self.spacegroupoption=""
       	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 22  F222 #22
    if(self.sg_num==22) :
      self.spacegroup="F222"
      self.spacegrouplabel="#22"
      self.spacegroupnumber=22
      self.spacegroupoption=""
      	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 23  I222 #23
    if(self.sg_num==23) :
      self.spacegroup="I222"
      self.spacegrouplabel="#23"
      self.spacegroupnumber=23
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 24  I2_:12_:12_:1 #24
    if(self.sg_num==24) :
      self.spacegroup="I2_:12_:12_:1"
      self.spacegrouplabel="#24"
      self.spacegroupnumber=24
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 25  Pmm2 #25
    if(self.sg_num==25) :
      self.spacegroup="Pmm2"
      self.spacegrouplabel="#25"
      self.spacegroupnumber=25
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 26  Pmc2_:1 #26
    if(self.sg_num==26) :
      self.spacegroup="Pmc2_:1"
      self.spacegrouplabel="#26"
      self.spacegroupnumber=26
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 27  Pcc2 #27
    if(self.sg_num==27) :
      self.spacegroup="Pcc2"
      self.spacegrouplabel="#27"
      self.spacegroupnumber=27
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 28  Pma2 #28
    if(self.sg_num==28) :
      self.spacegroup="Pma2"
      self.spacegrouplabel="#28"
      self.spacegroupnumber=28
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x+1/2,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 29  Pca2_:1 #29
    if(self.sg_num==29) :
      self.spacegroup="Pca2_:1"
      self.spacegrouplabel="#29"
      self.spacegroupnumber=29
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x+1/2,y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 30  Pnc2 #30
    if(self.sg_num==30) :
      self.spacegroup="Pnc2"
      self.spacegrouplabel="#30"
      self.spacegroupnumber=30
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 31  Pmn2_:1 #31
    if(self.sg_num==31) :
      self.spacegroup="Pmn2_:1"
      self.spacegrouplabel="#31"
      self.spacegroupnumber=31
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 32  Pba2 #32
    if(self.sg_num==32) :
      self.spacegroup="Pba2"
      self.spacegrouplabel="#32"
      self.spacegroupnumber=32
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 33  Pna2_:1 #33
    if(self.sg_num==33) :
      self.spacegroup="Pna2_:1"
      self.spacegrouplabel="#33"
      self.spacegroupnumber=33
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 34  Pnn2 #34
    if(self.sg_num==34) :
      self.spacegroup="Pnn2"
      self.spacegrouplabel="#34"
      self.spacegroupnumber=34
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 35  Cmm2 #35
    if(self.sg_num==35) :
      self.spacegroup="Cmm2"
      self.spacegrouplabel="#35"
      self.spacegroupnumber=35
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 36  Cmc2_:1 #36
    if(self.sg_num==36) :
      self.spacegroup="Cmc2_:1"
      self.spacegrouplabel="#36"
      self.spacegroupnumber=36
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 37  Ccc2 #37
    if(self.sg_num==37) :
      self.spacegroup="Ccc2"
      self.spacegrouplabel="#37"
      self.spacegroupnumber=37
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 38  Amm2 #38
    if(self.sg_num==38) :
      self.spacegroup="Amm2"
      self.spacegrouplabel="#38"
      self.spacegroupnumber=38
      self.spacegroupoption=""
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 39  Aem2 #39
    if(self.sg_num==39) :
      self.spacegroup="Aem2"
      self.spacegrouplabel="#39"
      self.spacegroupnumber=39
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x,y+1/2,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 40  Ama2 #40
    if(self.sg_num==40) :
      self.spacegroup="Ama2"
      self.spacegrouplabel="#40"
      self.spacegroupnumber=40
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x+1/2,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 41  Aea2 #41
    if(self.sg_num==41) :
      self.spacegroup="Aea2"
      self.spacegrouplabel="#41"
      self.spacegroupnumber=41
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 42  Fmm2 #42
    if(self.sg_num==42) :
      self.spacegroup="Fmm2"
      self.spacegrouplabel="#42"
      self.spacegroupnumber=42
      self.spacegroupoption=""
      	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 43  Fdd2 #43
    if(self.sg_num==43) :
      self.spacegroup="Fdd2"
      self.spacegrouplabel="#43"
      self.spacegroupnumber=43
      self.spacegroupoption=""
      				
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/4,-y+1/4,z+1/4")
      self.sym_ops.append("-x+1/4,y+1/4,z+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 44  Imm2 #44
    if(self.sg_num==44) :
      self.spacegroup="Imm2"
      self.spacegrouplabel="#44"
      self.spacegroupnumber=44
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 45  Iba2 #45
    if(self.sg_num==45) :
      self.spacegroup="Iba2"
      self.spacegrouplabel="#45"
      self.spacegroupnumber=45
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 46  Ima2 #46
    if(self.sg_num==46) :
      self.spacegroup="Ima2"
      self.spacegrouplabel="#46"
      self.spacegroupnumber=46
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x+1/2,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 47  Pmmm #47
    if(self.sg_num==47) :
      self.spacegroup="Pmmm"
      self.spacegrouplabel="#47"
      self.spacegroupnumber=47
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 48  Pnnn #48
    if(self.sg_num==48 and self.option==1) :
      self.spacegroup="Pnnn"
      self.spacegrouplabel="#48"
      self.spacegroupnumber=48
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
    
    if(self.sg_num==48 and self.option==2) :
      self.spacegroup="Pnnn"
      self.spacegrouplabel="#48"
      self.spacegroupnumber=48
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 49  Pccm #49
    if(self.sg_num==49) :
      self.spacegroup="Pccm"
      self.spacegrouplabel="#49"
      self.spacegroupnumber=49
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 50  Pban #50
    if(self.sg_num==50 and self.option==1) :
      self.spacegroup="Pban"
      self.spacegrouplabel="#50"
      self.spacegroupnumber=50
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
    
    if(self.sg_num==50 and self.option==2) :
      self.spacegroup="Pban"
      self.spacegrouplabel="#50"
      self.spacegroupnumber=50
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z")
      self.sym_ops.append("x,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x,y+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 51  Pmma #51
    if(self.sg_num==51) :
      self.spacegroup="Pmma"
      self.spacegrouplabel="#51"
      self.spacegroupnumber=51
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x+1/2,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x+1/2,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 52  Pnna #52
    if(self.sg_num==52) :
      self.spacegroup="Pnna"
      self.spacegrouplabel="#52"
      self.spacegroupnumber=52
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 53  Pmna #53
    if(self.sg_num==53) :
      self.spacegroup="Pmna"
      self.spacegrouplabel="#53"
      self.spacegroupnumber=53
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 54  Pcca #54
    if(self.sg_num==54) :
      self.spacegroup="Pcca"
      self.spacegrouplabel="#54"
      self.spacegroupnumber=54
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x+1/2,y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 55  Pbam #55
    if(self.sg_num==55) :
      self.spacegroup="Pbam"
      self.spacegrouplabel="#55"
      self.spacegroupnumber=55
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 56  Pccn #56
    if(self.sg_num==56) :
      self.spacegroup="Pccn"
      self.spacegrouplabel="#56"
      self.spacegroupnumber=56
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 57  Pbcm #57
    if(self.sg_num==57) :
      self.spacegroup="Pbcm"
      self.spacegrouplabel="#57"
      self.spacegroupnumber=57
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 58  Pnnm #58
    if(self.sg_num==58) :
      self.spacegroup="Pnnm"
      self.spacegrouplabel="#58"
      self.spacegroupnumber=58
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 59  Pmmn #59
    if(self.sg_num==59 and self.option==1) :
      self.spacegroup="Pmmn"
      self.spacegrouplabel="#59"
      self.spacegroupnumber=59
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
    
    if(self.sg_num==59 and self.option==2) :
      self.spacegroup="Pmmn"
      self.spacegrouplabel="#59"
      self.spacegroupnumber=59
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 60  Pbcn #60
    if(self.sg_num==60) :
      self.spacegroup="Pbcn"
      self.spacegrouplabel="#60"
      self.spacegroupnumber=60
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 61  Pbca #61
    if(self.sg_num==61) :
      self.spacegroup="Pbca"
      self.spacegrouplabel="#61"
      self.spacegroupnumber=61
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 62  Pnma #62
    if(self.sg_num==62) :
      self.spacegroup="Pnma"
      self.spacegrouplabel="#62"
      self.spacegroupnumber=62
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 63  Cmcm #63
    if(self.sg_num==63) :
      self.spacegroup="Cmcm"
      self.spacegrouplabel="#63"
      self.spacegroupnumber=63
      self.spacegroupoption=""

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z")
        
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 64  Cmce #64
    if(self.sg_num==64) :
      self.spacegroup="Cmce"
      self.spacegrouplabel="#64"
      self.spacegroupnumber=64
      self.spacegroupoption=""

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y+1/2,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 65  Cmmm #65
    if(self.sg_num==65) :
      self.spacegroup="Cmmm"
      self.spacegrouplabel="#65"
      self.spacegroupnumber=65
      self.spacegroupoption=""
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 66  Cccm #66
    if(self.sg_num==66) :
      self.spacegroup="Cccm"
      self.spacegrouplabel="#66"
      self.spacegroupnumber=66
      self.spacegroupoption=""
		
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 67  Cmme #67
    if(self.sg_num==67) :
      self.spacegroup="Cmme"
      self.spacegrouplabel="#67"
      self.spacegroupnumber=67
      self.spacegroupoption=""
		
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y+1/2,z")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y+1/2,-z")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 68  Ccce #68
    if(self.sg_num==68 and self.option==1) :
      self.spacegroup="Ccce"
      self.spacegrouplabel="#68"
      self.spacegroupnumber=68
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y,z+1/2")
      
    
    if(self.sg_num==68 and self.option==2) :
      self.spacegroup="Ccce"
      self.spacegrouplabel="#68"
      self.spacegroupnumber=68
      self.spacegroupoption="origin choice 2"

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x+1/2,y,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 69  Fmmm #69
    if(self.sg_num==69) :
      self.spacegroup="Fmmm"
      self.spacegrouplabel="#69"
      self.spacegroupnumber=69
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 70  Fddd #70
    if(self.sg_num==70 and self.option==1) :
      self.spacegroup="Fddd"
      self.spacegrouplabel="#70"
      self.spacegroupnumber=70
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x+1/4,-y+1/4,-z+1/4")
      self.sym_ops.append("x+1/4,y+1/4,-z+1/4")
      self.sym_ops.append("x+1/4,-y+1/4,z+1/4")
      self.sym_ops.append("-x+1/4,y+1/4,z+1/4")
      
    
    if(self.sg_num==70 and self.option==2) :
      self.spacegroup="Fddd"
      self.spacegrouplabel="#70"
      self.spacegroupnumber=70
      self.spacegroupoption="origin choice 2"
		
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+3/4,-y+3/4,z")
      self.sym_ops.append("-x+3/4,y,-z+3/4")
      self.sym_ops.append("x,-y+3/4,-z+3/4")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/4,y+1/4,-z")
      self.sym_ops.append("x+1/4,-y,z+1/4")
      self.sym_ops.append("-x,y+1/4,z+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 71  Immm #71
    if(self.sg_num==71) :
      self.spacegroup="Immm"
      self.spacegrouplabel="#71"
      self.spacegroupnumber=71
      self.spacegroupoption=""	
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 72  Ibam #72
    if(self.sg_num==72) :
      self.spacegroup="Ibam"
      self.spacegrouplabel="#72"
      self.spacegroupnumber=72
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 73  Ibca #73
    if(self.sg_num==73) :
      self.spacegroup="Ibca"
      self.spacegrouplabel="#73"
      self.spacegroupnumber=73
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 74  Imma #74
    if(self.sg_num==74) :
      self.spacegroup="Imma"
      self.spacegrouplabel="#74"
      self.spacegroupnumber=74
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y+1/2,z")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y+1/2,-z")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x,y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 75  P4 #75
    if(self.sg_num==75) :
      self.spacegroup="P4"
      self.spacegrouplabel="#75"
      self.spacegroupnumber=75
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 76  P4_:1 #76
    if(self.sg_num==76) :
      self.spacegroup="P4_:1"
      self.spacegrouplabel="#76"
      self.spacegroupnumber=76
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-y,x,z+1/4")
      self.sym_ops.append("y,-x,z+3/4")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 77  P4_:2 #77
    if(self.sg_num==77) :
      self.spacegroup="P4_:2"
      self.spacegrouplabel="#77"
      self.spacegroupnumber=77
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 78  P4_:3 #78
    if(self.sg_num==78) :
      self.spacegroup="P4_:3"
      self.spacegrouplabel="#78"
      self.spacegroupnumber=78
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-y,x,z+3/4")
      self.sym_ops.append("y,-x,z+1/4")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 79  I4 #79
    if(self.sg_num==79) :
      self.spacegroup="I4"
      self.spacegrouplabel="#79"
      self.spacegroupnumber=79
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 80  I4_:1 #80
    if(self.sg_num==80) :
      self.spacegroup="I4_:1"
      self.spacegrouplabel="#80"
      self.spacegroupnumber=80
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 81  P-4 #81
    if(self.sg_num==81) :
      self.spacegroup="P-4"
      self.spacegrouplabel="#81"
      self.spacegroupnumber=81
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 82  I-4 #82
    if(self.sg_num==82) :
      self.spacegroup="I-4"
      self.spacegrouplabel="#82"
      self.spacegroupnumber=82
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 83  P4/m #83
    if(self.sg_num==83) :
      self.spacegroup="P4/m"
      self.spacegrouplabel="#83"
      self.spacegroupnumber=83
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 84  P4_:2/m #84
    if(self.sg_num==84) :
      self.spacegroup="P4_:2/m"
      self.spacegrouplabel="#84"
      self.spacegroupnumber=84
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z+1/2")
      self.sym_ops.append("-y,x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 85  P4/n #85
    if(self.sg_num==85 and self.option==1) :
      self.spacegroup="P4/n"
      self.spacegrouplabel="#85"
      self.spacegroupnumber=85
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z")
      self.sym_ops.append("y+1/2,-x+1/2,z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
    
    if(self.sg_num==85 and self.option==2) :
      self.spacegroup="P4/n"
      self.spacegrouplabel="#85"
      self.spacegroupnumber=85
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z")
      self.sym_ops.append("y,-x+1/2,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z")
      self.sym_ops.append("-y,x+1/2,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 86  P4_:2/n #86
    if(self.sg_num==86 and self.option==1) :
      self.spacegroup="P4_:2/n"
      self.spacegrouplabel="#86"
      self.spacegroupnumber=86
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
    
    if(self.sg_num==86 and self.option==2) :
      self.spacegroup="P4_:2/n"
      self.spacegrouplabel="#86"
      self.spacegroupnumber=86
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x,z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y,-x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 87  I4/m #87
    if(self.sg_num==87) :
      self.spacegroup="I4/m"
      self.spacegrouplabel="#87"
      self.spacegroupnumber=87
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 88  I4_:1/a #88
    if(self.sg_num==88 and self.option==1) :
      self.spacegroup="I4_:1/a"
      self.spacegrouplabel="#88"
      self.spacegroupnumber=88
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      self.sym_ops.append("-x,-y+1/2,-z+1/4")
      self.sym_ops.append("x+1/2,y,-z+3/4")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      
    
    if(self.sg_num==88 and self.option==2) :
      self.spacegroup="I4_:1/a"
      self.spacegrouplabel="#88"
      self.spacegroupnumber=88
      self.spacegroupoption="origin choice 2"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-y+3/4,x+1/4,z+1/4")
      self.sym_ops.append("y+3/4,-x+3/4,z+3/4")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("y+1/4,-x+3/4,-z+3/4")
      self.sym_ops.append("-y+1/4,x+1/4,-z+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 89  P422 #89
    if(self.sg_num==89) :
      self.spacegroup="P422"
      self.spacegrouplabel="#89"
      self.spacegroupnumber=89
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 90  P42_:12 #90
    if(self.sg_num==90) :
      self.spacegroup="P42_:12"
      self.spacegrouplabel="#90"
      self.spacegroupnumber=90
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z")
      self.sym_ops.append("y+1/2,-x+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 91  P4_:122 #91
    if(self.sg_num==91) :
      self.spacegroup="P4_:122"
      self.spacegrouplabel="#91"
      self.spacegroupnumber=91
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-y,x,z+1/4")
      self.sym_ops.append("y,-x,z+3/4")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y,x,-z+3/4")
      self.sym_ops.append("-y,-x,-z+1/4")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 92  P4_:12_:12 #92
    if(self.sg_num==92) :
      self.spacegroup="P4_:12_:12"
      self.spacegrouplabel="#92"
      self.spacegroupnumber=92
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x+1/2,z+3/4")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/4")
      self.sym_ops.append("x+1/2,-y+1/2,-z+3/4")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 93  P4_:222 #93
    if(self.sg_num==93) :
      self.spacegroup="P4_:222"
      self.spacegrouplabel="#93"
      self.spacegroupnumber=93
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 94  P4_:22_:12 #94
    if(self.sg_num==94) :
      self.spacegroup="P4_:22_:12"
      self.spacegrouplabel="#94"
      self.spacegroupnumber=94
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 95  P4_:322 #95
    if(self.sg_num==95) :
      self.spacegroup="P4_:322"
      self.spacegrouplabel="#95"
      self.spacegroupnumber=95
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-y,x,z+3/4")
      self.sym_ops.append("y,-x,z+1/4")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y,x,-z+1/4")
      self.sym_ops.append("-y,-x,-z+3/4")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 96  P4_:32_:12 #96
    if(self.sg_num==96) :
      self.spacegroup="P4_:32_:12"
      self.spacegrouplabel="#96"
      self.spacegroupnumber=96
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,z+3/4")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/4")
      self.sym_ops.append("-x+1/2,y+1/2,-z+3/4")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/4")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 97  I422 #97
    if(self.sg_num==97) :
      self.spacegroup="I422"
      self.spacegrouplabel="#97"
      self.spacegroupnumber=97
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
        
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 98  I4_:122 #98
    if(self.sg_num==98) :
      self.spacegroup="I4_:122"
      self.spacegrouplabel="#98"
      self.spacegroupnumber=98
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      self.sym_ops.append("-x+1/2,y,-z+3/4")
      self.sym_ops.append("x,-y+1/2,-z+1/4")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y,-x,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 99  P4mm #99
    if(self.sg_num==99) :
      self.spacegroup="P4mm"
      self.spacegrouplabel="#99"
      self.spacegroupnumber=99
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 100  P4bm #100
    if(self.sg_num==100) :
      self.spacegroup="P4bm"
      self.spacegrouplabel="#100"
      self.spacegroupnumber=100
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 101  P4_:2cm #101
    if(self.sg_num==101) :
      self.spacegroup="P4_:2cm"
      self.spacegrouplabel="#101"
      self.spacegroupnumber=101
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 102  P4_:2nm #102
    if(self.sg_num==102) :
      self.spacegroup="P4_:2nm"
      self.spacegrouplabel="#102"
      self.spacegroupnumber=102
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 103  P4cc #103
    if(self.sg_num==103) :
      self.spacegroup="P4cc"
      self.spacegrouplabel="#103"
      self.spacegroupnumber=103
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 104  P4nc #104
    if(self.sg_num==104) :
      self.spacegroup="P4nc"
      self.spacegrouplabel="#104"
      self.spacegroupnumber=104
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 105  P4_:2mc #105
    if(self.sg_num==105) :
      self.spacegroup="P4_:2mc"
      self.spacegrouplabel="#105"
      self.spacegroupnumber=105
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 106  P4_:2bc #106
    if(self.sg_num==106) :
      self.spacegroup="P4_:2bc"
      self.spacegrouplabel="#106"
      self.spacegroupnumber=106
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 107  I4mm #107
    if(self.sg_num==107) :
      self.spacegroup="I4mm"
      self.spacegrouplabel="#107"
      self.spacegroupnumber=107
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 108  I4cm #108
    if(self.sg_num==108) :
      self.spacegroup="I4cm"
      self.spacegrouplabel="#108"
      self.spacegroupnumber=108
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 109  I4_:1md #109
    if(self.sg_num==109) :
      self.spacegroup="I4_:1md"
      self.spacegrouplabel="#109"
      self.spacegroupnumber=109
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y,-x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,x,z+3/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 110  I4_:1cd #110
    if(self.sg_num==110) :
      self.spacegroup="I4_:1cd"
      self.spacegrouplabel="#110"
      self.spacegroupnumber=110
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y,-x+1/2,z+3/4")
      self.sym_ops.append("y+1/2,x,z+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 111  P-42m #111
    if(self.sg_num==111) :
      self.spacegroup="P-42m"
      self.spacegrouplabel="#111"
      self.spacegroupnumber=111
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 112  P-42c #112
    if(self.sg_num==112) :
      self.spacegroup="P-42c"
      self.spacegrouplabel="#112"
      self.spacegroupnumber=112
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 113  P-42_:1m #113
    if(self.sg_num==113) :
      self.spacegroup="P-42_:1m"
      self.spacegrouplabel="#113"
      self.spacegroupnumber=113
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 114  P-42_:1c #114
    if(self.sg_num==114) :
      self.spacegroup="P-42_:1c"
      self.spacegrouplabel="#114"
      self.spacegroupnumber=114
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 115  P-4m2 #115
    if(self.sg_num==115) :
      self.spacegroup="P-4m2"
      self.spacegrouplabel="#115"
      self.spacegroupnumber=115
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 116  P-4c2 #116
    if(self.sg_num==116) :
      self.spacegroup="P-4c2"
      self.spacegrouplabel="#116"
      self.spacegroupnumber=116
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 117  P-4b2 #117
    if(self.sg_num==117) :
      self.spacegroup="P-4b2"
      self.spacegrouplabel="#117"
      self.spacegroupnumber=117
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 118  P-4n2 #118
    if(self.sg_num==118) :
      self.spacegroup="P-4n2"
      self.spacegrouplabel="#118"
      self.spacegroupnumber=118
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 119  I-4m2 #119
    if(self.sg_num==119) :
      self.spacegroup="I-4m2"
      self.spacegrouplabel="#119"
      self.spacegroupnumber=119
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 120  I-4c2 #120
    if(self.sg_num==120) :
      self.spacegroup="I-4c2"
      self.spacegrouplabel="#120"
      self.spacegroupnumber=120
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 121  I-42m #121
    if(self.sg_num==121) :
      self.spacegroup="I-42m"
      self.spacegrouplabel="#121"
      self.spacegroupnumber=121
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 122  I-42d #122
    if(self.sg_num==122) :
      self.spacegroup="I-42d"
      self.spacegrouplabel="#122"
      self.spacegroupnumber=122
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("-x+1/2,y,-z+3/4")
      self.sym_ops.append("x+1/2,-y,-z+3/4")
      self.sym_ops.append("-y+1/2,-x,z+3/4")
      self.sym_ops.append("y+1/2,x,z+3/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 123  P4/mmm #123
    if(self.sg_num==123) :
      self.spacegroup="P4/mmm"
      self.spacegrouplabel="#123"
      self.spacegroupnumber=123
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 124  P4/mcc #124
    if(self.sg_num==124) :
      self.spacegroup="P4/mcc"
      self.spacegrouplabel="#124"
      self.spacegroupnumber=124
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 125  P4/nbm #125
    if(self.sg_num==125 and self.option==1) :
      self.spacegroup="P4/nbm"
      self.spacegrouplabel="#125"
      self.spacegroupnumber=125
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x+1/2,-z")
      self.sym_ops.append("-y+1/2,x+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    if(self.sg_num==125 and self.option==2) :
      self.spacegroup="P4/nbm"
      self.spacegrouplabel="#125"
      self.spacegroupnumber=125
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z")
      self.sym_ops.append("y,-x+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z")
      self.sym_ops.append("x,-y+1/2,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z")
      self.sym_ops.append("-y,x+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x,y+1/2,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 126  P4/nnc #126
    if(self.sg_num==126 and self.option==1) :
      self.spacegroup="P4/nnc"
      self.spacegrouplabel="#126"
      self.spacegroupnumber=126
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    if(self.sg_num==126 and self.option==2) :
      self.spacegroup="P4/nnc"
      self.spacegrouplabel="#126"
      self.spacegroupnumber=126
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z")
      self.sym_ops.append("y,-x+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z")
      self.sym_ops.append("-y,x+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 127  P4/mbm #127
    if(self.sg_num==127) :
      self.spacegroup="P4/mbm"
      self.spacegrouplabel="#127"
      self.spacegroupnumber=127
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 128  P4/mnc #128
    if(self.sg_num==128) :
      self.spacegroup="P4/mnc"
      self.spacegrouplabel="#128"
      self.spacegroupnumber=128
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 129  P4/nmm #129
    if(self.sg_num==129 and self.option==1) :
      self.spacegroup="P4/nmm"
      self.spacegrouplabel="#129"
      self.spacegroupnumber=129
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z")
      self.sym_ops.append("y+1/2,-x+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    if(self.sg_num==129 and self.option==2) :
      self.spacegroup="P4/nmm"
      self.spacegrouplabel="#129"
      self.spacegroupnumber=129
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z")
      self.sym_ops.append("y,-x+1/2,z")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,-z")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z")
      self.sym_ops.append("-y,x+1/2,-z")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 130  P4/ncc #130
    if(self.sg_num==130 and self.option==1) :
      self.spacegroup="P4/ncc"
      self.spacegrouplabel="#130"
      self.spacegroupnumber=130
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z")
      self.sym_ops.append("y+1/2,-x+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x+1/2,-y+1/2,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    if(self.sg_num==130 and self.option==2) :
      self.spacegroup="P4/ncc"
      self.spacegrouplabel="#130"
      self.spacegroupnumber=130
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z")
      self.sym_ops.append("y,-x+1/2,z")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z")
      self.sym_ops.append("-y,x+1/2,-z")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 131  P4_:2/mmc #131
    if(self.sg_num==131) :
      self.spacegroup="P4_:2/mmc"
      self.spacegrouplabel="#131"
      self.spacegroupnumber=131
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z+1/2")
      self.sym_ops.append("-y,x,-z+1/2")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 132  P4_:2/mcm #132
    if(self.sg_num==132) :
      self.spacegroup="P4_:2/mcm"
      self.spacegrouplabel="#132"
      self.spacegroupnumber=132
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z+1/2")
      self.sym_ops.append("-y,x,-z+1/2")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 133  P4_:2/nbc #133
    if(self.sg_num==133 and self.option==1) :
      self.spacegroup="P4_:2/nbc"
      self.spacegrouplabel="#133"
      self.spacegroupnumber=133
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    if(self.sg_num==133 and self.option==2) :
      self.spacegroup="P4_:2/nbc"
      self.spacegrouplabel="#133"
      self.spacegroupnumber=133
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z+1/2")
      self.sym_ops.append("y,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y,-z")
      self.sym_ops.append("x,-y+1/2,-z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z+1/2")
      self.sym_ops.append("-y,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x,y+1/2,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 134  P4_:2/nnm #134
    if(self.sg_num==134 and self.option==1) :
      self.spacegroup="P4_:2/nnm"
      self.spacegrouplabel="#134"
      self.spacegroupnumber=134
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    if(self.sg_num==134 and self.option==2) :
      self.spacegroup="P4_:2/nnm"
      self.spacegrouplabel="#134"
      self.spacegroupnumber=134
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z+1/2")
      self.sym_ops.append("y,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z+1/2")
      self.sym_ops.append("-y,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 135  P4_:2/mbc #135
    if(self.sg_num==135) :
      self.spacegroup="P4_:2/mbc"
      self.spacegrouplabel="#135"
      self.spacegroupnumber=135
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z+1/2")
      self.sym_ops.append("y,-x,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z+1/2")
      self.sym_ops.append("-y,x,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 136  P4_:2/mnm #136
    if(self.sg_num==136) :
      self.spacegroup="P4_:2/mnm"
      self.spacegrouplabel="#136"
      self.spacegroupnumber=136
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 137  P4_:2/nmc #137
    if(self.sg_num==137 and self.option==1) :
      self.spacegroup="P4_:2/nmc"
      self.spacegrouplabel="#137"
      self.spacegroupnumber=137
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
    
    if(self.sg_num==137 and self.option==2) :
      self.spacegroup="P4_:2/nmc"
      self.spacegrouplabel="#137"
      self.spacegroupnumber=137
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z+1/2")
      self.sym_ops.append("y,-x+1/2,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,-z")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z+1/2")
      self.sym_ops.append("-y,x+1/2,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,z")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 138  P4_:2/ncm #138
    if(self.sg_num==138 and self.option==1) :
      self.spacegroup="P4_:2/ncm"
      self.spacegrouplabel="#138"
      self.spacegroupnumber=138
      self.spacegroupoption="origin choice 1"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z")
    
    if(self.sg_num==138 and self.option==2) :
      self.spacegroup="P4_:2/ncm"
      self.spacegrouplabel="#138"
      self.spacegroupnumber=138
      self.spacegroupoption="origin choice 2"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-y+1/2,x,z+1/2")
      self.sym_ops.append("y,-x+1/2,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z+1/2")
      self.sym_ops.append("-y,x+1/2,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 139  I4/mmm #139
    if(self.sg_num==139) :
      self.spacegroup="I4/mmm"
      self.spacegrouplabel="#139"
      self.spacegroupnumber=139
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 140  I4/mcm #140
    if(self.sg_num==140) :
      self.spacegroup="I4/mcm"
      self.spacegrouplabel="#140"
      self.spacegroupnumber=140
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-x,y,-z+1/2")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,-y,z+1/2")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 141  I4_:1/amd #141
    if(self.sg_num==141 and self.option==1) :
      self.spacegroup="I4_:1/amd"
      self.spacegrouplabel="#141"
      self.spacegroupnumber=141
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      self.sym_ops.append("-x+1/2,y,-z+3/4")
      self.sym_ops.append("x,-y+1/2,-z+1/4")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-y+1/2,-z+1/4")
      self.sym_ops.append("x+1/2,y,-z+3/4")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y+1/2,-x,z+3/4")
      self.sym_ops.append("y,x+1/2,z+1/4")
      
    
    if(self.sg_num==141 and self.option==2) :
      self.spacegroup="I4_:1/amd"
      self.spacegrouplabel="#141"
      self.spacegroupnumber=141
      self.spacegroupoption="origin choice 2"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-y+1/4,x+3/4,z+1/4")
      self.sym_ops.append("y+1/4,-x+1/4,z+3/4")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("y+1/4,x+3/4,-z+1/4")
      self.sym_ops.append("-y+1/4,-x+1/4,-z+3/4")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("y+3/4,-x+1/4,-z+3/4")
      self.sym_ops.append("-y+3/4,x+3/4,-z+1/4")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-y+3/4,-x+1/4,z+3/4")
      self.sym_ops.append("y+3/4,x+3/4,z+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 142  I4_:1/acd #142
    if(self.sg_num==142 and self.option==1) :
      self.spacegroup="I4_:1/acd"
      self.spacegrouplabel="#142"
      self.spacegroupnumber=142
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/4")
      self.sym_ops.append("y+1/2,-x,z+3/4")
      self.sym_ops.append("-x+1/2,y,-z+1/4")
      self.sym_ops.append("x,-y+1/2,-z+3/4")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x,-y+1/2,-z+1/4")
      self.sym_ops.append("x+1/2,y,-z+3/4")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y+1/2,-x,z+1/4")
      self.sym_ops.append("y,x+1/2,z+3/4")
          
    if(self.sg_num==142 and self.option==2) :
      self.spacegroup="I4_:1/acd"
      self.spacegrouplabel="#142"
      self.spacegroupnumber=142
      self.spacegroupoption="origin choice 2"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-y+1/4,x+3/4,z+1/4")
      self.sym_ops.append("y+1/4,-x+1/4,z+3/4")
      self.sym_ops.append("-x+1/2,y,-z")
      self.sym_ops.append("x,-y,-z+1/2")
      self.sym_ops.append("y+1/4,x+3/4,-z+3/4")
      self.sym_ops.append("-y+1/4,-x+1/4,-z+1/4")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("y+3/4,-x+1/4,-z+3/4")
      self.sym_ops.append("-y+3/4,x+3/4,-z+1/4")
      self.sym_ops.append("x+1/2,-y,z")
      self.sym_ops.append("-x,y,z+1/2")
      self.sym_ops.append("-y+3/4,-x+1/4,z+1/4")
      self.sym_ops.append("y+3/4,x+3/4,z+3/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 143  P3 #143
    if(self.sg_num==143) :
      self.spacegroup="P3"
      self.spacegrouplabel="#143"
      self.spacegroupnumber=143
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 144  P3_:1 #144
    if(self.sg_num==144) :
      self.spacegroup="P3_:1"
      self.spacegrouplabel="#144"
      self.spacegroupnumber=144
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 145  P3_:2 #145
    if(self.sg_num==145) :
      self.spacegroup="P3_:2"
      self.spacegrouplabel="#145"
      self.spacegroupnumber=145
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 146  R3 #146
    if(self.sg_num==146 and self.option==1) :
      self.spacegroup="R3"
      self.spacegrouplabel="#146"
      self.spacegroupnumber=146
      self.spacegroupoption="rhombohedral axes"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
    
    if(self.sg_num==146 and self.option==2) :
      self.spacegroup="R3"
      self.spacegrouplabel="#146"
      self.spacegroupnumber=146
      self.spacegroupoption="hexagonal axes"
      	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 147  P-3 #147
    if(self.sg_num==147) :
      self.spacegroup="P-3"
      self.spacegrouplabel="#147"
      self.spacegroupnumber=147
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 148  R-3 #148
    if(self.sg_num==148 and self.option==1) :
      self.spacegroup="R-3"
      self.spacegrouplabel="#148"
      self.spacegroupnumber=148
      self.spacegroupoption="rhombohedral axes"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-y,-z,-x")
    
    if(self.sg_num==148 and self.option==2) :
      self.spacegroup="R-3"
      self.spacegrouplabel="#148"
      self.spacegroupnumber=148
      self.spacegroupoption="hexagonal axes"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 149  P312 #149
    if(self.sg_num==149) :
      self.spacegroup="P312"
      self.spacegrouplabel="#149"
      self.spacegroupnumber=149
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 150  P321 #150
    if(self.sg_num==150) :
      self.spacegroup="P321"
      self.spacegrouplabel="#150"
      self.spacegroupnumber=150
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 151  P3_:112 #151
    if(self.sg_num==151) :
      self.spacegroup="P3_:112"
      self.spacegrouplabel="#151"
      self.spacegroupnumber=151
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
      self.sym_ops.append("-y,-x,-z+2/3")
      self.sym_ops.append("-x+y,y,-z+1/3")
      self.sym_ops.append("x,x-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 152  P3_:121 #152
    if(self.sg_num==152) :
      self.spacegroup="P3_:121"
      self.spacegrouplabel="#152"
      self.spacegroupnumber=152
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z+2/3")
      self.sym_ops.append("-x,-x+y,-z+1/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 153  P3_:212 #153
    if(self.sg_num==153) :
      self.spacegroup="P3_:212"
      self.spacegrouplabel="#153"
      self.spacegroupnumber=153
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
      self.sym_ops.append("-y,-x,-z+1/3")
      self.sym_ops.append("-x+y,y,-z+2/3")
      self.sym_ops.append("x,x-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 154  P3_:221 #154
    if(self.sg_num==154) :
      self.spacegroup="P3_:221"
      self.spacegrouplabel="#154"
      self.spacegroupnumber=154
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z+1/3")
      self.sym_ops.append("-x,-x+y,-z+2/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 155  R32 #155
    if(self.sg_num==155 and self.option==1) :
      self.spacegroup="R32"
      self.spacegrouplabel="#155"
      self.spacegroupnumber=155
      self.spacegroupoption="rhombohedral axes"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-z,-y")
    
    if(self.sg_num==155 and self.option==2) :
      self.spacegroup="R32"
      self.spacegrouplabel="#155"
      self.spacegroupnumber=155
      self.spacegroupoption="hexagonal axes"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 156  P3m1 #156
    if(self.sg_num==156) :
      self.spacegroup="P3m1"
      self.spacegrouplabel="#156"
      self.spacegroupnumber=156
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 157  P31m #157
    if(self.sg_num==157) :
      self.spacegroup="P31m"
      self.spacegrouplabel="#157"
      self.spacegroupnumber=157
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 158  P3c1 #158
    if(self.sg_num==158) :
      self.spacegroup="P3c1"
      self.spacegrouplabel="#158"
      self.spacegroupnumber=158
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 159  P31c #159
    if(self.sg_num==159) :
      self.spacegroup="P31c"
      self.spacegrouplabel="#159"
      self.spacegroupnumber=159
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 160  R3m #160
    if(self.sg_num==160 and self.option==1) :
      self.spacegroup="R3m"
      self.spacegrouplabel="#160"
      self.spacegroupnumber=160
      self.spacegroupoption="rhombohedral axes"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("z,y,x")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x,z,y")
    
    if(self.sg_num==160 and self.option==2) :
      self.spacegroup="R3m"
      self.spacegrouplabel="#160"
      self.spacegroupnumber=160
      self.spacegroupoption="hexagonal axes"
      	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 161  R3c #161
    if(self.sg_num==161 and self.option==1) :
      self.spacegroup="R3c"
      self.spacegrouplabel="#161"
      self.spacegroupnumber=161
      self.spacegroupoption="rhombohedral axes"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
    
    if(self.sg_num==161 and self.option==2) :
      self.spacegroup="R3c"
      self.spacegrouplabel="#161"
      self.spacegroupnumber=161
      self.spacegroupoption="hexagonal axes"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 162  P-31m #162
    if(self.sg_num==162) :
      self.spacegroup="P-31m"
      self.spacegrouplabel="#162"
      self.spacegroupnumber=162
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 163  P-31c #163
    if(self.sg_num==163) :
      self.spacegroup="P-31c"
      self.spacegrouplabel="#163"
      self.spacegroupnumber=163
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x+y,y,-z+1/2")
      self.sym_ops.append("x,x-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 164  P-3m1 #164
    if(self.sg_num==164) :
      self.spacegroup="P-3m1"
      self.spacegrouplabel="#164"
      self.spacegroupnumber=164
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 165  P-3c1 #165
    if(self.sg_num==165) :
      self.spacegroup="P-3c1"
      self.spacegrouplabel="#165"
      self.spacegroupnumber=165
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("x-y,-y,-z+1/2")
      self.sym_ops.append("-x,-x+y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 166  R-3m #166
    if(self.sg_num==166 and self.option==1) :
      self.spacegroup="R-3m"
      self.spacegrouplabel="#166"
      self.spacegroupnumber=166
      self.spacegroupoption="rhombohedral axes"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("z,y,x")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x,z,y")
    
    if(self.sg_num==166 and self.option==2) :
      self.spacegroup="R-3m"
      self.spacegrouplabel="#166"
      self.spacegroupnumber=166
      self.spacegroupoption="hexagonal axes"
     	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 167  R-3c #167
    if(self.sg_num==167 and self.option==1) :
      self.spacegroup="R-3c"
      self.spacegrouplabel="#167"
      self.spacegroupnumber=167
      self.spacegroupoption="rhombohedral axes"
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
    
    if(self.sg_num==167 and self.option==2) :
      self.spacegroup="R-3c"
      self.spacegrouplabel="#167"
      self.spacegroupnumber=167
      self.spacegroupoption="hexagonal axes"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("x-y,-y,-z+1/2")
      self.sym_ops.append("-x,-x+y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 168  P6 #168
    if(self.sg_num==168) :
      self.spacegroup="P6"
      self.spacegrouplabel="#168"
      self.spacegroupnumber=168
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 169  P6_:1 #169
    if(self.sg_num==169) :
      self.spacegroup="P6_:1"
      self.spacegrouplabel="#169"
      self.spacegroupnumber=169
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+5/6")
      self.sym_ops.append("x-y,x,z+1/6")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 170  P6_:5 #170
    if(self.sg_num==170) :
      self.spacegroup="P6_:5"
      self.spacegrouplabel="#170"
      self.spacegroupnumber=170
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/6")
      self.sym_ops.append("x-y,x,z+5/6")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 171  P6_:2 #171
    if(self.sg_num==171) :
      self.spacegroup="P6_:2"
      self.spacegrouplabel="#171"
      self.spacegroupnumber=171
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z+2/3")
      self.sym_ops.append("x-y,x,z+1/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 172  P6_:4 #172
    if(self.sg_num==172) :
      self.spacegroup="P6_:4"
      self.spacegrouplabel="#172"
      self.spacegroupnumber=172
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z+1/3")
      self.sym_ops.append("x-y,x,z+2/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 173  P6_:3 #173
    if(self.sg_num==173) :
      self.spacegroup="P6_:3"
      self.spacegrouplabel="#173"
      self.spacegroupnumber=173
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 174  P-6 #174
    if(self.sg_num==174) :
      self.spacegroup="P-6"
      self.spacegrouplabel="#174"
      self.spacegroupnumber=174
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("-y,x-y,-z")
      self.sym_ops.append("-x+y,-x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 175  P6/m #175
    if(self.sg_num==175) :
      self.spacegroup="P6/m"
      self.spacegrouplabel="#175"
      self.spacegroupnumber=175
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("-y,x-y,-z")
      self.sym_ops.append("-x+y,-x,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 176  P6_:3/m #176
    if(self.sg_num==176) :
      self.spacegroup="P6_:3/m"
      self.spacegrouplabel="#176"
      self.spacegroupnumber=176
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("-y,x-y,-z+1/2")
      self.sym_ops.append("-x+y,-x,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 177  P622 #177
    if(self.sg_num==177) :
      self.spacegroup="P622"
      self.spacegrouplabel="#177"
      self.spacegroupnumber=177
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 178  P6_:122 #178
    if(self.sg_num==178) :
      self.spacegroup="P6_:122"
      self.spacegrouplabel="#178"
      self.spacegroupnumber=178
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+5/6")
      self.sym_ops.append("x-y,x,z+1/6")
      self.sym_ops.append("y,x,-z+1/3")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z+2/3")
      self.sym_ops.append("-y,-x,-z+5/6")
      self.sym_ops.append("-x+y,y,-z+1/2")
      self.sym_ops.append("x,x-y,-z+1/6")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 179  P6_:522 #179
    if(self.sg_num==179) :
      self.spacegroup="P6_:522"
      self.spacegrouplabel="#179"
      self.spacegroupnumber=179
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/6")
      self.sym_ops.append("x-y,x,z+5/6")
      self.sym_ops.append("y,x,-z+2/3")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z+1/3")
      self.sym_ops.append("-y,-x,-z+1/6")
      self.sym_ops.append("-x+y,y,-z+1/2")
      self.sym_ops.append("x,x-y,-z+5/6")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 180  P6_:222 #180
    if(self.sg_num==180) :
      self.spacegroup="P6_:222"
      self.spacegrouplabel="#180"
      self.spacegroupnumber=180
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+2/3")
      self.sym_ops.append("-x+y,-x,z+1/3")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z+2/3")
      self.sym_ops.append("x-y,x,z+1/3")
      self.sym_ops.append("y,x,-z+2/3")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z+1/3")
      self.sym_ops.append("-y,-x,-z+2/3")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z+1/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 181  P6_:422 #181
    if(self.sg_num==181) :
      self.spacegroup="P6_:422"
      self.spacegrouplabel="#181"
      self.spacegroupnumber=181
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z+1/3")
      self.sym_ops.append("-x+y,-x,z+2/3")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z+1/3")
      self.sym_ops.append("x-y,x,z+2/3")
      self.sym_ops.append("y,x,-z+1/3")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z+2/3")
      self.sym_ops.append("-y,-x,-z+1/3")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z+2/3")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 182  P6_:322 #182
    if(self.sg_num==182) :
      self.spacegroup="P6_:322"
      self.spacegrouplabel="#182"
      self.spacegroupnumber=182
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x+y,y,-z+1/2")
      self.sym_ops.append("x,x-y,-z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 183  P6mm #183
    if(self.sg_num==183) :
      self.spacegroup="P6mm"
      self.spacegrouplabel="#183"
      self.spacegroupnumber=183
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 184  P6cc #184
    if(self.sg_num==184) :
      self.spacegroup="P6cc"
      self.spacegrouplabel="#184"
      self.spacegroupnumber=184
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 185  P6_:3cm #185
    if(self.sg_num==185) :
      self.spacegroup="P6_:3cm"
      self.spacegrouplabel="#185"
      self.spacegroupnumber=185
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 186  P6_:3mc #186
    if(self.sg_num==186) :
      self.spacegroup="P6_:3mc"
      self.spacegrouplabel="#186"
      self.spacegroupnumber=186
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 187  P-6m2 #187
    if(self.sg_num==187) :
      self.spacegroup="P-6m2"
      self.spacegrouplabel="#187"
      self.spacegroupnumber=187
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("-y,x-y,-z")
      self.sym_ops.append("-x+y,-x,-z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 188  P-6c2 #188
    if(self.sg_num==188) :
      self.spacegroup="P-6c2"
      self.spacegrouplabel="#188"
      self.spacegroupnumber=188
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("-y,x-y,-z+1/2")
      self.sym_ops.append("-x+y,-x,-z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 189  P-62m #189
    if(self.sg_num==189) :
      self.spacegroup="P-62m"
      self.spacegrouplabel="#189"
      self.spacegroupnumber=189
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("-y,x-y,-z")
      self.sym_ops.append("-x+y,-x,-z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 190  P-62c #190
    if(self.sg_num==190) :
      self.spacegroup="P-62c"
      self.spacegrouplabel="#190"
      self.spacegroupnumber=190
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("-y,x-y,-z+1/2")
      self.sym_ops.append("-x+y,-x,-z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 191  P6/mmm #191
    if(self.sg_num==191) :
      self.spacegroup="P6/mmm"
      self.spacegrouplabel="#191"
      self.spacegroupnumber=191
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("-y,x-y,-z")
      self.sym_ops.append("-x+y,-x,-z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 192  P6/mcc #192
    if(self.sg_num==192) :
      self.spacegroup="P6/mcc"
      self.spacegrouplabel="#192"
      self.spacegroupnumber=192
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("y,-x+y,z")
      self.sym_ops.append("x-y,x,z")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("x-y,-y,-z+1/2")
      self.sym_ops.append("-x,-x+y,-z+1/2")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x+y,y,-z+1/2")
      self.sym_ops.append("x,x-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("-y,x-y,-z")
      self.sym_ops.append("-x+y,-x,-z")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 193  P6_:3/mcm #193
    if(self.sg_num==193) :
      self.spacegroup="P6_:3/mcm"
      self.spacegrouplabel="#193"
      self.spacegroupnumber=193
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("x-y,-y,-z+1/2")
      self.sym_ops.append("-x,-x+y,-z+1/2")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("-x+y,y,-z")
      self.sym_ops.append("x,x-y,-z")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("-y,x-y,-z+1/2")
      self.sym_ops.append("-x+y,-x,-z+1/2")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("-x+y,y,z+1/2")
      self.sym_ops.append("x,x-y,z+1/2")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("x-y,-y,z")
      self.sym_ops.append("-x,-x+y,z")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 194  P6_:3/mmc #194
    if(self.sg_num==194) :
      self.spacegroup="P6_:3/mmc"
      self.spacegrouplabel="#194"
      self.spacegroupnumber=194
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-y,x-y,z")
      self.sym_ops.append("-x+y,-x,z")
      self.sym_ops.append("-x,-y,z+1/2")
      self.sym_ops.append("y,-x+y,z+1/2")
      self.sym_ops.append("x-y,x,z+1/2")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("x-y,-y,-z")
      self.sym_ops.append("-x,-x+y,-z")
      self.sym_ops.append("-y,-x,-z+1/2")
      self.sym_ops.append("-x+y,y,-z+1/2")
      self.sym_ops.append("x,x-y,-z+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("y,-x+y,-z")
      self.sym_ops.append("x-y,x,-z")
      self.sym_ops.append("x,y,-z+1/2")
      self.sym_ops.append("-y,x-y,-z+1/2")
      self.sym_ops.append("-x+y,-x,-z+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("-x+y,y,z")
      self.sym_ops.append("x,x-y,z")
      self.sym_ops.append("y,x,z+1/2")
      self.sym_ops.append("x-y,-y,z+1/2")
      self.sym_ops.append("-x,-x+y,z+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 195  P23 #195
    if(self.sg_num==195) :
      self.spacegroup="P23"
      self.spacegrouplabel="#195"
      self.spacegroupnumber=195
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 196  F23 #196
    if(self.sg_num==196) :
      self.spacegroup="F23"
      self.spacegrouplabel="#196"
      self.spacegroupnumber=196
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 197  I23 #197
    if(self.sg_num==197) :
      self.spacegroup="I23"
      self.spacegrouplabel="#197"
      self.spacegroupnumber=197
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 198  P2_:13 #198
    if(self.sg_num==198) :
      self.spacegroup="P2_:13"
      self.spacegrouplabel="#198"
      self.spacegroupnumber=198
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 199  I2_:13 #199
    if(self.sg_num==199) :
      self.spacegroup="I2_:13"
      self.spacegrouplabel="#199"
      self.spacegroupnumber=199
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 200  Pm-3 #200
    if(self.sg_num==200) :
      self.spacegroup="Pm-3"
      self.spacegrouplabel="#200"
      self.spacegroupnumber=200
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 201  Pn-3 #201
    if(self.sg_num==201 and self.option==1) :
      self.spacegroup="Pn-3"
      self.spacegrouplabel="#201"
      self.spacegroupnumber=201
      self.spacegroupoption="origin choice 1"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-z+1/2,-x+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,x+1/2,y+1/2")
      self.sym_ops.append("z+1/2,x+1/2,-y+1/2")
      self.sym_ops.append("z+1/2,-x+1/2,y+1/2")
      self.sym_ops.append("-y+1/2,-z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,x+1/2")
      self.sym_ops.append("-y+1/2,z+1/2,x+1/2")
      self.sym_ops.append("y+1/2,z+1/2,-x+1/2")
    
    if(self.sg_num==201 and self.option==2) :
      self.spacegroup="Pn-3"
      self.spacegrouplabel="#201"
      self.spacegroupnumber=201
      self.spacegroupoption="origin choice 2"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-x+1/2,y")
      self.sym_ops.append("-z+1/2,x,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/2,z,-x+1/2")
      self.sym_ops.append("y,-z+1/2,-x+1/2")
      self.sym_ops.append("-y+1/2,-z+1/2,x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x+1/2,y+1/2")
      self.sym_ops.append("z+1/2,x+1/2,-y")
      self.sym_ops.append("z+1/2,-x,y+1/2")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y+1/2,-z,x+1/2")
      self.sym_ops.append("-y,z+1/2,x+1/2")
      self.sym_ops.append("y+1/2,z+1/2,-x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 202  Fm-3 #202
    if(self.sg_num==202) :
      self.spacegroup="Fm-3"
      self.spacegrouplabel="#202"
      self.spacegroupnumber=202
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 203  Fd-3 #203
    if(self.sg_num==203 and self.option==1) :
      self.spacegroup="Fd-3"
      self.spacegrouplabel="#203"
      self.spacegroupnumber=203
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("-x+1/4,-y+1/4,-z+1/4")
      self.sym_ops.append("x+1/4,y+1/4,-z+1/4")
      self.sym_ops.append("x+1/4,-y+1/4,z+1/4")
      self.sym_ops.append("-x+1/4,y+1/4,z+1/4")
      self.sym_ops.append("-z+1/4,-x+1/4,-y+1/4")
      self.sym_ops.append("-z+1/4,x+1/4,y+1/4")
      self.sym_ops.append("z+1/4,x+1/4,-y+1/4")
      self.sym_ops.append("z+1/4,-x+1/4,y+1/4")
      self.sym_ops.append("-y+1/4,-z+1/4,-x+1/4")
      self.sym_ops.append("y+1/4,-z+1/4,x+1/4")
      self.sym_ops.append("-y+1/4,z+1/4,x+1/4")
      self.sym_ops.append("y+1/4,z+1/4,-x+1/4")
      
    
    if(self.sg_num==203 and self.option==2) :
      self.spacegroup="Fd-3"
      self.spacegrouplabel="#203"
      self.spacegroupnumber=203
      self.spacegroupoption="origin choice 2"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+3/4,-y+3/4,z")
      self.sym_ops.append("-x+3/4,y,-z+3/4")
      self.sym_ops.append("x,-y+3/4,-z+3/4")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x+3/4,-y+3/4")
      self.sym_ops.append("-z+3/4,-x+3/4,y")
      self.sym_ops.append("-z+3/4,x,-y+3/4")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+3/4,z,-x+3/4")
      self.sym_ops.append("y,-z+3/4,-x+3/4")
      self.sym_ops.append("-y+3/4,-z+3/4,x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/4,y+1/4,-z")
      self.sym_ops.append("x+1/4,-y,z+1/4")
      self.sym_ops.append("-x,y+1/4,z+1/4")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x+1/4,y+1/4")
      self.sym_ops.append("z+1/4,x+1/4,-y")
      self.sym_ops.append("z+1/4,-x,y+1/4")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y+1/4,-z,x+1/4")
      self.sym_ops.append("-y,z+1/4,x+1/4")
      self.sym_ops.append("y+1/4,z+1/4,-x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 204  Im-3 #204
    if(self.sg_num==204) :
      self.spacegroup="Im-3"
      self.spacegrouplabel="#204"
      self.spacegroupnumber=204
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 205  Pa-3 #205
    if(self.sg_num==205) :
      self.spacegroup="Pa-3"
      self.spacegrouplabel="#205"
      self.spacegroupnumber=205
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z+1/2,x+1/2,y")
      self.sym_ops.append("z+1/2,x,-y+1/2")
      self.sym_ops.append("z,-x+1/2,y+1/2")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z+1/2,x+1/2")
      self.sym_ops.append("-y+1/2,z+1/2,x")
      self.sym_ops.append("y+1/2,z,-x+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 206  Ia-3 #206
    if(self.sg_num==206) :
      self.spacegroup="Ia-3"
      self.spacegrouplabel="#206"
      self.spacegroupnumber=206
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z+1/2,x+1/2,y")
      self.sym_ops.append("z+1/2,x,-y+1/2")
      self.sym_ops.append("z,-x+1/2,y+1/2")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z+1/2,x+1/2")
      self.sym_ops.append("-y+1/2,z+1/2,x")
      self.sym_ops.append("y+1/2,z,-x+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 207  P432 #207
    if(self.sg_num==207) :
      self.spacegroup="P432"
      self.spacegrouplabel="#207"
      self.spacegroupnumber=207
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 208  P4_:232 #208
    if(self.sg_num==208) :
      self.spacegroup="P4_:232"
      self.spacegrouplabel="#208"
      self.spacegroupnumber=208
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 209  F432 #209
    if(self.sg_num==209) :
      self.spacegroup="F432"
      self.spacegrouplabel="#209"
      self.spacegroupnumber=209
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 210  F4_:132 #210
    if(self.sg_num==210) :
      self.spacegroup="F4_:132"
      self.spacegrouplabel="#210"
      self.spacegroupnumber=210
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x,-y+1/2")
      self.sym_ops.append("-z,-x+1/2,y+1/2")
      self.sym_ops.append("-z+1/2,x+1/2,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/2,z+1/2,-x")
      self.sym_ops.append("y+1/2,-z,-x+1/2")
      self.sym_ops.append("-y,-z+1/2,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+3/4")
      self.sym_ops.append("-y+1/4,-x+1/4,-z+1/4")
      self.sym_ops.append("y+1/4,-x+3/4,z+3/4")
      self.sym_ops.append("-y+3/4,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+3/4")
      self.sym_ops.append("-x+3/4,z+3/4,y+1/4")
      self.sym_ops.append("-x+1/4,-z+1/4,-y+1/4")
      self.sym_ops.append("x+1/4,-z+3/4,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+3/4")
      self.sym_ops.append("z+1/4,-y+3/4,x+3/4")
      self.sym_ops.append("-z+3/4,y+3/4,x+1/4")
      self.sym_ops.append("-z+1/4,-y+1/4,-x+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 211  I432 #211
    if(self.sg_num==211) :
      self.spacegroup="I432"
      self.spacegrouplabel="#211"
      self.spacegroupnumber=211
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 212  P4_:332 #212
    if(self.sg_num==212) :
      self.spacegroup="P4_:332"
      self.spacegrouplabel="#212"
      self.spacegroupnumber=212
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("y+1/4,x+3/4,-z+3/4")
      self.sym_ops.append("-y+1/4,-x+1/4,-z+1/4")
      self.sym_ops.append("y+3/4,-x+3/4,z+1/4")
      self.sym_ops.append("-y+3/4,x+1/4,z+3/4")
      self.sym_ops.append("x+1/4,z+3/4,-y+3/4")
      self.sym_ops.append("-x+3/4,z+1/4,y+3/4")
      self.sym_ops.append("-x+1/4,-z+1/4,-y+1/4")
      self.sym_ops.append("x+3/4,-z+3/4,y+1/4")
      self.sym_ops.append("z+1/4,y+3/4,-x+3/4")
      self.sym_ops.append("z+3/4,-y+3/4,x+1/4")
      self.sym_ops.append("-z+3/4,y+1/4,x+3/4")
      self.sym_ops.append("-z+1/4,-y+1/4,-x+1/4")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 213  P4_:132 #213
    if(self.sg_num==213) :
      self.spacegroup="P4_:132"
      self.spacegrouplabel="#213"
      self.spacegroupnumber=213
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+1/4")
      self.sym_ops.append("-y+3/4,-x+3/4,-z+3/4")
      self.sym_ops.append("y+1/4,-x+1/4,z+3/4")
      self.sym_ops.append("-y+1/4,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+1/4")
      self.sym_ops.append("-x+1/4,z+3/4,y+1/4")
      self.sym_ops.append("-x+3/4,-z+3/4,-y+3/4")
      self.sym_ops.append("x+1/4,-z+1/4,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+1/4")
      self.sym_ops.append("z+1/4,-y+1/4,x+3/4")
      self.sym_ops.append("-z+1/4,y+3/4,x+1/4")
      self.sym_ops.append("-z+3/4,-y+3/4,-x+3/4")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 214  I4_:132 #214
    if(self.sg_num==214) :
      self.spacegroup="I4_:132"
      self.spacegrouplabel="#214"
      self.spacegroupnumber=214
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+1/4")
      self.sym_ops.append("-y+3/4,-x+3/4,-z+3/4")
      self.sym_ops.append("y+1/4,-x+1/4,z+3/4")
      self.sym_ops.append("-y+1/4,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+1/4")
      self.sym_ops.append("-x+1/4,z+3/4,y+1/4")
      self.sym_ops.append("-x+3/4,-z+3/4,-y+3/4")
      self.sym_ops.append("x+1/4,-z+1/4,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+1/4")
      self.sym_ops.append("z+1/4,-y+1/4,x+3/4")
      self.sym_ops.append("-z+1/4,y+3/4,x+1/4")
      self.sym_ops.append("-z+3/4,-y+3/4,-x+3/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 215  P-43m #215
    if(self.sg_num==215) :
      self.spacegroup="P-43m"
      self.spacegrouplabel="#215"
      self.spacegroupnumber=215
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("z,y,x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("-z,-y,x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 216  F-43m #216
    if(self.sg_num==216) :
      self.spacegroup="F-43m"
      self.spacegrouplabel="#216"
      self.spacegroupnumber=216
      self.spacegroupoption=""	
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("z,y,x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("-z,-y,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 217  I-43m #217
    if(self.sg_num==217) :
      self.spacegroup="I-43m"
      self.spacegrouplabel="#217"
      self.spacegroupnumber=217
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("z,y,x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("-z,-y,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 218  P-43n #218
    if(self.sg_num==218) :
      self.spacegroup="P-43n"
      self.spacegrouplabel="#218"
      self.spacegroupnumber=218
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,x+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 219  F-43c #219
    if(self.sg_num==219) :
      self.spacegroup="F-43c"
      self.spacegrouplabel="#219"
      self.spacegroupnumber=219
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,x+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 220  I-43d #220
    if(self.sg_num==220) :
      self.spacegroup="I-43d"
      self.spacegrouplabel="#220"
      self.spacegroupnumber=220
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("y+1/4,x+1/4,z+1/4")
      self.sym_ops.append("-y+1/4,-x+3/4,z+3/4")
      self.sym_ops.append("y+3/4,-x+1/4,-z+3/4")
      self.sym_ops.append("-y+3/4,x+3/4,-z+1/4")
      self.sym_ops.append("x+1/4,z+1/4,y+1/4")
      self.sym_ops.append("-x+3/4,z+3/4,-y+1/4")
      self.sym_ops.append("-x+1/4,-z+3/4,y+3/4")
      self.sym_ops.append("x+3/4,-z+1/4,-y+3/4")
      self.sym_ops.append("z+1/4,y+1/4,x+1/4")
      self.sym_ops.append("z+3/4,-y+1/4,-x+3/4")
      self.sym_ops.append("-z+3/4,y+3/4,-x+1/4")
      self.sym_ops.append("-z+1/4,-y+3/4,x+3/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 221  Pm-3m #221
    if(self.sg_num==221) :
      self.spacegroup="Pm-3m"
      self.spacegrouplabel="#221"
      self.spacegroupnumber=221
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-z,-y,x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("z,y,x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 222  Pn-3n #222
    if(self.sg_num==222 and self.option==1) :
      self.spacegroup="Pn-3n"
      self.spacegrouplabel="#222"
      self.spacegroupnumber=222
      self.spacegroupoption="origin choice 1"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-z+1/2,-x+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,x+1/2,y+1/2")
      self.sym_ops.append("z+1/2,x+1/2,-y+1/2")
      self.sym_ops.append("z+1/2,-x+1/2,y+1/2")
      self.sym_ops.append("-y+1/2,-z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,x+1/2")
      self.sym_ops.append("-y+1/2,z+1/2,x+1/2")
      self.sym_ops.append("y+1/2,z+1/2,-x+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
    
    if(self.sg_num==222 and self.option==2) :
      self.spacegroup="Pn-3n"
      self.spacegrouplabel="#222"
      self.spacegroupnumber=222
      self.spacegroupoption="origin choice 2"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-x+1/2,y")
      self.sym_ops.append("-z+1/2,x,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/2,z,-x+1/2")
      self.sym_ops.append("y,-z+1/2,-x+1/2")
      self.sym_ops.append("-y+1/2,-z+1/2,x")
      self.sym_ops.append("y,x,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("y,-x+1/2,z")
      self.sym_ops.append("-y+1/2,x,z")
      self.sym_ops.append("x,z,-y+1/2")
      self.sym_ops.append("-x+1/2,z,y")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x,-z+1/2,y")
      self.sym_ops.append("z,y,-x+1/2")
      self.sym_ops.append("z,-y+1/2,x")
      self.sym_ops.append("-z+1/2,y,x")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x+1/2,y+1/2")
      self.sym_ops.append("z+1/2,x+1/2,-y")
      self.sym_ops.append("z+1/2,-x,y+1/2")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y+1/2,-z,x+1/2")
      self.sym_ops.append("-y,z+1/2,x+1/2")
      self.sym_ops.append("y+1/2,z+1/2,-x")
      self.sym_ops.append("-y,-x,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y,x+1/2,-z")
      self.sym_ops.append("y+1/2,-x,-z")
      self.sym_ops.append("-x,-z,y+1/2")
      self.sym_ops.append("x+1/2,-z,-y")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x,z+1/2,-y")
      self.sym_ops.append("-z,-y,x+1/2")
      self.sym_ops.append("-z,y+1/2,-x")
      self.sym_ops.append("z+1/2,-y,-x")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 223  Pm-3n #223
    if(self.sg_num==223) :
      self.spacegroup="Pm-3n"
      self.spacegrouplabel="#223"
      self.spacegroupnumber=223
      self.spacegroupoption=""
      self.origin=np.array([0,0,0])
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 224  Pn-3m #224
    if(self.sg_num==224 and self.option==1) :
      self.spacegroup="Pn-3m"
      self.spacegrouplabel="#224"
      self.spacegroupnumber=224
      self.spacegroupoption="origin choice 1"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-x+1/2,-y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z+1/2")
      self.sym_ops.append("-z+1/2,-x+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,x+1/2,y+1/2")
      self.sym_ops.append("z+1/2,x+1/2,-y+1/2")
      self.sym_ops.append("z+1/2,-x+1/2,y+1/2")
      self.sym_ops.append("-y+1/2,-z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,x+1/2")
      self.sym_ops.append("-y+1/2,z+1/2,x+1/2")
      self.sym_ops.append("y+1/2,z+1/2,-x+1/2")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-z,-y,x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("z,y,x")
    
    if(self.sg_num==224 and self.option==2) :
      self.spacegroup="Pn-3m"
      self.spacegrouplabel="#224"
      self.spacegroupnumber=224
      self.spacegroupoption="origin choice 2"
      
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y+1/2,z")
      self.sym_ops.append("-x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,-z+1/2")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-x+1/2,y")
      self.sym_ops.append("-z+1/2,x,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/2,z,-x+1/2")
      self.sym_ops.append("y,-z+1/2,-x+1/2")
      self.sym_ops.append("-y+1/2,-z+1/2,x")
      self.sym_ops.append("y+1/2,x+1/2,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y+1/2,-x,z+1/2")
      self.sym_ops.append("-y,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,-y")
      self.sym_ops.append("-x,z+1/2,y+1/2")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x+1/2,-z,y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,-x")
      self.sym_ops.append("z+1/2,-y,x+1/2")
      self.sym_ops.append("-z,y+1/2,x+1/2")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,z+1/2")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x+1/2,y+1/2")
      self.sym_ops.append("z+1/2,x+1/2,-y")
      self.sym_ops.append("z+1/2,-x,y+1/2")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y+1/2,-z,x+1/2")
      self.sym_ops.append("-y,z+1/2,x+1/2")
      self.sym_ops.append("y+1/2,z+1/2,-x")
      self.sym_ops.append("-y+1/2,-x+1/2,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y+1/2,x,-z+1/2")
      self.sym_ops.append("y,-x+1/2,-z+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,y")
      self.sym_ops.append("x,-z+1/2,-y+1/2")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x+1/2,z,-y+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,x")
      self.sym_ops.append("-z+1/2,y,-x+1/2")
      self.sym_ops.append("z,-y+1/2,-x+1/2")
      self.sym_ops.append("z,y,x")
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 225  Fm-3m #225
    if(self.sg_num==225) :
      self.spacegroup="Fm-3m"
      self.spacegrouplabel="#225"
      self.spacegroupnumber=225
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-z,-y,x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("z,y,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 226  Fm-3c #226
    if(self.sg_num==226) :
      self.spacegroup="Fm-3c"
      self.spacegrouplabel="#226"
      self.spacegroupnumber=226
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      self.sym_ops.append("-y+1/2,-x+1/2,z+1/2")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("-x+1/2,-z+1/2,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,z+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-y+1/2,x+1/2")
      self.sym_ops.append("-z+1/2,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 227  Fd-3m #227
    if(self.sg_num==227 and self.option==1) :
      self.spacegroup="Fd-3m"
      self.spacegrouplabel="#227"
      self.spacegroupnumber=227
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x,-y+1/2")
      self.sym_ops.append("-z,-x+1/2,y+1/2")
      self.sym_ops.append("-z+1/2,x+1/2,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/2,z+1/2,-x")
      self.sym_ops.append("y+1/2,-z,-x+1/2")
      self.sym_ops.append("-y,-z+1/2,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+3/4")
      self.sym_ops.append("-y+1/4,-x+1/4,-z+1/4")
      self.sym_ops.append("y+1/4,-x+3/4,z+3/4")
      self.sym_ops.append("-y+3/4,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+3/4")
      self.sym_ops.append("-x+3/4,z+3/4,y+1/4")
      self.sym_ops.append("-x+1/4,-z+1/4,-y+1/4")
      self.sym_ops.append("x+1/4,-z+3/4,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+3/4")
      self.sym_ops.append("z+1/4,-y+3/4,x+3/4")
      self.sym_ops.append("-z+3/4,y+3/4,x+1/4")
      self.sym_ops.append("-z+1/4,-y+1/4,-x+1/4")
      self.sym_ops.append("-x+1/4,-y+1/4,-z+1/4")
      self.sym_ops.append("x+1/4,y+3/4,-z+3/4")
      self.sym_ops.append("x+3/4,-y+3/4,z+1/4")
      self.sym_ops.append("-x+3/4,y+1/4,z+3/4")
      self.sym_ops.append("-z+1/4,-x+1/4,-y+1/4")
      self.sym_ops.append("-z+3/4,x+1/4,y+3/4")
      self.sym_ops.append("z+1/4,x+3/4,-y+3/4")
      self.sym_ops.append("z+3/4,-x+3/4,y+1/4")
      self.sym_ops.append("-y+1/4,-z+1/4,-x+1/4")
      self.sym_ops.append("y+3/4,-z+3/4,x+1/4")
      self.sym_ops.append("-y+3/4,z+1/4,x+3/4")
      self.sym_ops.append("y+1/4,z+3/4,-x+3/4")
      self.sym_ops.append("-y+1/2,-x,z+1/2")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,x+1/2,-z+1/2")
      self.sym_ops.append("y+1/2,-x+1/2,-z")
      self.sym_ops.append("-x+1/2,-z,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/2,-y")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z+1/2,-y+1/2")
      self.sym_ops.append("-z+1/2,-y,x+1/2")
      self.sym_ops.append("-z,y+1/2,-x+1/2")
      self.sym_ops.append("z+1/2,-y+1/2,-x")
      self.sym_ops.append("z,y,x")
      
    
    if(self.sg_num==227 and self.option==2) :
      ##   cerr << "aflow_wyckoff.cpp: spacegroup=" << spacegroup << " option=" << option << endl
      self.spacegroup="Fd-3m"
      self.spacegrouplabel="#227"
      self.spacegroupnumber=227
      self.spacegroupoption="origin choice 2"
		
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+3/4,-y+1/4,z+1/2")
      self.sym_ops.append("-x+1/4,y+1/2,-z+3/4")
      self.sym_ops.append("x+1/2,-y+3/4,-z+1/4")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+3/4,-y+1/4")
      self.sym_ops.append("-z+3/4,-x+1/4,y+1/2")
      self.sym_ops.append("-z+1/4,x+1/2,-y+3/4")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/4,z+1/2,-x+3/4")
      self.sym_ops.append("y+1/2,-z+3/4,-x+1/4")
      self.sym_ops.append("-y+3/4,-z+1/4,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+1/2")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y+1/4,-x+1/2,z+3/4")
      self.sym_ops.append("-y+1/2,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+1/2")
      self.sym_ops.append("-x+1/2,z+3/4,y+1/4")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x+1/4,-z+1/2,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+1/2")
      self.sym_ops.append("z+1/4,-y+1/2,x+3/4")
      self.sym_ops.append("-z+1/2,y+3/4,x+1/4")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/4,y+3/4,-z+1/2")
      self.sym_ops.append("x+3/4,-y+1/2,z+1/4")
      self.sym_ops.append("-x+1/2,y+1/4,z+3/4")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z+1/2,x+1/4,y+3/4")
      self.sym_ops.append("z+1/4,x+3/4,-y+1/2")
      self.sym_ops.append("z+3/4,-x+1/2,y+1/4")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y+3/4,-z+1/2,x+1/4")
      self.sym_ops.append("-y+1/2,z+1/4,x+3/4")
      self.sym_ops.append("y+1/4,z+3/4,-x+1/2")
      self.sym_ops.append("-y+1/4,-x+3/4,z+1/2")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y+3/4,x+1/2,-z+1/4")
      self.sym_ops.append("y+1/2,-x+1/4,-z+3/4")
      self.sym_ops.append("-x+1/4,-z+3/4,y+1/2")
      self.sym_ops.append("x+1/2,-z+1/4,-y+3/4")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x+3/4,z+1/2,-y+1/4")
      self.sym_ops.append("-z+1/4,-y+3/4,x+1/2")
      self.sym_ops.append("-z+3/4,y+1/2,-x+1/4")
      self.sym_ops.append("z+1/2,-y+1/4,-x+3/4")
      self.sym_ops.append("z,y,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 228  Fd-3c #228
    if(self.sg_num==228 and self.option==1) :
      self.spacegroup="Fd-3c"
      self.spacegrouplabel="#228"
      self.spacegroupnumber=228
      self.spacegroupoption="origin choice 1"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,-z")
      self.sym_ops.append("x+1/2,-y,-z+1/2")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x,-y+1/2")
      self.sym_ops.append("-z,-x+1/2,y+1/2")
      self.sym_ops.append("-z+1/2,x+1/2,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+1/2,z+1/2,-x")
      self.sym_ops.append("y+1/2,-z,-x+1/2")
      self.sym_ops.append("-y,-z+1/2,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+3/4")
      self.sym_ops.append("-y+1/4,-x+1/4,-z+1/4")
      self.sym_ops.append("y+1/4,-x+3/4,z+3/4")
      self.sym_ops.append("-y+3/4,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+3/4")
      self.sym_ops.append("-x+3/4,z+3/4,y+1/4")
      self.sym_ops.append("-x+1/4,-z+1/4,-y+1/4")
      self.sym_ops.append("x+1/4,-z+3/4,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+3/4")
      self.sym_ops.append("z+1/4,-y+3/4,x+3/4")
      self.sym_ops.append("-z+3/4,y+3/4,x+1/4")
      self.sym_ops.append("-z+1/4,-y+1/4,-x+1/4")
      self.sym_ops.append("-x+3/4,-y+3/4,-z+3/4")
      self.sym_ops.append("x+3/4,y+1/4,-z+1/4")
      self.sym_ops.append("x+1/4,-y+1/4,z+3/4")
      self.sym_ops.append("-x+1/4,y+3/4,z+1/4")
      self.sym_ops.append("-z+3/4,-x+3/4,-y+3/4")
      self.sym_ops.append("-z+1/4,x+3/4,y+1/4")
      self.sym_ops.append("z+3/4,x+1/4,-y+1/4")
      self.sym_ops.append("z+1/4,-x+1/4,y+3/4")
      self.sym_ops.append("-y+3/4,-z+3/4,-x+3/4")
      self.sym_ops.append("y+1/4,-z+1/4,x+3/4")
      self.sym_ops.append("-y+1/4,z+3/4,x+1/4")
      self.sym_ops.append("y+3/4,z+1/4,-x+1/4")
      self.sym_ops.append("-y,-x+1/2,z")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+1/2,x,-z")
      self.sym_ops.append("y,-x,-z+1/2")
      self.sym_ops.append("-x,-z+1/2,y")
      self.sym_ops.append("x,-z,-y+1/2")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+1/2,z,-y")
      self.sym_ops.append("-z,-y+1/2,x")
      self.sym_ops.append("-z+1/2,y,-x")
      self.sym_ops.append("z,-y,-x+1/2")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      
    
    if(self.sg_num==228 and self.option==2) :
      self.spacegroup="Fd-3c"
      self.spacegrouplabel="#228"
      self.spacegroupnumber=228
      self.spacegroupoption="origin choice 2"
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/4,-y+3/4,z+1/2")
      self.sym_ops.append("-x+3/4,y+1/2,-z+1/4")
      self.sym_ops.append("x+1/2,-y+1/4,-z+3/4")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/4,-y+3/4")
      self.sym_ops.append("-z+1/4,-x+3/4,y+1/2")
      self.sym_ops.append("-z+3/4,x+1/2,-y+1/4")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y+3/4,z+1/2,-x+1/4")
      self.sym_ops.append("y+1/2,-z+1/4,-x+3/4")
      self.sym_ops.append("-y+1/4,-z+3/4,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z")
      self.sym_ops.append("-y+1/2,-x+1/2,-z+1/2")
      self.sym_ops.append("y+1/4,-x,z+3/4")
      self.sym_ops.append("-y,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y")
      self.sym_ops.append("-x,z+3/4,y+1/4")
      self.sym_ops.append("-x+1/2,-z+1/2,-y+1/2")
      self.sym_ops.append("x+1/4,-z,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x")
      self.sym_ops.append("z+1/4,-y,x+3/4")
      self.sym_ops.append("-z,y+3/4,x+1/4")
      self.sym_ops.append("-z+1/2,-y+1/2,-x+1/2")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+3/4,y+1/4,-z+1/2")
      self.sym_ops.append("x+1/4,-y+1/2,z+3/4")
      self.sym_ops.append("-x+1/2,y+3/4,z+1/4")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z+1/2,x+3/4,y+1/4")
      self.sym_ops.append("z+3/4,x+1/4,-y+1/2")
      self.sym_ops.append("z+1/4,-x+1/2,y+3/4")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y+1/4,-z+1/2,x+3/4")
      self.sym_ops.append("-y+1/2,z+3/4,x+1/4")
      self.sym_ops.append("y+3/4,z+1/4,-x+1/2")
      self.sym_ops.append("-y+1/4,-x+3/4,z")
      self.sym_ops.append("y+1/2,x+1/2,z+1/2")
      self.sym_ops.append("-y+3/4,x,-z+1/4")
      self.sym_ops.append("y,-x+1/4,-z+3/4")
      self.sym_ops.append("-x+1/4,-z+3/4,y")
      self.sym_ops.append("x,-z+1/4,-y+3/4")
      self.sym_ops.append("x+1/2,z+1/2,y+1/2")
      self.sym_ops.append("-x+3/4,z,-y+1/4")
      self.sym_ops.append("-z+1/4,-y+3/4,x")
      self.sym_ops.append("-z+3/4,y,-x+1/4")
      self.sym_ops.append("z,-y+1/4,-x+3/4")
      self.sym_ops.append("z+1/2,y+1/2,x+1/2")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 229  Im-3m #229
    if(self.sg_num==229) :
      self.spacegroup="Im-3m"
      self.spacegrouplabel="#229"
      self.spacegroupnumber=229
      self.spacegroupoption=""
	
      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x,-y,z")
      self.sym_ops.append("-x,y,-z")
      self.sym_ops.append("x,-y,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z,-x,-y")
      self.sym_ops.append("-z,-x,y")
      self.sym_ops.append("-z,x,-y")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z,-x")
      self.sym_ops.append("y,-z,-x")
      self.sym_ops.append("-y,-z,x")
      self.sym_ops.append("y,x,-z")
      self.sym_ops.append("-y,-x,-z")
      self.sym_ops.append("y,-x,z")
      self.sym_ops.append("-y,x,z")
      self.sym_ops.append("x,z,-y")
      self.sym_ops.append("-x,z,y")
      self.sym_ops.append("-x,-z,-y")
      self.sym_ops.append("x,-z,y")
      self.sym_ops.append("z,y,-x")
      self.sym_ops.append("z,-y,x")
      self.sym_ops.append("-z,y,x")
      self.sym_ops.append("-z,-y,-x")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x,y,-z")
      self.sym_ops.append("x,-y,z")
      self.sym_ops.append("-x,y,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z,x,y")
      self.sym_ops.append("z,x,-y")
      self.sym_ops.append("z,-x,y")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z,x")
      self.sym_ops.append("-y,z,x")
      self.sym_ops.append("y,z,-x")
      self.sym_ops.append("-y,-x,z")
      self.sym_ops.append("y,x,z")
      self.sym_ops.append("-y,x,-z")
      self.sym_ops.append("y,-x,-z")
      self.sym_ops.append("-x,-z,y")
      self.sym_ops.append("x,-z,-y")
      self.sym_ops.append("x,z,y")
      self.sym_ops.append("-x,z,-y")
      self.sym_ops.append("-z,-y,x")
      self.sym_ops.append("-z,y,-x")
      self.sym_ops.append("z,-y,-x")
      self.sym_ops.append("z,y,x")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## 230  Ia-3d #230
    if(self.sg_num==230) :
      self.spacegroup="Ia-3d"
      self.spacegrouplabel="#230"
      self.spacegroupnumber=230
      self.spacegroupoption=""

      self.sym_ops.append("x,y,z")
      self.sym_ops.append("-x+1/2,-y,z+1/2")
      self.sym_ops.append("-x,y+1/2,-z+1/2")
      self.sym_ops.append("x+1/2,-y+1/2,-z")
      self.sym_ops.append("z,x,y")
      self.sym_ops.append("z+1/2,-x+1/2,-y")
      self.sym_ops.append("-z+1/2,-x,y+1/2")
      self.sym_ops.append("-z,x+1/2,-y+1/2")
      self.sym_ops.append("y,z,x")
      self.sym_ops.append("-y,z+1/2,-x+1/2")
      self.sym_ops.append("y+1/2,-z+1/2,-x")
      self.sym_ops.append("-y+1/2,-z,x+1/2")
      self.sym_ops.append("y+3/4,x+1/4,-z+1/4")
      self.sym_ops.append("-y+3/4,-x+3/4,-z+3/4")
      self.sym_ops.append("y+1/4,-x+1/4,z+3/4")
      self.sym_ops.append("-y+1/4,x+3/4,z+1/4")
      self.sym_ops.append("x+3/4,z+1/4,-y+1/4")
      self.sym_ops.append("-x+1/4,z+3/4,y+1/4")
      self.sym_ops.append("-x+3/4,-z+3/4,-y+3/4")
      self.sym_ops.append("x+1/4,-z+1/4,y+3/4")
      self.sym_ops.append("z+3/4,y+1/4,-x+1/4")
      self.sym_ops.append("z+1/4,-y+1/4,x+3/4")
      self.sym_ops.append("-z+1/4,y+3/4,x+1/4")
      self.sym_ops.append("-z+3/4,-y+3/4,-x+3/4")
      self.sym_ops.append("-x,-y,-z")
      self.sym_ops.append("x+1/2,y,-z+1/2")
      self.sym_ops.append("x,-y+1/2,z+1/2")
      self.sym_ops.append("-x+1/2,y+1/2,z")
      self.sym_ops.append("-z,-x,-y")
      self.sym_ops.append("-z+1/2,x+1/2,y")
      self.sym_ops.append("z+1/2,x,-y+1/2")
      self.sym_ops.append("z,-x+1/2,y+1/2")
      self.sym_ops.append("-y,-z,-x")
      self.sym_ops.append("y,-z+1/2,x+1/2")
      self.sym_ops.append("-y+1/2,z+1/2,x")
      self.sym_ops.append("y+1/2,z,-x+1/2")
      self.sym_ops.append("-y+1/4,-x+3/4,z+3/4")
      self.sym_ops.append("y+1/4,x+1/4,z+1/4")
      self.sym_ops.append("-y+3/4,x+3/4,-z+1/4")
      self.sym_ops.append("y+3/4,-x+1/4,-z+3/4")
      self.sym_ops.append("-x+1/4,-z+3/4,y+3/4")
      self.sym_ops.append("x+3/4,-z+1/4,-y+3/4")
      self.sym_ops.append("x+1/4,z+1/4,y+1/4")
      self.sym_ops.append("-x+3/4,z+3/4,-y+1/4")
      self.sym_ops.append("-z+1/4,-y+3/4,x+3/4")
      self.sym_ops.append("-z+3/4,y+3/4,-x+1/4")
      self.sym_ops.append("z+3/4,-y+1/4,-x+3/4")
      self.sym_ops.append("z+1/4,y+1/4,x+1/4")
      
    
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
    ## ----------------------------------------------------------------------------
  
  ## cycl atoms
  ## cerr << "DONE atoms=" << self.atoms.size() << endl
  ##DONE


##----------------------------------------------------------------------------

#endif##_WYCKOFF_IMPLEMENTATIONS_

## ***************************************************************************
## *           *
## *       Aflow STEFANO CURTAROLO - Duke University 2003-2008       *
## *           *
## ***************************************************************************


