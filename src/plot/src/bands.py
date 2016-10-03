import AFLOWpi
import numpy
import matplotlib
matplotlib.use('PDF')
from matplotlib import pylab
from matplotlib import pyplot
import os
import logging
import StringIO
import glob
import re
import cPickle


def bands(calcs,yLim=[-10,10],DOSPlot='',runlocal=False,postfix='',tight_banding=False):
	'''
	Generates electronic band structure plots for the calculations in the dictionary of dictionaries
	of calculations with the option to have a DOS or PDOS plot to accompany it.

	Arguments:
              calcs (dict): dictionary of dictionaries representing the set of calculations
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              DOSPlot (str): DOS or the PDOS plots next to the eletronic band structure plot (assuming you ran 
                             either ppDOS or ppPDOS) (default: NONE) 
              postfix (str): Postfix to the filename of the plot
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x

        Returns:
              None

	'''

	for oneCalc in calcs.items():
		if runlocal:
			__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding)
		else:
			AFLOWpi.prep._addToBlock(oneCalc[1],oneCalc[0],'PLOT',"AFLOWpi.plot.__bands(oneCalc,ID,yLim=[%s,%s],DOSPlot='%s',postfix='%s',tight_banding=%s)" % (yLim[0],yLim[1],DOSPlot,postfix,tight_banding))	


def __bands(oneCalc,ID,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False):
	'''
        Wrapper function for AFLOWpi.plot.__bandPlot 
        OBSOLETE. NEEDS REMOVAL


	Arguments:
              oneCalc (dict): Single calculation that is the value of the dictionary of dictionaries of calculations
              ID (str): ID of the calculation
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              DOSPlot (str): DOS or the PDOS plots next to the eletronic band structure plot (assuming you ran 
                             either ppDOS or ppPDOS) (default: NONE) 
              postfix (str): Postfix to the filename of the plot
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x

        Returns:
              None

	'''

        oneCalc = (ID,oneCalc)
        __bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding)





def __getPath_WanT(oneCalc,ID):
    '''


    Arguments:


    Keyword Arguments:

        
    Returns:


    '''
    try:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/*_WanT_bands.out')[-1]
    except:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/*_WanT_bands_up.out')[-1]

    with open(want_stdout_path,'r') as in_file_obj:
        in_string = in_file_obj.read()

    path = AFLOWpi.retr._getPath(0.01,oneCalc,ID=ID)

    plot_bool=[]
    path_name=[]
    path_split = [x for x in  path.split('\n')[1:] if len(x.strip())]

    for i in path_split:
        path_name.append(i.split()[-1])
        if  int(i.split()[3]):

            plot_bool.append(True)
        else:
            plot_bool.append(False)

    gg = re.findall("  Number of kpts in each segment\n((?:.*:\W+(?:\d*)\n)*)",in_string)
    num = [int(x) for x in re.findall('line.*:\W+(\d+)',gg[0])]

    total = 0
    include=[]

    output_path_string = ''
    for i in range(len(num)):
        total+=num[i]+1
        if plot_bool[i]:
            if i==0:
                output_path_string+='%s %s %s %s ! %s\n' %(0.0,0.0,0.0,num[i],path_name[i])
            else:
                output_path_string+='%s %s %s %s ! %s\n' %(0.0,0.0,0.0,num[i]+1,path_name[i],)
        else:
            output_path_string+='%s %s %s %s ! %s\n' %(0.0,0.0,0.0,0,path_name[i])
        for j in range(num[i]):
            include.append(plot_bool[i])

        include.append(True)

    output_path_string+='%s %s %s %s ! %s' %(0.0,0.0,0.0,0,path_name[-1])+'\n' 

#    return 
    print output_path_string
    return  output_path_string



def _clean_bands_data_qe(filebands,Efermi_shift):
     	x = []
	y = []
	k_x = []
	k_y = []

	'''
        looks at the filebands file and reads in the columns for the energy of the
        band and the value of the k point to make the path. this looks at the output
        of the plot_bands.x script and when it finds a blank line in the data it 
        appends a list of k point values and energy values for each band into a list
        of the bands
	'''
	try:

		with open(filebands,'r') as datafile:
			data =datafile.readlines()

		for line in data:
			if line:
				try:
					x_val = float(line.split()[0])
                                        #to shift all the y values with respect to the Fermi level
					y_val = float(line.split()[1])-Efermi_shift
					x.append(x_val)
					y.append(y_val)
				except Exception,e:
					pass
			if not line.split():

				k_x.append(x)
				k_y.append(y)
				x=[]
				y=[]


	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
		logging.warning("output from bands calculation not found. Are you sure you ran ppBands and it completed properly?")
		print "Are you sure you ran ppBands and it completed properly?"
		return
	'''
        Eliminating the gaps between the paths in band structure plots by looking 
        for gaps in the data in the k point values that are greater than some value
	'''
 	gapThreshold = 2.00

 	for band in range(len(k_x)):
 		for kpoint in range(2,len(k_x[band])):
			if k_x[band][kpoint] - k_x[band][kpoint-1] > (k_x[band][kpoint-1] - k_x[band][kpoint-2])*gapThreshold:
				if k_x[band][kpoint-1] - k_x[band][kpoint-2]!=0:
					difference = k_x[band][kpoint] - k_x[band][kpoint-1]
					higher_vals = k_x[band][kpoint:]
				else:
					difference=0
			
				k_x[band][kpoint:] = [x - difference for x in k_x[band][kpoint:]]
		
				




	return k_x,k_y

def __bandPlot(oneCalc,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False): 
	'''
        Function to take the data files generated by the sumpdos ppBands functions and plots 
        the electronic band structure and the projected density of states with energy shifted 
        relative to the Fermi Energy.


	Arguments:
              oneCalc (tuple): Single calculation dictionary and the ID of the calculation 
                               NEEDS TO BE CHANGED TO oneCalc,ID FORMAT SOON
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              DOSPlot (str): DOS or the PDOS plots next to the eletronic band structure plot (assuming you ran 
                             either ppDOS or ppPDOS) (default: NONE) 
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x
              postfix (str): Postfix to the filename of the plot

        Returns:
              None

	'''
	
	calcCopy = oneCalc
	calcID = oneCalc[0]
	oneCalc = oneCalc[1]
        if tight_banding==True:
            calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='PAO-TB',last=True)
        else:
            calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='bands',last=True)

        LSDA=False
        nspin = int(AFLOWpi.scfuj.chkSpinCalc(oneCalc,calcID))
        if nspin!=1:
            LSDA=True

        
	if DOSPlot != '' and DOSPlot != 'APDOS' and DOSPlot != 'DOS':
		print "Not a valid choice for DOSPlot. Valid options are:'APDOS','DOS'"
		return

	if postfix!='':
		postfix='_'+postfix

	if DOSPlot=='APDOS' or DOSPlot=='DOS':
		AFLOWpi.plot.__sumpdos(oneCalc,calcID,TB=tight_banding)	
	

	if tight_banding:
		bandSym = AFLOWpi.plot.__getPath_WanT(oneCalc,calcID)
	else:
		bandSym = AFLOWpi.retr._getPathFromFile(calcCopy)
	
	if bandSym==None:
		print 'ERRORRRR!'
		return
	try:
                if tight_banding:
                    Efermi=0.0
                else:
                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,calcID)
                    print 'EFERMI BANDS NO TB',Efermi

		if type(Efermi)!=type(0.5):
			Efermi=Efermi[0]

	except Exception,e:
            print e
            Efermi=0.0

	subdir=oneCalc['_AFLOWPI_FOLDER_']
       	fileplot = os.path.join(subdir,'BANDS_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix))

	"""get the path to the subdirectory of the calc that you are making plots for"""

	if tight_banding==True:
		filebands = os.path.join(subdir,'%s_bands_want.dat'%calcID)
		if not os.path.exists(filebands):
			filebands_up = os.path.join(subdir,'%s_bands_want_up.dat'%calcID)
			filebands_dn = os.path.join(subdir,'%s_bands_want_down.dat'%calcID)
                Efermi_shift=0.0

	else:
		filebands = os.path.join(subdir,'%s_bands.xmgr'%calcID)
                Efermi_shift=Efermi



        dos_ID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='bands',last=True)
	filedos = os.path.join(subdir,'%s_dos.dat'%dos_ID)

	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	

       #to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=10)             #to set the font size

        
        width = 20
	height = 14
	pylab.figure(figsize=(width, height))#to adjust the figure size
	
     #to do the gaussian smoothing               
 
     #################################
	if nspin==2 and tight_banding==True:
		k_x_up,k_y_up = AFLOWpi.plot._clean_bands_data_qe(filebands_up,Efermi_shift)
		k_x_dn,k_y_dn = AFLOWpi.plot._clean_bands_data_qe(filebands_dn,Efermi_shift)
		k_x,k_y=k_x_up,k_y_up
	else:
		k_x,k_y = AFLOWpi.plot._clean_bands_data_qe(filebands,Efermi_shift)

	a=k_x[1]   # a set of k point values for one band for axis scaling purposes
       	b=k_y[1]



	if DOSPlot != '':
		ax1=pylab.subplot(121)

	if DOSPlot == '':
		ax1=pylab.subplot(111)	
		print 'Plotting electronic band structure of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))

	'''
        Plot each band (k_x[i]),(k_y[i]) on the band structure plot from list of
        the values of energy and position from their respective list by k points
	'''
	if nspin==2 and tight_banding==True:
		for i in range(len(k_x_up)):
			try:
				pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=0.5,label="$\uparrow$")
				pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=0.5,label="$\downarrow$")
			except:
				pass
		handles, labels = ax1.get_legend_handles_labels()
		ax1.legend(handles[-2:], labels[-2:],numpoints=1)

	else:
		if tight_banding==True:
			for i in range(len(k_x)):
				pylab.plot((k_x[i]),(k_y[i]),'k')			
		else:
			for i in range(len(k_x)):
				pylab.plot((k_x[i]),(k_y[i]),'k')


	pylab.ylabel('E(eV)')
	pylab.xlim(min(k_x[1]),max(k_x[1])) 
	pylab.ylim(yLim[0],yLim[1])    
	pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,2))

	'''
        takes in a list of k points that was used as pw.x input for the 'bands'
        calculation as a string. It puts parts of that string into lists and 
        manipulates them to display the symmetry point boundary lines on the 
        band structure plot
	'''
          
	bandSymSplit =  bandSym.split()
	HSPList = []
	HSPSymList = []
	buf = StringIO.StringIO(bandSym)



	for line in buf:
		splitLine = line.split()
		if len(splitLine)==2: # new style kpoint path input case
			HSPList.append(splitLine[1])
			specialPointName = splitLine[-1].rstrip()
			
 	                #renames gG to greek letter for capital gamma
			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG':
				specialPointName = r"$\Gamma$"

			elif len(specialPointName) != 1:
				specialPointName = "$"+specialPointName[0]+r'_{'+specialPointName[1]+'}$' #if there is a subscript it makes the print out on the plot have the number subscripted 
			else:
				specialPointName = "$"+specialPointName[0]+"$" #formats with internal math renderer so all the labels look the same
			HSPSymList.append(specialPointName)
			
		elif len(splitLine)==6: # old style kpoint path input case with kpoint names
			HSPList.append(splitLine[3])
			specialPointName = splitLine[-1].strip()

			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG': #renames gG to greek letter for capital gamma
				specialPointName = r"$\Gamma$"
			elif len(specialPointName) != 1:
				specialPointName = "$"+specialPointName[0]+r'_{'+specialPointName[1]+'}$' #if there is a subscript it makes the print out on the plot have the number subscripted 
			else:
				specialPointName = "$"+specialPointName[0]+"$"  #formats with internal math renderer so all the labels look the same
			HSPSymList.append(specialPointName)                     

		elif len(splitLine)==5: # old style kpoint path input case without kpoint names
			try:
				
				if HSPSymList[-1]!=splitLine[4]:
					HSPList.append(counter)
					HSPSymList.append(splitLine[4])
					counter=1

				else:
					counter+=1
			except Exception,e:
				print e
				counter=1
				HSPSymList.append(splitLine[4])


	'''
        takes the number of k points between each point in the k point paths and
        figures out if they are separate paths (where the end of one path and the 
        next begin have zero k points between them). it also takes the labels for 
        the k points that were in the bands calculation input and creates a list 
        of them and makes a special label for the path boundary e.g. X|Q. All of
        these symmetry lines and symmetry path symbols are put into their own list
        symIndex: for the symmetry line's index in the data set
        symPrint: for the symbols that represent the special symmetry path k points 
	'''

	symIndex = [0]
	totalX =0
	SymPrint = []

	for i in range(len(HSPList)-1):
             try:
		if i==0: # for the first k point in the first path
			SymPrint.append(HSPSymList[i])
		if int(HSPList[i]) == 0 and i<len(HSPList)-2: # for the end of a path (where the number of k points between one point and another is zero)
			continue
		elif int(HSPList[i+1]) == 0 and i!=len(HSPList)-2: # for the point that begins a new path where the end of the last path (which has zero k points from it to this k point)
			totalX +=(int(HSPList[i])+1)
			symIndex.append(totalX)
			mid = '|'
			pathBetweenString = HSPSymList[i+1]+mid+HSPSymList[i+2]
			SymPrint.append(pathBetweenString)
		elif int(HSPList[i+1]) != 0: # for kpoints that are not at the beginning or end of paths
			SymPrint.append(HSPSymList[i+1])
			totalX +=int(HSPList[i])
			symIndex.append(totalX)
		elif i==len(HSPList)-2: # for the end of the last path
			symIndex.append(totalX+int(HSPList[i]))
			SymPrint.append(HSPSymList[i+1])
                elif int(HSPList[i-1]) == 0 and int(HSPList[i]) == 0 and i!=len(HSPList)-2:
			logging.debug('can not find HSP in __bandPlot. This shouldnt be able to be tripped')
	     except:
		     pass
     #add symmetry lines to the band structure plot


	for sym in symIndex:
            try:
                pylab.axvline(a[sym], color = 'k')
            except Exception,e:
                print sym
                pylab.axvline(a[-1], color = 'k')
		
                pass
     #Print path labels to band structure x-axis
	try:
		bars=[]
		for index in symIndex:
			try:
				bars.append(a[index] )
			except:
				bars.append(a[-1] )
		pylab.xticks(bars,SymPrint)
	except Exception,e:
                print e
		pylab.xticks([a[-1] for index in symIndex],SymPrint)

	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Femi level line
	locs, labels = pylab.xticks()

##########################################################################################################


	if DOSPlot == 'APDOS':
		print 'Plotting electronic band structure and projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure and projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))

		fileplot = os.path.join(subdir,'BANDPDOS_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix))	
		ax2=pylab.subplot(122)

		def getPlotData(sumpdosFile):
			try:
				with open(sumpdosFile,'rb') as dataFile:
					data = cPickle.load(dataFile)

			except Exception,e:
				print e
				with open(sumpdosFile,'r') as dataFile:
					data = dataFile.read()
					data = data.split('\n')
					for entry in range(len(data)):
						data[entry]=data[entry].split()


			en = []
			pdos = []
			ldos = []
			ldosDOWN = []
			for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
				try:
					#to shift all the y values with respect to the Fermi level
					en.append(float(data[i][0])-Efermi) 
					ldos.append(float(data[i][1]))
					pdos.append(float(data[i][2]))
					ldosDOWN.append(-1*float(data[i][2]))
					
				except Exception, e:
					pass

			endos=map(float,en)  #to convert the list x from float to numbers
			floatdos=map(float,ldos)
			floatdosDOWN=map(float,ldosDOWN)
			enshift = numpy.array(endos) #to treat the list b as an array?

			return enshift,floatdos,floatdosDOWN

		maxDOS=0
		minDOS=0
		atomName= []
		pDOSNoPath = []
		atomList=[]

		if os.path.exists(os.path.join(subdir,'%s_dos.dat'%dos_ID)):
			pDOSNoPath.append(os.path.join(subdir,'%s_dos.dat'%dos_ID))


		atomList = list(set(AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])))

		for species in atomList:
			filePath = os.path.join(subdir,'%s_All.sumpdos' % species)
			if os.path.exists(filePath):
				pDOSNoPath.append(filePath)
		color = 'k'
		ax2.set_color_cycle(['r','g','b','c', 'm', 'y', 'k'])
		if LSDA:
			ax2.set_color_cycle(['r','r','g','g','b','b','c','c','m','m', 'y','y','k','k'])
		
		for filepath in pDOSNoPath:
			filename =  filepath.split('/')[-1]

			try:
				species = re.findall(r'(.+?)_.+',filename)[0]
			except:
				if filename=='%s_dos.dat'%dos_ID:
					species='TOTAL'
			
			'''gets the energy and the DOS for the orbital for the atom'''
			enshift, floatdos,floatdosDOWN = getPlotData(filepath)
			"""makes a sum of all the pdos for a total"""
			floatdos=AFLOWpi.plot.__smoothGauss(floatdos)
			floatdosDOWN=AFLOWpi.plot.__smoothGauss(floatdosDOWN)
			enshift=AFLOWpi.plot.__smoothGauss(enshift)


			''' scales DOS to larges value of DOS in the given energy range and finds the largest DOS between the different orbitals'''
			try:
				if max([floatdos[k] for k in range(len(floatdos)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ]) > maxDOS:
					maxDOS = max([floatdos[k] for k in range(len(floatdos)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ])
			
				if min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ]) < minDOS:
					minDOS = min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ])
			except:
                            pass
			if species=='TOTAL':	
				pylab.plot(floatdos,enshift,'-',label=species,color='k')		
			else:
				pylab.plot(floatdos,enshift,'-',label=species)							

			if LSDA:
				if species=='TOTAL':
					pylab.plot(floatdosDOWN,enshift,'-',label=species,color='k')		
				else:
					pylab.plot(floatdosDOWN,enshift,'-',label=species)						

		handles, labels = ax2.get_legend_handles_labels()

		if LSDA:
			ax2.legend(handles[::-2], labels[::-2])
			pylab.xlim(1.1*minDOS,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range 
			pylab.axvline(0.0, color = 'k', linewidth = 1.3) #line separating up and down spin
		else:
			ax2.legend(handles[::-1], labels[::-1])
			pylab.xlim(0,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range

		pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,2))

		ax2.spines['bottom'].set_linewidth(1.5)
		ax2.spines['left'].set_linewidth(1.5)
		ax2.spines['right'].set_linewidth(1.5)
		ax2.spines['top'].set_linewidth(1.5)

		ax2.set_yticklabels([])	     #to hide ticklabels
		ax1.set_position([0.07,0.1,0.67,0.8]) #[left,bottom,width,height]
		ax2.set_position([0.75,0.1,0.20,0.8]) #other useful options for the frame! :D

		ax2.yaxis.set_ticks([])
		ax2.yaxis.set_ticks_position('left')
		pylab.xlabel('Density of States (States/eV)')
		ax2.axes.xaxis.set_label_position('top')
		locs, labels = pylab.xticks()
		
		for item in range(len(labels)):
			if item == len(labels)/2:
				labels[item]='arbitrary units'
			else:
				labels[item]=''
		
		ax2.set_xticklabels(labels)

##########################################################################################################
	     #to plot the DOS
	if DOSPlot == 'DOS':
		print 'Plotting electronic band structure and DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure and DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		fileplot = os.path.join(subdir,'BANDDOS_%s%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),oneCalc['_AFLOWPI_PREFIX_'],postfix))
		ax2=pylab.subplot(122)

		try:
			data = open(filedos,'r').readlines()
		except Exception:
			logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
			print "Are you sure you ran ppDOS and it completed properly?"
			return

		en = []
		dos = []
		dosdw = []

		for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
			try:
				#to shift all the y values with respect to the Fermi level 
				en.append(float(data[i].split()[0])-Efermi) 
				dos.append(float(data[i].split()[1]))
				dosdw.append(-1*float(data[i].split()[2]))

			except Exception, e:
				pass

		endos=map(float,en)  #to convert the list x from float to numbers
		floatdos=map(float,dos)
		floatdosDOWN=map(float,dosdw)
		enshift = numpy.array(endos) #to treat the list b as an array?

		floatdos=AFLOWpi.plot.__smoothGauss(floatdos)
		floatdosDOWN=AFLOWpi.plot.__smoothGauss(floatdosDOWN)
		enshift=AFLOWpi.plot.__smoothGauss(enshift)		

		ax2=pylab.subplot(122)

		dosMAX = max([dos[k] for k in range(len(dos)) if en[k] < yLim[1] and en[k] > yLim[0] ])
		pylab.plot(floatdos,enshift,'k') #to plot the smoothed data
		
		if LSDA:
			dosMIN = min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if en[k] < yLim[1] and en[k] > yLim[0] ])
			pylab.plot(floatdosDOWN,enshift,'k') #to plot the smoothed data
			pylab.xlim(1.1*dosMIN,1.1*dosMAX) # scales DOS to larges value of DOS in the given energy range
			pylab.axvline(0.0, color = 'k', linewidth = 1.3) #line separating up and down spin
		else:
			pylab.xlim(0,1.1*dosMAX) # scales DOS to larges value of DOS in the given energy range

		ax2.spines['bottom'].set_linewidth(1.5)
		ax2.spines['left'].set_linewidth(1.5)
		ax2.spines['right'].set_linewidth(1.5)
		ax2.spines['top'].set_linewidth(1.5)

		locs, labels = pylab.xticks()
		
		for item in range(len(labels)):
			if item == len(labels)/2:
				labels[item]='arbitrary units'
			else:
				labels[item]=''
		
		ax2.set_xticklabels(labels)
		ax1.set_position([0.07,0.1,0.67,0.8]) #[left,bottom,width,height]
		ax2.set_position([0.75,0.1,0.20,0.8])

	     #other useful options for the frame! :D
		ax2.yaxis.set_ticks([])
		ax2.yaxis.set_ticks_position('left')
		pylab.xlabel('Density of States (States/eV)')
		ax2.axes.xaxis.set_label_position('top')
###########################################################################################

	ax1.spines['bottom'].set_linewidth(1.5)
	ax1.spines['left'].set_linewidth(1.5)
	ax1.spines['right'].set_linewidth(1.5)
	ax1.spines['top'].set_linewidth(1.5)
	ax1.set_frame_on(True) #or False 

	pylab.ylim(yLim[0],yLim[1])
	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line

     ##############################
     #to increase the linewidth of the axis, on both subplots
        description='Electronic Band Structure'
        if DOSPlot=='APDOS':
            description+=' and Atom Projected DOS'
        if DOSPlot=='DOS':
            description+=' and DOS'
        
	'''gives the name of the compound in the calculation in the name of the file for the band structure plot'''
	figtitle = ''
        compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)
	figtitle = '%s: %s' % (description,compoundNameLatex) 
	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=14,horizontalalignment='center') #[x,y]

	matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

	try:
		AFLOWpi.retr._moveToSavedir(fileplot)
	except Exception,e:
		pass

        pyplot.cla()
	pyplot.clf()
	pyplot.close()

