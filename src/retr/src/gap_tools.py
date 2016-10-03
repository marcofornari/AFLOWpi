import AFLOWpi
import os
import numpy

def gap_size(oneCalc,ID):
    try:
        bandgap_type=gap_type(oneCalc,ID)

        if bandgap_type not in ['p-type','n-type','insulator']:
            print 'conductor..no gap'
            return 0.0
        elif bandgap_type=='p-type':
            en,dos=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[0.1,40.0])
            start=0.1
            end=0.1
            for i in range(len(dos))[:-2]:
#            for i in range(len(dos)):
                if dos[i]>0.00001:
                    end=en[i]
                    break
            return end-start


        elif bandgap_type=='n-type':
            en,dos=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[-40.0,-0.1])
            start=-0.1
            end=0.1
            for i in reversed(range(len(dos))[:-2]):
#            for i in reversed(range(len(dos))):
                if dos[i]>0.00001:
                    end=en[i]
                    break
            return numpy.abs(end-start)


        elif bandgap_type=='insulator':
            en,dos_down=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[-40.0,-0.1])
            en,dos_up=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[0.1,40.0])
            start=-0.1
            end=0.1
            for i in reversed(range(len(dos_down))[:-2]):
                if dos_down[i]>0.00001:
                    end=en[i]
                    break
            for i in range(len(dos_up)[-2]):
                if dos_down[i]>0.00001:
                    start=en[i]
                    break

            return numpy.abs(end-start)



                       
    except Exception,e:
        print e
        return 0.0

def gap_type(oneCalc,ID):
    try:
        dos_range = 0.1
        range_up=[0.05,dos_range]
        range_down=[-1.0*dos_range,-0.05]
        en_above,dos_above=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=range_up)
        en_below,dos_below=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=range_down)
        dos_above_found=False
        dos_below_found=False

        for i in reversed(range(len(dos_above))[2:]):

            if dos_above[i]>0.0001:
                dos_above_found=True
        
        for i in reversed(range(len(dos_below))[2:]):

            if dos_below[i]>0.0001:
                dos_below_found=True

        if dos_above_found==True and dos_below_found==True:
            return 'conductor'
        elif dos_above_found==False and dos_below_found==True:
            return 'p-type'
        elif dos_above_found==True and dos_below_found==False:
            return 'n-type'
        elif dos_above_found==False and dos_below_found==False:
            return 'insulator'
    except:
        return 'None'


def _get_dos(oneCalc,ID,LSDA=False,dos_range=[-0.1,0,1],normalize=True):
    try:
        dos_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='dos',last=True)

	'''extracts HOMO from nscf calculation output file as input to the plotting'''

	try:
		Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
		if type(Efermi)!=type(0.5):
			LSDA=True
			
	except:
		Efermi=0.0
	subdir=oneCalc['_AFLOWPI_FOLDER_']

	"""get the path to the subdirectory of the calc that you are making plots for"""


	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	fileplot = os.path.join(subdir,'DOS_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),oneCalc['_AFLOWPI_PREFIX_']))

	#to set figure size and default fonts


	"""get the path to the subdirectory of the calc that you are making plots for"""
        filedos = os.path.join(subdir,'%s_dos.dat'%dos_ID)
        
	try:
		data = open(filedos,'r').readlines()
	except Exception:
            pass

	en = []
	enup = []
	endown = []
	dos = []
	dosdw = []
        scaling_dosup=[]
        scaling_downdw=[]
        scaling_dos=[]
	for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
		try:
			if LSDA==True:
                            val_up   = float(data[i].split()[0])-Efermi[0]
                            val_down = float(data[i].split()[0])-Efermi[1]
 #                           if val_up < dos_range[1] and valZ_up > dos_range[0]:
                            enup.append(val_up)
#                            if val_down <dos_range[1] and val_down > dos_range[0]:
                            endown.append(val_down)
				
			else:
                            val=float(data[i].split()[0])-Efermi

                        try:
                            scaling_dosdw.append(-1*float(data[i].split()[2]))
                            scaling_dosup.append(float(data[i].split()[1]))

                            if val_down <dos_range[1] and val_down > dos_range[0]:
                                dosdw.append(-1*float(data[i].split()[2]))
                                en_down.append(val_down) #to shift all the y values with respect to the Fermi level
                            if val_up <dos_range[1] and val_up > dos_range[0]:
                                dos.append(float(data[i].split()[1]))
                                en_up.append(val_up) #to shift all the y values with respect to the Fermi level
                        except:
                            scaling_dos.append(float(data[i].split()[1]))
                            if val > dos_range[0] and val <dos_range[1]:
                                dos.append(float(data[i].split()[1]))
                                en.append(val) #to shift all the y values with respect to the Fermi level
		
		except Exception, e:
			pass
	if LSDA==True:
		enup  = map(float,enup)
		endown= map(float,endown)

                floatdosDOWN=map(float,dosdw)
                floatdos=map(float,dos)

#                renormalize_dw=1.0/sum(scaling_dosdw)
#                renormalize_up=1.0/sum(scaling_dosup)
                renormalize_dw=1.0/sum(dosdw)
                renormalize_up=1.0/sum(dosup)

                array_dosup = numpy.asarray(floatdos)*renormalize_up
                floatdos=array_dos.tolist()

                array_dosdw = numpy.asarray(floatdosDOWN)*renormalize_dw
                floatdosDOWN=array_dosdw.tolist()

	else:
		endos=map(float,en)  #to convert the list x from float to numbers
#                renormalize=1.0/sum(scaling_dos)
                renormalize=1.0/sum(dos)                
                floatdos=map(float,dos)
                array_dos = numpy.asarray(floatdos)*renormalize
                floatdos=array_dos.tolist()

	enshift = numpy.array(endos) #to treat the list b as an array?

        
        if LSDA==True:
            return enshift,floatdos,floatdosDOWN
        else:
            return enshift,floatdos

    except:
        return None
