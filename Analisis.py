#ipython Analisis.py

import numpy as np
import pyneb as pn 

##########################################################################################################

DataFileDict = {
'O2': {'atom': 'o_ii_atom_FFT04.dat', 'coll': 'o_ii_coll_Kal09.dat'},
'O2r': {'atom': 'o_ii_rec_SSB17-B-opt.hdf5'},
'O3': {'atom': 'o_iii_atom_SZ00-WFD96.dat', 'coll': 'o_iii_coll_SSB14.dat'},
'N2': {'atom': 'n_ii_atom_FFT04.dat', 'coll': 'n_ii_coll_T11.dat'},
'S2': {'atom': 's_ii_atom_IFF05.dat', 'coll': 's_ii_coll_TZ10.dat'},
'S3': {'atom': 's_iii_atom_FFTI06.dat', 'coll': 's_iii_coll_GRHK14.dat'},
'Ne3':{'atom': 'ne_iii_atom_McL11.dat', 'coll': 'ne_iii_coll_McL11.dat'},
'Ar3':{'atom': 'ar_iii_atom_M83-KS86.dat', 'coll': 'ar_iii_coll_GMZ95.dat'},
'Ar4':{'atom': 'ar_iv_atom_MZ82.dat', 'coll': 'ar_iv_coll_RB97.dat'},
'Cl2':{'atom': 'cl_ii_atom_MZ83.dat', 'coll': 'cl_ii_coll_T04.dat'},
'Cl3':{'atom': 'cl_iii_atom_Fal99.dat', 'coll': 'cl_iii_coll_BZ89.dat'},
}

pn.atomicData.setDataFileDict(DataFileDict)

########################## Definimos los iones a usar. 
O3 = pn.Atom('O',3)
N2 = pn.Atom('N',2)
O2 = pn.Atom('O',2)
S2 = pn.Atom('S',2)
Cl3 = pn.Atom('Cl',3)
Ar4 = pn.Atom('Ar',4)
S3 = pn.Atom('S',3)
Ne3 = pn.Atom('Ne',3)
Ar3 = pn.Atom('Ar',3)
He1 = pn.RecAtom('He',1)
He2 = pn.RecAtom('He',2)
H1 = pn.RecAtom('H', 1)

########################## archivos de salida
f1 = open('Ab_ion25_22.txt', 'w')
f2 = open('Ab_total25_22.txt', 'w')
f4 = open('Ab_XO25_22.txt', 'w')
f10 = open('Ab_He25_22.txt', 'w')
f11 = open('Diagn25_22.txt', 'w')

""" En caso de querer anadir mas objetos sin modificar los archivos existentes cambiar 'w' por 'a' y comentar los 4 renglones de abajo """
f1.write(str('#') + '\t' +str('#') + '\t' +str('O+b') + '\t' + str('O+r') +'\t' + str('O++')+'\t'+ ('N+')+'\t' + str('Ne++')+'\t' + str('S+')+'\t'+ str('S++')+'\t' + str('Ar++')+'\t'+ str('Cl++')+'\n')
f2.write(str('#') + '\t' +str('#') + '\t' +str('O/Hb') + '\t' + str('N/Hb')+'\t'+ str('Ne/Hb') +'\t'+ str('S/Hb') +'\t'+ str('Ar/Hb') +'\t'+ str('Cl/Hb')+'\t'	
	+ str('O/Hr') + '\t' + str('N/Hr')+'\t'+ str('Ne/Hr') +'\t'+ str('S/Hr') +'\t'+ str('Ar/Hr') +'\t'+ str('Cl/Hr')+'\n')	
f4.write(str('#') + '\t' +str('#') + '\t' +str('N/Ob') +'\t'+ str('Ne/Ob') +'\t'+ str('S/Ob') +'\t'+ str('Ar/Ob') + '\t' + str('Cl/Ob')
	+ '\t' +str('N/Or') +'\t'+ str('Ne/Or') +'\t'+ str('S/Or') +'\t'+ str('Ar/Or') + '\t' + str('Cl/Or')+'\n')
f10.write(str('#') + '\t' +str('#') +'\t'+ str('4026') +'\t'+ str('4388') +'\t'+ str('4922') + '\t' +str('4471') +'\t'+ str('5876') +'\t'+ str('6678') +'\t'+ str('He+') +'\t'+ str('He++') + '\t'+ str('He/H') +'\n')
f11.write(str('#') + '\t' +str('#') + '\t' +str('neS') + '\t' + str('neO') + '\t' +  str('neCl') + '\t' + str('neAr') + '\t' +  str('ne') + '\t' + str('TeN') + '\t' + str('TeO') + '\n')

########################## lista de espectros
NP_au_list=[]
CCD=[]
f = open('/Users/mrg/PNeSpec/GalacticObjects.dat',"r")
lines = f.readlines()
for row in lines:
	x = row.split()
	NP_au_list.append(str(x[0]))
f.close()

########################## intensidades para cada espectro (CCD o no con su score)
for jj in range(len(NP_au_list)):
	
	data_folder = "/Users/mrg/PNeSpec/GalacticSpectra/"
	obs_data = data_folder + NP_au_list[jj] +'.dat'	
	print(NP_au_list[jj])
	datos_obs=open(obs_data, 'r')
	a=[]
	b=[]
	for line in datos_obs:
		datos=line.split()
		a.append(str(datos[0]))
		b.append((datos[1]))
	
	bb=[]
	for i in range (len(b)): #### intensidades
		if b[i] !='NA':
			bb.append(b[i])
		else:
			bb.append(0)
	for i in range(len(a)): #### lineas de emision
		
		if a[i] == 'O2_3726A':
			i3726 = float(bb[i])
		if a[i] == 'O2_3729A':
			i3729 = float(bb[i])
		if a[i] == 'O2_3727A+':
			i3727b = float(bb[i])
		if a[i] == 'Ne3_3868A':
			i3868 = float(bb[i])
		if a[i] == 'Ne3_3967A':
			i3967 = float(bb[i])
		if a[i] == 'He1_4026A':
			i4026 = float(bb[i])
		if a[i] == 'O3_4363A':
			i4363 = float(bb[i])
		if a[i] == 'He1_4388A':
			i4388 = float(bb[i])
		if a[i] == 'He1_4471A':
			i4471 = float(bb[i])
		if a[i] == 'He2_4686A':
			i4686= float(bb[i])
		if a[i] == 'Ar4_4711A':
			i4711 = float(bb[i])
		if a[i] == 'Ar4_4740A':
			i4740 = float(bb[i])
		if a[i] == 'He1_4922A':
			i4922 = float(bb[i])
		if a[i] == 'O3_4959A':
			i4959 = float(bb[i])
		if a[i] == 'O3_5007A':
			i5007 = float(bb[i])
		if a[i] == 'Cl3_5518A':
			i5518 = float(bb[i])
		if a[i] == 'Cl3_5538A':
			i5538 = float(bb[i])
		if a[i] == 'N2_5755A':
			i5755 = float(bb[i])
		if a[i] == 'He1_5876A':
			i5876 = float(bb[i])
		if a[i] == 'S3_6312A':
			i6312 = float(bb[i])
		if a[i] == 'N2_6548A':
			i6548 = float(bb[i])
		if a[i] == 'N2_6584A':
			i6584 = float(bb[i])
		if a[i] == 'He1_6678A':
			i6678 = float(bb[i])
		if a[i] == 'S2_6716A':
			i6716 = float(bb[i])
		if a[i] == 'S2_6731A':
			i6731 = float(bb[i])
		if a[i] == 'Ar3_7135A':
			i7135 = float(bb[i])
		if a[i] == 'O2_7319A+':
			i7319b = float(bb[i])
		if a[i] == 'O2_7330A+':
			i7330b = float(bb[i])
		if a[i] == 'O2_7325A+':
			i7325b = float(bb[i])

	######################################################################################################
	########################## Calculo ne, Te y O++/O+
	######################################################################################################

	temN=10000.0
	temO=10000.0

	for i in range(4):
		if i6716 !=0.0 and i6731 !=0.0:
			ne_s2 = S2.getTemDen(i6716/i6731, tem=temN, to_eval = 'L(6716)/L(6731)')
		else:
			ne_s2 = -1
	
		if i3726 !=0.0 and i3729 !=0.0:
			ne_o2 = O2.getTemDen(i3726/i3729, tem=temN, to_eval = 'L(3726)/L(3729)')
		else:
			ne_o2 = -1

		if i5518 !=0.0 and i5538 !=0.0:
			ne_cl3 = Cl3.getTemDen(i5518/i5538, tem=temN, to_eval = 'L(5518)/L(5538)')
		else:
			ne_cl3 = -1
			
		if i4711 !=0.0 and i4740 !=0.0:
			ne_ar4 = Ar4.getTemDen(i4711/i4740, tem=temO, to_eval = 'L(4711)/L(4740)')
		else:
			ne_ar4 = -1

		ane = np.array([ne_s2, ne_o2, ne_cl3, ne_ar4])
		lne = np.log10(ane)
		mne = np.nanmean(lne)
		den = 10**mne
													
		if i5755 !=0.0:
			if i6548 !=0.0 and i6584 !=0.0:
				temN = N2.getTemDen(i5755/(i6548+i6584), den=den, to_eval = 'L(5755)/(L(6548)+L(6584))')
			if i6548 !=0.0 and i6584 ==0.0:
				temN = N2.getTemDen(i5755/i6548, den=den, to_eval = 'L(5755)/L(6548)')
			if i6548 ==0.0 and i6584 !=0.0:
				temN = N2.getTemDen(i5755/i6584, den=den, to_eval = 'L(5755)/L(6584)')
			
		if i4363 !=0.0:
			if i5007 !=0.0 and i4959 !=0.0:
				temO = O3.getTemDen(i4363/(i5007+i4959), den=den, to_eval = 'L(4363)/(L(4959)+L(5007))') 
			if i5007 !=0.0 and i4959 ==0.0:
				temO = O3.getTemDen(i4363/i5007, den=den, to_eval = 'L(4363)/L(5007)') 
			if i5007 ==0.0 and i4959 !=0.0:
				temO= O3.getTemDen(i4363/i4959, den=den, to_eval = 'L(4363)/L(4959)') 
				
	######################################################################################################
	########################## Ab ionicas de O++, O+, N+
	######################################################################################################

	if i4959 !=0.0 and i5007 !=0.0:
		abO3 = O3.getIonAbundance(int_ratio=(i4959+i5007), tem=temO, den=den, to_eval='L(4959)+L(5007)')
	if i4959 !=0.0 and i5007 ==0.0:
		abO3 = O3.getIonAbundance(int_ratio=(i4959), tem=temO, den=den, to_eval='L(4959)')
	if i4959 ==0.0 and i5007 !=0.0:
		abO3 = O3.getIonAbundance(int_ratio=(i5007), tem=temO, den=den, to_eval='L(5007)')
	labO3 = 12.0 + np.log10(abO3)
	
	if i3726 !=0.0:
		abO2 = O2.getIonAbundance(int_ratio=(i3726+i3729), tem=temN, den=den, to_eval='L(3726)+L(3729)')
	else:
		abO2 = O2.getIonAbundance(int_ratio=(i3727b), tem=temN, den=den, to_eval='L(3726)+L(3729)')
	labO2 = 12.0 + np.log10(abO2)

	if i7319b > 0.0 and i7330b > 0.0:
		abO2r = O2.getIonAbundance(int_ratio=(i7319b+i7330b), tem=temN, den=den, to_eval='L(7319)+L(7320)+L(7330)+L(7331)')
		labO2r = 12.0 + np.log10(abO2r)
	if i7325b > 0.0:
		abO2r = O2.getIonAbundance(int_ratio=i7325b, tem=temN, den=den, to_eval='L(7319)+L(7320)+L(7330)+L(7331)')
		labO2r = 12.0 + np.log10(abO2r)
	if i7319b == 0.0 and i7330b == 0.0 and i7325b == 0.0:
		abO2r = 0.0
		labO2r = 0.0		

	if i6548 !=0.0 and i6584 !=0.0:
		abN2 = N2.getIonAbundance(int_ratio=(i6548+i6584), tem=temN, den=den, to_eval='L(6548)+L(6584)')
	if i6548 !=0.0 and i6584 ==0.0:
		abN2 = N2.getIonAbundance(int_ratio=(i6548), tem=temN, den=den, to_eval='L(6548)')
	if i6548 ==0.0 and i6584 !=0.0:
		abN2 = N2.getIonAbundance(int_ratio=(i6584), tem=temN, den=den, to_eval='L(6584)')
	labN2 = 12.0 + np.log10(abN2)

	######################################################################################################
	########################## archivo con ne, Te
	######################################################################################################

	if np.isnan(ne_s2): ne_s2=-1.0
	if np.isnan(ne_o2): ne_o2=-1.0
	if np.isnan(ne_cl3): ne_cl3=-1.0
	if np.isnan(ne_ar4): ne_ar4=-1.0

	ne_s2 = round(ne_s2)
	ne_o2 = round(ne_o2)
	ne_cl3 = round(ne_cl3)
	ne_ar4 = round(ne_ar4)
	den = round(den)

	f11.write(str(NP_au_list[jj].split('_')[0]) + '\t'+ str(NP_au_list[jj].split('_')[1]) + '\t' + str(ne_s2) + '\t' + str(ne_o2) + '\t' +  str(ne_cl3) + '\t' + str(ne_ar4) + '\t' +  str(den) + '\t' + str(round(temN)) + '\t' + str(round(temO)) + '\n')

	######################################################################################################
	########################## Abundancias ionicas y totales He
	######################################################################################################

	if i4026 != 0.0:
		abHe1a0 = He1.getIonAbundance(i4026, tem=temO, den=den, wave=4026)
		labHe1a0 = 12.0+np.log10(abHe1a0)
	else:
		abHe1a0 = 0
		labHe1a0 = 0

	if i4388 != 0.0:
		abHe1a1 = He1.getIonAbundance(i4388, tem=temO, den=den, wave=4388)
		labHe1a1 = 12.0+np.log10(abHe1a1)
	else:
		abHe1a1 = 0
		labHe1a1 = 0

	if i4922 != 0.0:
		abHe1a2 = He1.getIonAbundance(i4922, tem=temO, den=den, wave=4922)
		labHe1a2 = 12.0+np.log10(abHe1a2)
	else:
		abHe1a2 = 0
		labHe1a2 = 0

	if i4471 != 0.0:
		abHe1a = He1.getIonAbundance(i4471, tem=temO, den=den, wave=4471)
		labHe1a = 12.0+np.log10(abHe1a)
	else:
		abHe1a = 0
		labHe1a = 0
												
	if i5876 != 0.0:
		abHe1b = He1.getIonAbundance(i5876, tem=temO, den=den, wave=5876)
		labHe1b = 12.0+np.log10(abHe1b)
	else:
		abHe1b = 0
		labHe1b = 0
												
	if i6678 != 0.0:
		abHe1c = He1.getIonAbundance(i6678, tem=temO, den=den, wave=6678)
		labHe1c = 12.0+np.log10(abHe1c)
	else:
		abHe1c = 0
		labHe1c = 0
	
	listHe1 = np.array([labHe1a, labHe1b, labHe1c])
	li1 = listHe1[listHe1 > 0]
	labHe1 = np.mean(li1)
	abHe1 = np.power(10,labHe1 - 12.0)


	if i4686 != 0.0:
		abHe2 = He2.getIonAbundance(i4686, tem=temO, den=den, to_eval='I(4,3)')
		labHe2 = 12.0+np.log10(abHe2)
	else:
		abHe2 = 0.0
		labHe2 = 0.0
		
	abHe = (abHe1 + abHe2)	
	labHe = 12+np.log10(abHe1 + abHe2)	

	######################################################################################################
	########################## Resto de abundancias ionicas
	######################################################################################################

	if i5518 !=0.0 and i5538 !=0.0:
		abCl3 = Cl3.getIonAbundance((i5518+i5538), tem=temN, den=den, to_eval='L(5518)+L(5538)')
	if i5518 !=0.0 and i5538 ==0.0:
		abCl3 =Cl3.getIonAbundance(i5518, tem=temN, den=den, to_eval='L(5518)')
	if i5518 ==0.0 and i5538 !=0.0:
		abCl3 = Cl3.getIonAbundance(i5538, tem=temN, den=den, to_eval='L(5538)')
	if i5518 ==0.0 and i5538 ==0.0:
		abCl3 = 0.0
	if abCl3 > 0: labCl3 = 12.0+np.log10(abCl3)
	if abCl3 == 0: labCl3 = 0


	if i7135 !=0.0:
		abAr3 = Ar3.getIonAbundance(i7135, tem=temO, den=den, to_eval='L(7135)')
	else:								
		abAr3 = 0.0
	if abAr3 > 0: labAr3 = 12.0 + np.log10(abAr3)
	if abAr3 == 0: labAr3 = 0


	if i3868 !=0.0 and i3967 !=0.0:
		abNe3 = Ne3.getIonAbundance((i3868+i3967), tem=temO, den=den, to_eval='L(3869)+L(3968)')
	if i3868 !=0.0 and i3967 ==0.0:
		abNe3 = Ne3.getIonAbundance(i3868, tem=temO, den=den, to_eval='L(3869)')
	if i3868 ==0.0 and i3967 ==0.0:
		abNe3 = 0.0
	if abNe3 > 0: labNe3 = 12.0 + np.log10(abNe3)
	if abNe3 == 0: labNe3 = 0


	if i6716 !=0.0 and i6731 !=0.0:
		abS2 = S2.getIonAbundance((i6716+i6731), tem=temN, den=den, to_eval='L(6717)+L(6731)')
	if i6716 !=0.0 and i6731 ==0.0:
		abS2 = S2.getIonAbundance(i6716, tem=temN, den=den, to_eval='L(6717)')																								
	if i6716 ==0.0 and i6731 !=0.0:
		abS2 = S2.getIonAbundance(i6731, tem=temN, den=den, to_eval='L(6731)')												
	if abS2 > 0: labS2 = 12.0+np.log10(abS2)
	if abS2 == 0: labS2 = 0

	if i6312 !=0.0:
		abS3 = S3.getIonAbundance(i6312, tem=temO, den=den, to_eval='L(6312)')
	else:
		abS3 = 0.0
	if abS3 > 0: labS3 = 12.0 + np.log10(abS3)
	if abS3 == 0: labS3 = 0

	f1.write(str(NP_au_list[jj].split('_')[0]) + '\t'+ str(NP_au_list[jj].split('_')[1]) + '\t'+"{0:.3f}".format(labO2) + '\t'+"{0:.3f}".format(labO2r) +'\t' + "{0:.3f}".format(labO3)+'\t' + "{0:.3f}".format(labN2)+'\t' + "{0:.3f}".format(labNe3)+'\t'+ "{0:.3f}".format(labS2)+'\t'+ "{0:.3f}".format(labS3)+'\t' + "{0:.3f}".format(labAr3)+'\t'+ "{0:.3f}".format(labCl3)+'\n')

	######################################################################################################
	########################## Abundancias totales X/H
	######################################################################################################

	########## Oxigeno ############
	if abHe2 > 0.0:
		v = abHe2/(abHe1+abHe2)
		logICF_oxigeno = ((0.08*v)+0.006*(v**2))/(0.34-(0.27*v))
		ICF_oxigeno = 10**(logICF_oxigeno)
		abO = (abO3+abO2)*ICF_oxigeno
		labO = 12+np.log10(abO)
		if abO2r > 0.0:
			abOr = (abO3+abO2r)*ICF_oxigeno
			labOr = 12+np.log10(abOr)
		else:
			abOr = 0
			labOr = 0
	else:
		v=0
		abO=(abO3+abO2)
		labO=12+np.log10(abO)
		if abO2r > 0.0:
			abOr=(abO3+abO2r)
			labOr=12+np.log10(abOr)
		else:
			abOr = 0
			labOr = 0


	w = abO3/(abO2+abO3)
	if abO2r > 0.0: wr = abO3/(abO2r+abO3)

	########## Cloro ############
	if abCl3 > 0:
		ICF_cloro = (4.1620-4.1622*(w**0.21))**0.75
		abCl = (abCl3*ICF_cloro)*(abO/abO2)
		labCl = 12+np.log10(abCl)
	else:
		abCl = 0
		labCl = 0

	if abCl3 > 0 and abO2r > 0:
		ICF_cloror = (4.1620-4.1622*(wr**0.21))**0.75
		abClr = (abCl3*ICF_cloror)*(abOr/abO2r)
		labClr = 12+np.log10(abClr)
	else:
		abClr = 0
		labClr = 0


	########## Nitrogeno ############
	abN = (abN2/abO2)*abO				
	labN = 12+np.log10(abN)	

	if abO2r > 0.0:
		abNr = (abN2/abO2r)*abOr				
		labNr = 12+np.log10(abNr)		
	else:
		abNr = 0				
		labNr = 0	

	########## Argon ############
	if w <= 0.5:
		ICF_argon = (10**(0.05/(0.06 + w)-0.07)) 
	else:
		ICF_argon = (10**((0.03*w)/(0.4-0.3*w)-0.05)) 

	if abAr3 > 0:
		abAr = abAr3*ICF_argon*abO/(abO2+abO3)
		labAr = 12 + np.log10(abAr)
	else:
		abAr = 0
		labAr = 0

	if abAr3 > 0.0 and abO2r > 0.0:
		if wr <= 0.5:
			ICF_argonr = (10**(0.05/(0.06+wr)-0.07)) 
		else:
			ICF_argonr = (10**((0.03*wr)/(0.4-0.3*wr)-0.05)) 
		abArr = (abAr3)*ICF_argonr*abOr/(abO2r+abO3)
		labArr = 12 + np.log10(abArr)	
	else:
		abArr = 0	
		labArr = 0	

	########## Neon ############
	if abNe3 > 0 and i4686 == 0.0:
		ICF_neon = w + ((0.014/0.01)+2*(0.01**2.7))**3 * (0.7+0.2*w-0.8*(w**2))
		abNe = abNe3*ICF_neon*abO/abO3
		labNe = 12.0 + np.log10(abNe)
	if abNe3 > 0 and i4686 > 0.0:
		if v < 0.015:
			ICF_neon = w + ((0.014/0.015)+2*(0.015**2.7))**3 * (0.7+0.2*w-0.8*(w**2))
			abNe = abNe3*ICF_neon*abO/abO3
		if v >= 0.015:
			ICF_neon = w + ((0.014/v)+2*(v**2.7))**3 * (0.7+0.2*w-0.8*(w**2))
			abNe = abNe3*ICF_neon*abO/abO3
		labNe = 12.0 + np.log10(abNe)
	if abNe3 == 0:
		abNe = 0
		labNe = 0
	
	if abNe3 > 0 and i4686 == 0.0 and abO2r > 0:
		ICF_neonr = wr + ((0.014/0.01)+2*(0.01**2.7))**3 * (0.7+0.2*wr-0.8*(wr**2))
		abNer = abNe3*ICF_neonr*abOr/abO3
		labNer = 12.0 + np.log10(abNer)
	if abNe3 > 0 and i4686 > 0.0 and abO2r > 0:
		if v < 0.015:
			ICF_neonr = wr + ((0.014/0.015)+2*(0.015**2.7))**3 * (0.7+0.2*wr-0.8*(wr**2))
			abNer = abNe3*ICF_neonr*abOr/abO3
		if v >= 0.015:
			ICF_neonr = wr + ((0.014/v)+2*(v**2.7))**3 * (0.7+0.2*wr-0.8*(wr**2))
			abNer = abNe3*ICF_neonr*abOr/abO3
		labNer = 12.0 + np.log10(abNer)
	if abNe3 == 0 or abO2r == 0:
		abNer = 0
		labNer = 0
	
	########## Azufre ############	
	if i6312!=0.0:
		ICF_azufre = 10**((-0.02-0.03*w-2.31*(w**2)+2.19*(w**3))/(0.69+2.09*w-2.69*(w**2)))
		abS = (abS2+abS3)*ICF_azufre*abO/abO2
		labS = 12 + np.log10(abS)
	else:
		abS = 0
		labS = 0
	
	if abO2r > 0.0:
		if i6312!=0.0:
			ICF_azufrer = 10**((-0.02-0.03*wr-2.31*(wr**2)+2.19*(wr**3))/(0.69+2.09*wr-2.69*(wr**2)))
			abSr = (abS2+abS3)*ICF_azufrer*abOr/abO2r
			labSr = 12 + np.log10(abSr)
		else:
			abSr = 0
			labSr = 0
	else:
		abSr = 0
		labSr = 0

	f2.write(str(NP_au_list[jj].split('_')[0]) + '\t'+ str(NP_au_list[jj].split('_')[1]) + '\t'+ "{0:.3f}".format(labO) + '\t' + "{0:.3f}".format(labN)+'\t'+ "{0:.3f}".format(labNe) +'\t'+ "{0:.3f}".format(labS) +'\t'+ "{0:.3f}".format(labAr) +'\t'+ "{0:.3f}".format(labCl)+'\t'+ 
		 "{0:.3f}".format(labOr) +  '\t'+ "{0:.3f}".format(labNr)+'\t'+ "{0:.3f}".format(labNer) +'\t'+ "{0:.3f}".format(labSr) +'\t'+ "{0:.3f}".format(labArr) +'\t'+ "{0:.3f}".format(labClr)+'\n')
	f10.write(str(NP_au_list[jj].split('_')[0]) + '\t'+ str(NP_au_list[jj].split('_')[1]) +'\t'+ "{0:.3f}".format(labHe1a0) +'\t'+ "{0:.3f}".format(labHe1a1) +'\t'+ "{0:.3f}".format(labHe1a2) +'\t'+ "{0:.3f}".format(labHe1a) + '\t' + "{0:.3f}".format(labHe1b) + '\t' + "{0:.3f}".format(labHe1c) + '\t' + "{0:.3f}".format(labHe1) + '\t' + "{0:.3f}".format(labHe2)+ '\t' + "{0:.3f}".format(labHe) +'\n')

	######################################################################################################
	########################## Abundancias totales X/O
	######################################################################################################

	if labN > 0: lNO = labN-labO
	if labN == 0: lNO = -10
	if labCl > 0: lClO = labCl-labO
	if labCl == 0: lClO = -10
	if labS > 0: lSO = labS-labO
	if labS == 0: lSO = -10
	if labAr > 0: lArO = labAr-labO
	if labAr == 0: lArO = -10
	if labNe > 0: lNeO = labNe-labO
	if labNe == 0: lNeO = -10

	if labNr > 0: lNOr = labNr-labOr
	if labNr == 0: lNOr = -10
	if labClr > 0: 	lClOr = labClr-labOr
	if labClr == 0: lClOr = -10
	if labSr > 0: lSOr = labSr-labOr
	if labSr == 0: lSOr = -10
	if labArr > 0: lArOr = labArr-labOr
	if labArr == 0: lArOr = -10
	if labNer > 0: lNeOr = labNer-labOr
	if labNer == 0: lNeOr = -10

	f4.write(str(NP_au_list[jj].split('_')[0]) + '\t'+ str(NP_au_list[jj].split('_')[1]) + '\t'+ "{0:.3f}".format(lNO) +'\t'+ "{0:.3f}".format(lNeO) +'\t'+ "{0:.3f}".format(lSO) +'\t'+ "{0:.3f}".format(lArO)+'\t'+ "{0:.3f}".format(lClO)+
		'\t'+ "{0:.3f}".format(lNOr) +'\t'+ "{0:.3f}".format(lNeOr) +'\t'+ "{0:.3f}".format(lSOr) +'\t'+ "{0:.3f}".format(lArOr)+'\t'+ "{0:.3f}".format(lClOr)+'\n')


f1.close()
f2.close()
f4.close()
f10.close()
f11.close()
