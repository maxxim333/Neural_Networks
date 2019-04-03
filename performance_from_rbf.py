
import pandas as pd
import math
import os
import glob
import numpy as np
from pandas import DataFrame  
import statistics 

def obtain_rbfDB(rbfPath,trainingLOOCV=False):

	rbfDB=[]
	with open(rbfPath) as fh:
		for i in range(5): next(fh)
		for line in fh:
			#Read line
			line=line.rstrip('\n').split()
			error=''
			if len(line)==5: instance,actual,predicted,distribution,attributes=line
			elif len(line)==6: instance,actual,predicted,error,distribution,attributes=line
			else: continue

			#Clean variant,uniprot
			variant,uniprot=attributes.replace('(','').replace(')','').split(',')

			#Clean training
			if actual=='1:?' or actual=='2:?': actual='?'
			elif actual=='1:1': actual='P'
			elif actual=='2:0': actual='N'
			else: sys.exit('ERROR: Which training rbf is this '+actual)

			#Clean label
			if predicted=='1:1': predicted='Pathogenic'
			elif predicted=='2:0': predicted='Neutral'
			else: sys.exit('ERROR: Which predicted rbf is this '+predicted)
			
			#Clean score
			distribution=distribution.split(',')[0].replace('*','')

			rbfDB.append([variant,uniprot,actual,predicted,distribution])

			
	if trainingLOOCV: columns=['variant','protein_uniprot','prediction_training','prediction_label_loocv','prediction_score_loocv']
	else: columns=['variant','protein_uniprot','prediction_training','prediction_label','prediction_score']
	rbfDB=pd.DataFrame(rbfDB,columns=columns)

	return rbfDB

def calculate_contingency_value(row):
		if row.prediction_training=='P' and row.prediction_label=='Pathogenic': return 'TP'
		elif row.prediction_training=='N' and row.prediction_label=='Pathogenic': return 'FP'
		elif row.prediction_training=='P' and row.prediction_label=='Neutral': return 'FN'
		elif row.prediction_training=='N' and row.prediction_label=='Neutral': return 'TN'
		else: return ''

allfiles= os.listdir('weka_output')

average_array_homologs_siftpoly_ohl =[]			
average_array_homologs_siftpoly_2hl=[]
average_array_homologs_siftpoly_params_ohl=[]
average_array_homologs_siftpoly_params_2hl=[]
average_array_orthologs_siftpoly_params_ohl=[]
average_array_orthologs_siftpoly_params_ohl=[]




for file in allfiles:
	if file.endswith(".rbf"):
		print ""
		print "working on file" + file
		myrbf=obtain_rbfDB("weka_output/{0}".format(file),trainingLOOCV=False)
		
		
		tp=0
		fp=0
		tn=0
		fn=0
		for index, row in myrbf.iterrows():
			if calculate_contingency_value(row) == "TP": tp+=1
			elif calculate_contingency_value(row) == "FP": fp+=1
			elif calculate_contingency_value(row) == "TN": tn+=1
			elif calculate_contingency_value(row) =="FN":fn+=1
		
		
		print tn, fp, tp, fn
		
		performanceValues= tp,tn,fp,fn
		
		
		#Metrics
		p=tp+fn
		n=tn+fp
		tot=p+n
		
		#Sensitivity
		sen=float(tp)/p
		
		#Specificity
		denominator=n
		if denominator!=0: spc=float(tn)/denominator
		else: spc=0
		
		#Accuracy
		denominator=tot
		if denominator!=0: acc=(float(tn)+tp)/denominator
		else: acc=0
		
		#MCC
		denominator=math.sqrt((tn+fp)*(tn+fn)*(tp+fp)*(tp+fn))
		if denominator!=0: mcc=(tp*tn-fp*fn)/denominator
		else: mcc=0
		
		#PerformanceDB
		columns=['MCC','SN','SP','ACC','TP','TN','FP','FN','P','N','TOTAL']
		parameters= "MCC: " + str(round(mcc,3))  + " \nSensitivity: " + str(round(sen,3)) +  "\nSpecificity: " + str(round(spc,3)) + "\nACC: " + str(round(acc,3)) + "\nTrue Negative / False Positive /False Negative / Positive / Negative / Total: " + str(tn) + " "+str(fp) +" "+ str(fn)+" " + str(p)+" " + str(n)+" " + str(tot)
		print parameters
		
		if str("_homologs_siftpoly_ohl_") in str(file):
			average_array_homologs_siftpoly_ohl.extend(("SIFT/cadd 0 HL", float(mcc), sen, spc, acc, tn, fp, tp, fn))
		if str("_homologs_siftpoly_2hl_") in str(file):
			average_array_homologs_siftpoly_2hl.extend(("SIFT/cadd 2 HL", mcc, sen, spc, acc, tn, fn, tp, fp, fn))	
		if str("_homologs_siftpoly_pssmnatentropy_ohl_") in str(file):
			average_array_homologs_siftpoly_params_ohl.extend(("Homologs SIFT/cadd + params 0 HL", mcc, sen, spc, acc, tn, fp, tp, fn))
		if str("_homologs_siftpoly_pssmnatentropy_2hl_") in str(file):
			average_array_homologs_siftpoly_params_2hl.extend(("Homologs SIFT/cadd + params 2 HL", mcc, sen, spc, acc, tn, fp, tp, fn))
		if str("_orthologs_siftpoly_pssmnatentropy_ohl_") in str(file):
			average_array_homologs_siftpoly_params_ohl.extend(("Orthologs SIFT/cadd + params 0 HL", mcc, sen, spc, acc, tn, fp, tp, fn))
		if str("_orthologs_siftpoly_pssmnatentropy_2hl_") in str(file):
			average_array_homologs_siftpoly_params_2hl.extend(("Orthologs SIFT/cadd + params 2 HL", mcc, sen, spc, acc, tn, fp, tp, fn))
											
