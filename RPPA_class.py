import os, sys, argparse, re
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats
from collections import OrderedDict
import fnmatch
import shutil
import xlsxwriter
from matplotlib.backends.backend_pdf import PdfPages
import Image


class RPPA:

    def __init__(self, slidefile, antibodyfile):
        self.slidefile = slidefile
        self.antibodyfile = antibodyfile
        self.slidedata = pd.read_table(self.slidefile, sep="\t")
        self.antibodydata = pd.read_table(self.antibodyfile, sep="\t",
                                          index_col=0)

    def normalize(self, path, start, end, pmt_settings, prot_data, debug=None):
        print "Start in-house normalization"
        self.path = path
        print "gpr files directory = "+str(os.listdir(self.path))
        self.start = start
        self.end = end

        self.dirs = sorted(os.listdir(self.path))
        master_normdata_dic = {}

        gpr_file_dic = {}
        for folder in self.dirs:
            if folder.startswith(self.start):
                negative_slide = {}
                self.folder_with_path = os.path.join(self.path, folder)
                self.path2 = self.folder_with_path
                self.files = sorted(os.listdir(self.path2))
                for filename in self.files:
                    if filename.endswith(self.end):
                        #print filename.strip()
                        if debug is not None:
                            #print "Printing gpr file"
                            print filename.strip()
                        # cols = filename.strip().split('_')
                        # pmt = cols[2]
                        # if pmt == pmt_int:
                        gpr_file_dic[filename] = os.path.join(self.path2,
                                                              filename)

        self.gpr_file_dic = gpr_file_dic

        for setting in pmt_settings:
            self.pmt_int = setting
            ant_sample_dict = OrderedDict()
            # data_list = []
            ab_file = {}
            batch = {}
            neg_file = {}
            total_prot_file = {}
            total_prot_file_cycle = {}
            for i in range(0, len(self.slidedata.index.values)):
                if self.slidedata.values[i, 5] in gpr_file_dic:
                    # print self.slidedata.values[i, 5]
                    cols = (self.slidedata.values[i, 5]).split('_')
                    pmt = cols[2]
                    # print pmt, self.pmt_int
                    if pmt == self.pmt_int:
                        # print self.slidedata.values[i, 5]
                        ab_file[self.slidedata.values[i, 1]] = self.slidedata.values[i, 5]
                        gpr_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,5]], sep="\t",skiprows=32)
                        gpr_flag = gpr_data["Flags"]
                        neg_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,7]], sep="\t",skiprows=32)
                        neg_flag = neg_data["Flags"]
                        print self.slidedata.values[i, 6]
                        # if ~np.isnan(self.slidedata.values[i, 6]):
                        if ~pd.isnull(self.slidedata.values[i, 6]):
                            total_protein_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,6]], sep="\t",skiprows=32)
                            total_flag = total_protein_data["Flags"]
                            fac = []
                            for item in total_protein_data["Name"]:
                                cols = item.split("_")
                                #print cols[0]
                                fac.append(cols[0])
                            grouped = ((total_protein_data["F532 Median - B532"]).astype(np.float64)).groupby([fac])
                            total_median = grouped.median()
                            transformed = grouped.transform(np.median)
                            transformed[pd.isnull((prot_data.ix[:,5]))] = 1
                            (total_protein_data["F532 Median - B532"])[pd.isnull((prot_data.ix[:,5]))] = 1
                            norm_value = ((gpr_data["F635 Median - B635"]).astype(np.float64)-(neg_data["F635 Median - B635"]).astype(np.float64))*transformed/(total_protein_data["F532 Median - B532"]).astype(np.float64)
                        else:
                            norm_value = ((gpr_data["F635 Median - B635"]).astype(np.float64)-(neg_data["F635 Median - B635"]).astype(np.float64))
                        norm_value[norm_value<=0] = 1
                        norm_value[gpr_flag==-100] = "NA"
                        if ~pd.isnull(self.slidedata.values[i, 6]):
                            norm_value[total_flag==-100] = "NA"
                        norm_value[neg_flag==-100] = "NA"
                        ant_sample_dict[self.slidedata.values[i,1]] = norm_value
                        #data_list.append(ant_sample_dict[self.slidedata.values[i,1]])
                        #print data.values[i,3]
                else:
                    print("No gpr file for "+self.slidedata.values[i, 5])                

            combined_data = pd.DataFrame(ant_sample_dict)
            print combined_data.columns
            sample_list = []
            sample_dic = {}
            j = 0
            for i in gpr_data["Name"].values:
                #print i
                if j == 0:
                    sam = i+".1"
                    j += 1
                elif j == 1:
                    sam = i+".2"
                    j += 1
                elif j == 2:
                    sam = i+".3"
                    j = 0
                sample_list.append(sam)
                sample_dic[i] = i
            
            #combined_data1 = combined_data.set_index(sample_list)
            combined_data1 = combined_data.copy()        
            combined_data1.index = sample_list
            master_data = combined_data1.transpose()
            ab_file_column = pd.DataFrame(pd.Series(ab_file))
            master_data.insert(0, "Slide_file", ab_file_column)
            #master_data.to_csv("RPPA0014_Norm_"+pmt_int+".xls",sep="\t", na_rep="NA")
            master_normdata_dic[self.pmt_int] = master_data
            #print master_data
            
        return master_normdata_dic
        
    def norm_data(self,qpath,qstart,pmt_settings):
        self.qpath = qpath
        #dirs = sorted(os.listdir(path))
        self.qstart = qstart
        #self.qend = qend
        self.dirs = sorted(os.listdir(self.qpath))
        master_normdata_dic = {}
        
        for setting in pmt_settings:            
            self.pmt_int = setting
            
            qend = self.pmt_int+".xls"
            
            i = 0
            for folder in self.dirs:
                if folder.startswith(self.qstart):
                    folder_with_path = os.path.join(self.qpath, folder)
                    path2 = folder_with_path
                    files = sorted(os.listdir(path2))
                    for filename in files:
                        if filename.endswith(qend):
                            file_with_path = os.path.join(path2, filename)
                            if i == 0:
                                df = pd.read_table(file_with_path, sep="\t", skiprows=1)
                                data = df
                            else:
                                df = pd.read_table(file_with_path, sep="\t", skiprows=1)
                                data = pd.merge(data,df)
                                #data = data1                    
                            #print filename        
                            i += 1
            self.normdata = (data.ix[:,7:]).transpose()

            sample_list = []
            sam = ""
            for i in data.ix[:,3].values:
                if sam == i+".1":
                    sam = i+".2"
                elif sam == i+".2":
                    sam = i+".3"
                else:
                    sam = i+".1"
                sample_list.append(sam)
            self.normdata.columns = sample_list
            
            
            slide_ab = {}
            slide_name = {}
            slide_pmt = {}
            for i in range(0,len(self.slidedata.index.values)):
                #print slide_data.values[i,5]
                cols = (self.slidedata.values[i,5]).rstrip().split("_")
                slide_ab[cols[3]] = self.slidedata.values[i,1]
                if cols[2] == self.pmt_int or cols[2] == self.pmt_int:
                    slide_name[slide_ab[cols[3]]] = self.slidedata.values[i,5]
                slide_pmt[slide_ab[cols[3]]] = cols[2]
                #ant_name[slide_ab[cols[3]]] =  self.antibodydata.ix[slide_ab[cols[3]],"PI_name"]
                total_cols = (self.slidedata.values[i,6]).rstrip().split("_")
                slide_ab[total_cols[3]] = "Total_"+total_cols[3]
                neg_cols = (self.slidedata.values[i,7]).rstrip().split("_")
                slide_ab[neg_cols[3]+neg_cols[2]] = "Neg_"+neg_cols[3]
                #self.antibodyfile = antibodyfile"""
                
            slide_column = pd.DataFrame(pd.Series(slide_name))
            #antibody_column = pd.DataFrame(pd.Series(ant_name))
            
            index_list = [] 
            for item in self.normdata.index:
                print item
                cols = item.strip().split("_")
                if (cols[3] or cols[0]) in slide_ab:
                    if item.startswith("G"):
                        index_list.append(slide_ab[cols[3]])
                    else:
                        index_list.append(slide_ab[cols[0]])
                else:
                    #print item
                    self.normdata = self.normdata.drop(item,axis=0)
                    #self.normdata340 = self.normdata340.drop(item,axis=0)
            self.normdata.index = index_list
            self.normdata.insert(0, "Slide_file", slide_column)
            master_normdata_dic[self.pmt_int] = self.normdata
            
        return master_normdata_dic
 
    
    def pi_data(self,masterfile,samplefile,experiment,drop_items,conf_file):
        
        self.masterfile = masterfile
        self.samplefile = samplefile
        
        pi_files_path = experiment+"_results/pi_files/"
        
        if not os.path.exists(pi_files_path+self.masterfile[:-4]):
                os.makedirs(pi_files_path+self.masterfile[:-4])
        
        self.masterdata_all = pd.read_table(self.masterfile, sep="\t", index_col=0)
        self.sampledata_all = pd.read_table(self.samplefile, sep="\t", index_col=0)
        self.masterdata_raw = self.masterdata_all.ix[:,1:]
        self.sampledata_raw = self.sampledata_all.ix[:,1:]
        
        ab_name = {}
        ab_gene = {}
        ab_swiss = {}
        for i in self.masterdata_all.index:
            if str(i) in self.antibodydata.index.astype(str):
                #print i
                if self.antibodydata.ix[int(i),"Host"] == "Mouse":
                    h = "M"
                elif self.antibodydata.ix[int(i),"Host"] == "Rabbit":
                    h = "R"
                elif self.antibodydata.ix[int(i),"Host"] == "Goat":
                    h = "G"
                else:
                    h = "N"
                    
                if self.antibodydata.ix[int(i),"Current_Validation_Status"] == "Validated":
                    s = "V"
                elif self.antibodydata.ix[int(i),"Current_Validation_Status"] == "Not_Valid":
                    s = "N"
                elif self.antibodydata.ix[int(i),"Current_Validation_Status"] == "Caution":
                    s = "C"
                elif self.antibodydata.ix[int(i),"Current_Validation_Status"] == "Progress":
                    s = "P"
                elif self.antibodydata.ix[int(i),"Current_Validation_Status"] == "Failed":
                    s = "F"
                else:
                    s = ""
                ab_name[i] = self.antibodydata.ix[int(i),"PI_name"]+"_"+h+"_"+s    
                ab_gene[i] = self.antibodydata.ix[int(i),"Gene_ID"]
                ab_swiss[i] = self.antibodydata.ix[int(i),"Swiss-Prot_Acc"]
            else:
                print i
                if str(i).startswith("Pr"):
                    ab_name[(i)] = "Tot_"+i
                else:
                    ab_name[(i)] = "Neg_"+i[-7]+"_"+i
                ab_gene[(i)] = "None"
                ab_swiss[(i)] = "None"           
            
        ab_name_column = pd.DataFrame(pd.Series(ab_name))
        ab_gene_column = pd.DataFrame(pd.Series(ab_gene))
        ab_swiss_column = pd.DataFrame(pd.Series(ab_swiss))
            
        
        drop_list = []
        if drop_items is not None:
            for item in (drop_items[0]):
                for i in self.masterdata_raw.columns:
                    #print drop_items[0],item,i
                    if re.match(item, i):
                        drop_list.append(i)
                        print i
        
        self.masterdata_clean = self.masterdata_raw.drop(drop_list, axis=1)
        self.sampledata_clean = self.sampledata_raw.drop(drop_list, axis=1)
            
        drop_max_list = []
        drop_items_for_max = drop_items[1]
        if drop_items_for_max is not None:
            for item in drop_items_for_max:
                for i in self.sampledata_clean.columns:
                    if re.match(item, i):
                        drop_max_list.append(i)
                        print i
        
        self.masterdata_clean_for_max = self.masterdata_raw.drop(drop_list, axis=1)
        self.sampledata_clean_for_max = self.sampledata_raw.drop(drop_list, axis=1)
            
        fac = []
        for i in self.masterdata_clean.columns.values:
            a = i.split('_')
            fac.append((a[0]))
            
        grouped = self.masterdata_clean.groupby([fac], axis=1)
        sample_grouped = self.sampledata_clean.groupby([fac], axis=1)
        
        data_org = self.masterdata_all
        data1_org = self.masterdata_all
        data_org['Slide_Mean'] = scipy.stats.nanmean(self.masterdata_clean, axis=1)
        data_org['Slide_Median'] = scipy.stats.nanmedian(self.masterdata_clean, axis=1)
        data_org['Slide_Max'] = np.nanmax(self.masterdata_clean_for_max, axis=1)
        data_org['Max_sample'] = self.masterdata_clean_for_max.idxmax(axis=1)
        
        qc_ind_list = []
        with open(conf_file) as f:
            while True:
                text1 = f.readline()
                if text1 == "":
                    break
                cols = text1.strip().split('\t')
                #print len(cols)
                if cols[0] == "pi_qc_list":
                    qc_col_num = -1*((len(cols)-1)+4)
                    for i in range(1, len(cols)):
                        qc_ind_list = []
                        terms = ((cols[i])[1:-1]).strip().split(',')
                        #print (terms)
                        for col in data_org.columns:
                            #print col
                            if re.match(terms[1], col):
                                print col
                                qc_ind_list.append(col)
                                data_org[terms[0]] = scipy.stats.nanmedian(data_org.ix[:,qc_ind_list], axis=1)
                                
        print qc_col_num
        quality_data = data_org.ix[:,qc_col_num:]
        print quality_data.columns
        
        data_org.insert(0, "AB_name", ab_name_column)
        quality_data.insert(0, "AB_name", ab_name_column)
        #pi_data.insert(1, "slide_file", self.masterdata_all["slide_file"])
        quality_data.insert(1, "Slide_file", self.masterdata_all["Slide_file"])
        data_org.insert(2, "Gene_ID", ab_gene_column)
        quality_data.insert(2, "Gene_ID", ab_gene_column)
        data_org.insert(3, "Swiss_ID", ab_swiss_column)
        quality_data.insert(3, "Swiss_ID", ab_swiss_column)
                
        data_org.to_csv(pi_files_path+self.masterfile[:-4]+"/"+self.masterfile[:-4]+"_master_table.xls",sep="\t", na_rep='NA')
        quality_data_no_max_sample = quality_data.drop('Max_sample', 1)
        quality_data_no_max_sample.to_csv(pi_files_path+self.masterfile[:-4]+"/"+self.masterfile[:-4]+"_quality_table.xls",sep="\t", na_rep='NA') 

        for name,group in grouped:
            pidata = grouped.get_group(name)
            pi_sdata = sample_grouped.get_group(name)
            pidata.insert(0, "AB_name", ab_name_column)
            pidata.insert(1, "Slide_file", self.masterdata_all["Slide_file"])
            pidata.insert(2, "Gene_ID", ab_gene_column)
            pidata.insert(3, "Swiss_ID", ab_swiss_column)
            
            pi_sdata.insert(0, "AB_name", "AB_name")
            pi_sdata.insert(1, "Slide_file", "Slide_file")
            pi_sdata.insert(2, "Gene_ID", "Gene_ID")
            pi_sdata.insert(3, "Swiss_ID", "Swiss_ID")
            
            temp = pi_sdata.append(pidata)
            temp.to_csv(pi_files_path+self.masterfile[:-4]+"/"+name+"_"+self.masterfile[:-4]+"_table.xls",sep="\t", na_rep='NA')
            
    
    def raw_data(self,path,start,end,pmt_settings):
        #start='RPPA0013_',end='.gpr'
               
        self.path = path
        self.start = start
        self.end = end
        self.dirs = sorted(os.listdir(self.path))
        master_rawdata_dic = {}
        
        gpr_file_dic = {}

        for folder in self.dirs:                
            if folder.startswith(self.start):
                negative_slide = {}
                self.folder_with_path = os.path.join(self.path, folder)
                self.path2 = self.folder_with_path
                self.files = sorted(os.listdir(self.path2))
                for filename in self.files:
                    if filename.endswith(self.end):
                        # print filename.strip()
                        cols = filename.strip().split('_')
                        pmt = cols[2]
                        # if pmt == pmt_int:
                        gpr_file_dic[filename] = os.path.join(self.path2, filename)

        # print gpr_file_dic
        self.gpr_file_dic = gpr_file_dic

        for setting in pmt_settings:
            self.pmt_int = setting

            ant_sample_dict = OrderedDict()
            #ant_sample_dict390 = OrderedDict()
            data_list = []
            ab_file = {}
            batch = {}
            neg_file = {}
            total_prot_file = {}
            total_prot_file_cycle = {}

            for i in range(0, len(self.slidedata.index.values)):
                # print i
                print (self.slidedata).values[i,5]
                if ((self.slidedata).values[i, 5]).strip() in (gpr_file_dic):
                    # print (self.slidedata).values[i,5]
                    cols = (self.slidedata.values[i, 5]).split('_')
                    pmt = cols[2]

                    if pmt == self.pmt_int:
                        ab_file[self.slidedata.values[i, 1]] = self.slidedata.values[i,5]
                        neg_file[self.slidedata.values[i, 7]] = self.slidedata.values[i,7]
                        total_prot_file[self.slidedata.values[i, 6]] = self.slidedata.values[i,6]
                        gpr_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,5]], sep="\t",skiprows=32)
                        raw_data = gpr_data["F635 Median - B635"].astype(str)
                        flag_data = gpr_data["Flags"]
                        raw_data[flag_data == -100] = "NA"
                        #if pmt == "340PMT":
                        ant_sample_dict[self.slidedata.values[i, 1]] = raw_data

                        #######################################################
                        neg_cols = (self.slidedata.values[i,7]).split("_")
                        neg_id = "Ne"+neg_cols[3]+neg_cols[2]
                        ab_file[neg_id] = self.slidedata.values[i,7]
                        gpr_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,7]], sep="\t",skiprows=32)
                        raw_data = gpr_data["F635 Median - B635"].astype(str)
                        flag_data = gpr_data["Flags"]
                        raw_data[flag_data==-100] = "NA"
                        ant_sample_dict[neg_id] = raw_data

                        if ~pd.isnull(self.slidedata.values[i, 6]):
                            tot_cols = (self.slidedata.values[i,6]).split("_")
                            tot_id = tot_cols[5]+"T"
                            ab_file[tot_id] = self.slidedata.values[i,6]
                            gpr_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,6]], sep="\t",skiprows=32)
                            try:
                                raw_data = gpr_data["F594 Median - B594"].astype(str)
                            except:
                                raw_data = gpr_data["F532 Median - B532"].astype(str)
                            flag_data = gpr_data["Flags"]
                            raw_data[flag_data == -100] = "NA"
                            ant_sample_dict[tot_id] = raw_data
                        ########################################################
                        """
                        #data_list.append(ant_sample_dict[self.slidedata.values[i,1]])
                        #print self.slidedata.values[i,5]
                        if not batch:
                            batch[self.slidedata.values[i,0]] = self.slidedata.values[i,0]
                        else:
                            if (i == len(self.slidedata.index.values)-1) or (self.slidedata.values[i,0] not in batch):
                                batch[self.slidedata.values[i,0]] = self.slidedata.values[i,0]
                                #print neg_file[self.slidedata.values[i,7]]
                                if (i != len(self.slidedata.index.values)-1):
                                    del neg_file[self.slidedata.values[i,7]]
                                for j in neg_file:
                                    #print j
                                    cols = (j).split('_')
                                    #ab_file[cols[5]] = j
                                    ab_file["Ne"+cols[3]] = j
                                    gpr_data = pd.read_table(gpr_file_dic[j], sep="\t",skiprows=32)
                                    raw_data = gpr_data["F635 Median - B635"].astype(str)
                                    flag_data = gpr_data["Flags"]
                                    raw_data[flag_data==-100] = "NA"
                                    ant_sample_dict["Ne"+cols[3]] = raw_data
                                    #ant_sample_dict[cols[5]] = raw_data
                                neg_file = {}
                        
                        if not total_prot_file_cycle:
                            total_prot_file_cycle[self.slidedata.values[i,6]] = self.slidedata.values[i,6]
                        else:
                            if (self.slidedata.values[i,6] not in total_prot_file_cycle) or (i == len(self.slidedata.index.values)-1):
                                total_prot_file_cycle[self.slidedata.values[i,6]] = self.slidedata.values[i,6]
                                if (i != len(self.slidedata.index.values)-1):
                                    del total_prot_file[self.slidedata.values[i,6]]
                                for j in total_prot_file:
                                    #print j
                                    cols = (j).split('_')
                                    ab_file[cols[5]+"T"] = j
                                    gpr_data = pd.read_table(gpr_file_dic[j], sep="\t",skiprows=32)
                                    try:
                                        raw_data = gpr_data["F594 Median - B594"].astype(str)
                                    except:
                                        raw_data = gpr_data["F532 Median - B532"].astype(str)
                                    flag_data = gpr_data["Flags"]
                                    raw_data[flag_data==-100] = "NA"
                                    ant_sample_dict[cols[5]+"T"] = raw_data"""
                                    
                                                                        
            combined_data = pd.DataFrame(ant_sample_dict)
            sample_list = []
            sample_dic = {}
            j = 0
            for i in (gpr_data["Name"]).values:
                #print i
                if j == 0:
                    sam = i+".1"
                    j += 1
                elif j == 1:
                    sam = i+".2"
                    j += 1
                elif j == 2:
                    sam = i+".3"
                    j = 0
                    
                sample_list.append(sam)
                sample_dic[i] = i
    
            combined_data1 = combined_data.copy()        
            combined_data1.index = sample_list
            master_data = combined_data1.transpose()
            ab_file_column = pd.DataFrame(pd.Series(ab_file))
            master_data.insert(0, "Slide_file", ab_file_column)
            master_rawdata_dic[self.pmt_int] = master_data

        return master_rawdata_dic
            
    def flag_data(self,path,start,end,pmt_settings):
            #start='RPPA0013_',end='.gpr'
                   
            self.path = path
            self.start = start
            self.end = end
            self.dirs = sorted(os.listdir(self.path))
            master_flagdata_dic = {}
            
            for setting in pmt_settings:            
                self.pmt_int = setting
            
                gpr_file_dic = {}
        
                for folder in self.dirs:                
                    if folder.startswith(self.start):
                        negative_slide = {}
                        self.folder_with_path = os.path.join(self.path, folder)
                        self.path2 = self.folder_with_path
                        self.files = sorted(os.listdir(self.path2))
                        for filename in self.files:
                            if filename.endswith(self.end):
                                #print filename
                                cols = filename.rstrip().split('_')
                                pmt = cols[2]
                                #if pmt == pmt_int:
                                gpr_file_dic[filename] = os.path.join(self.path2, filename)
                    
                ant_sample_dict = OrderedDict()
                #ant_sample_dict390 = OrderedDict()
                data_list = []
                ab_file = {}
                #ab_file390 = {}
                k = 0
                for i in range(0,len(self.slidedata.index.values)):
                    #print (self.slidedata).values[i,5]
                    if (self.slidedata).values[i,5] in gpr_file_dic:
                        cols = (self.slidedata.values[i,5]).split('_')
                        pmt = cols[2]
                        if pmt == self.pmt_int:
                            ab_file[self.slidedata.values[i,1]] = self.slidedata.values[i,5]                
                            gpr_data = pd.read_table(gpr_file_dic[self.slidedata.values[i,5]], sep="\t",skiprows=32)
                            raw_data = gpr_data["Flags"].astype(str)
                            #flag_data = gpr_data["Flags"]
                            #raw_data[flag_data==-100] = "NA"
                            #if pmt == "340PMT":
                            ant_sample_dict[self.slidedata.values[i,1]] = raw_data
                            
                combined_data = pd.DataFrame(ant_sample_dict)
                sample_list = []
                sample_dic = {}
                j = 0
                for i in gpr_data["Name"].values:
                    #print i
                    if j == 0:
                        sam = i+".1"
                        j += 1
                    elif j == 1:
                        sam = i+".2"
                        j += 1
                    elif j == 2:
                        sam = i+".3"
                        j = 0
                    sample_list.append(sam)
                    sample_dic[i] = i
        
                #combined_data1 = combined_data.set_index(sample_list)
                combined_data1 = combined_data.copy()        
                combined_data1.index = sample_list
                master_data = combined_data1.transpose()
                ab_file_column = pd.DataFrame(pd.Series(ab_file))
                master_data.insert(0, "Slide_file", ab_file_column)
                master_flagdata_dic[self.pmt_int] = master_data
                
            return master_flagdata_dic
        
    def sample_table(self,sample_data,sample_list,experiment):
        f1 =  open(experiment+'_sample_table.xls',"w")
        f1.write("Original")
        for i in (sample_data.columns):
            f1.write("\t"+i)
        f1.write("\n")
        f1.write("AB_ID")
        for i in sample_data.columns:
            a = re.split(r"(\..$)",i)
            if a[0] in (sample_list.index).values:
                #print sample_list.loc[a[0],"Annotation"]
                #print sample_list.loc[a[0],"Sample_ID"]
                f1.write("\t"+str(sample_list.loc[a[0],"Sample_ID"]))
            else:
                f1.write("\t"+i)
        f1.write("\n")
        f1.close()
        return experiment+'_sample_table.xls'

    def pi_directory(self,experiment,pmt_settings,pi_list,mouse_pi_list):

        self.pi_files_path = experiment+"_results/pi_files/"
        self.pi_dir_path = experiment+"_results/pi_dirs/"
    
        norm_data = {}
        raw_data = {}
        flag_data = {}
        quality_data = {}
        
        #pi_list = ["YL"]
        #mouse_pi_list = [""]
        #pi_list = ["AS","BWO","BY","Cal","Ctrl","DPE","HES","JR","SM","SPIKE","YC","YL"]
        #mouse_pi_list = ["BY","BWO","YL"]
        
        
        for pi in pi_list:
            print pi
            sf = open(pi+"_saturated_antibodies.xls","w")
            if not os.path.exists(self.pi_dir_path+pi+'_'+experiment):
                os.makedirs(self.pi_dir_path+pi+'_'+experiment)
            
            for setting in pmt_settings:            
                self.pmt_int = setting
                #for pmt in pmt_list:
                #configfiles = [os.path.join(dirpath, f)
                for dirpath, dirnames, files in os.walk(self.pi_files_path):
                    for f in fnmatch.filter(files, experiment+'*Norm*'+str(self.pmt_int)+'*quality_table*.xls'):
                        files = os.path.join(dirpath, f)
                        quality_data[self.pmt_int] = pd.read_table(files, sep="\t",index_col=0)
                        print f,dirpath
                        shutil.copy(files, self.pi_dir_path+pi+'_'+experiment)
                    
            for setting in pmt_settings:            
                self.pmt_int = setting
                #configfiles = [os.path.join(dirpath, f)
                for dirpath, dirnames, files in os.walk(self.pi_files_path):            
                    for f in fnmatch.filter(files, pi+'_*Norm*'+str(self.pmt_int)+'*.xls'):
                        files = os.path.join(dirpath, f)
                        norm_data[self.pmt_int] = pd.read_table(files, sep="\t",index_col=0)
                        print f,dirpath
                        shutil.copy(files, self.pi_dir_path+pi+'_'+experiment)
                    for f in fnmatch.filter(files, pi+'_*Raw*'+str(self.pmt_int)+'*.xls'):
                        files = os.path.join(dirpath, f)
                        raw_data[self.pmt_int] = pd.read_table(files, sep="\t",index_col=0)
                        print f,dirpath
                        shutil.copy(files, self.pi_dir_path+pi+'_'+experiment)
                    for f in fnmatch.filter(files, pi+'_*Flag*'+str(self.pmt_int)+'*.xls'):
                        files = os.path.join(dirpath, f)
                        flag_data[self.pmt_int] = pd.read_table(files, sep="\t",index_col=0)
                        print f,dirpath
                        shutil.copy(files, self.pi_dir_path+pi+'_'+experiment)
        
            rev_pmt_list = (np.sort(pmt_settings))[::-1]
                
            normdata = (norm_data[rev_pmt_list[0]]).ix[0:1,:]
            rawdata = (raw_data[rev_pmt_list[0]]).ix[0:1,:]
            flagdata = (flag_data[rev_pmt_list[0]]).ix[0:1,:]
            qualitydata = pd.DataFrame(columns=(quality_data[rev_pmt_list[0]]).columns)#Not correct as qulaity doesn't have 2nd raw with sample names
            #qualitydata = (quality_data[rev_pmt_list[0]]).ix[0:1,:]
            #rawdata = pd.DataFrame(columns=(raw_data[rev_pmt_list[0]]).columns)
            #flagdata = pd.DataFrame(columns=(flag_data[rev_pmt_list[0]]).columns)
            #qualitydata = pd.DataFrame(columns=(quality_data[rev_pmt_list[0]]).columns)

            #for ab_id in np.unique(self.slidedata.ix[:,1].values):
            finished_dic = {}

            fac = []
            for item in self.slidedata.index:
                cols = (self.slidedata.ix[item,"Antibody_file"]).split("_")
                fac.append((cols[2])[0:3])
            self.slidedata["PMT"] = fac

            sort_slidedata = self.slidedata.sort(["Experiment","Slide#","PMT"], ascending=[1,1,0])

            for slide in sort_slidedata.index:
                if (sort_slidedata.ix[slide,"Ab_ID"]) not in finished_dic:
                    slide_filt_data = sort_slidedata[sort_slidedata["Ab_ID"]==sort_slidedata.ix[slide,"Ab_ID"]]
                    for slide_id in slide_filt_data.index:
                        ab_id = slide_filt_data.loc[slide_id,"Ab_ID"]
                        print ab_id, pi
                        neg_cols = (slide_filt_data.ix[slide_id,"Negative_file"]).split("_")
                        tot_cols = (slide_filt_data.ix[slide_id,"Total_protein_file"]).split("_")
                        neg_id = "Ne"+neg_cols[3]+neg_cols[2]
                        tot_id = tot_cols[5]+"T"
                        pmt = str(slide_filt_data.ix[slide_id,"PMT"])+"PMT"
                        if (str(ab_id) in (norm_data[pmt]).index):
                            #print(raw_data[pmt]).index
                            if ((((raw_data[pmt].ix[str(ab_id),4:]).astype(float)>=55000).sum() < 6)|
                                (slide_id==slide_filt_data.index[-1])):
                                neg_id = "Ne"+neg_cols[3]+pmt
                                print "neg_id"+neg_id
                                #print "tot_id"+tot_id
                                #print rev_pmt_list[a]
                                #print len(normdata.columns)
                                #print len((norm_data[rev_pmt_list[a]]).columns)
                                #if len((norm_data[rev_pmt_list[a]]).ix[str(ab_id)]) == 1:
                                #print (norm_data[rev_pmt_list[a]]).ix[str(ab_id)]
                                normdata.loc[str(ab_id)] = (norm_data[pmt]).ix[str(ab_id)]
                                rawdata.loc[str(neg_id)] = raw_data[pmt].ix[str(neg_id)]
                                rawdata.loc[str(tot_id)] = raw_data[pmt].ix[str(tot_id)]
                                rawdata.loc[str(ab_id)] = raw_data[pmt].ix[str(ab_id)]
                                flagdata.loc[str(ab_id)] = flag_data[pmt].ix[str(ab_id)]
                                qualitydata.loc[str(ab_id)] = quality_data[pmt].ix[(ab_id)]
                                finished_dic[ab_id] = ab_id
                                break
            
            
            """for slide_id in self.slidedata.index:
                ab_cols = (self.slidedata.values[slide_id,5]).split("_")
                if ab_cols[2] == "460PMT":
                    if slide_id not in finished_dic:
                        finished_dic[slide_id] = slide_id
                        ab_id = self.slidedata.loc[slide_id,"Ab_ID"]
                        print ab_id, pi
                        neg_cols = (self.slidedata.values[slide_id,7]).split("_")
                        tot_cols = (self.slidedata.values[slide_id,6]).split("_")
                        neg_id = "Ne"+neg_cols[3]+neg_cols[2]
                        tot_id = tot_cols[5]+"T"
                        for a in range(0, len(rev_pmt_list)):

                            #print 
                            if (str(ab_id) in (norm_data[rev_pmt_list[a]]).index):
                                if (((raw_data[rev_pmt_list[a]].ix[str(ab_id),4:]).astype(float)>=55000).sum() < 6):
                                    #print str(ab_id)+" valid"
                                    neg_id = "Ne"+neg_cols[3]+rev_pmt_list[a]
                                    print "neg_id"+neg_id
                                    #print "tot_id"+tot_id
                                    #print rev_pmt_list[a]
                                    #print len(normdata.columns)
                                    #print len((norm_data[rev_pmt_list[a]]).columns)
                                    #if len((norm_data[rev_pmt_list[a]]).ix[str(ab_id)]) == 1:
                                    #print (norm_data[rev_pmt_list[a]]).ix[str(ab_id)]
                                    normdata.loc[str(ab_id)] = (norm_data[rev_pmt_list[a]]).ix[str(ab_id)]
                                    rawdata.loc[str(neg_id)] = raw_data[rev_pmt_list[a]].ix[str(neg_id)]
                                    rawdata.loc[str(tot_id)] = raw_data[rev_pmt_list[a]].ix[str(tot_id)]
                                    rawdata.loc[str(ab_id)] = raw_data[rev_pmt_list[a]].ix[str(ab_id)]
                                    flagdata.loc[str(ab_id)] = flag_data[rev_pmt_list[a]].ix[str(ab_id)]
                                    qualitydata.loc[str(ab_id)] = quality_data[rev_pmt_list[a]].ix[(ab_id)]
                                    break

                            if a == len(rev_pmt_list)-1:
                                print str(ab_id)+" saturated"
                                sf.write(str(ab_id)+" saturated\n")
                                for b in range((len(rev_pmt_list)-1),-1,-1):
                                    if (str(ab_id) in (norm_data[rev_pmt_list[b]]).index):
                                        neg_id = "Ne"+neg_cols[3]+rev_pmt_list[b]
                                        #if len((norm_data[rev_pmt_list[b]]).ix[str(ab_id)]) == 1:
                                        normdata.loc[str(ab_id)] = (norm_data[rev_pmt_list[b]]).ix[str(ab_id)]
                                        rawdata.loc[str(neg_id)] = raw_data[rev_pmt_list[b]].ix[str(neg_id)]
                                        rawdata.loc[str(tot_id)] = raw_data[rev_pmt_list[b]].ix[str(tot_id)]
                                        rawdata.loc[str(ab_id)] = raw_data[rev_pmt_list[b]].ix[str(ab_id)]
                                        flagdata.loc[str(ab_id)] = flag_data[rev_pmt_list[b]].ix[str(ab_id)]
                                        qualitydata.loc[str(ab_id)] = quality_data[rev_pmt_list[b]].ix[(ab_id)]
                                        break"""
            sf.close()
        
            if pi not in mouse_pi_list:
                normdata.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Norm_final.xls',sep="\t", na_rep='NA')
                rawdata.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Raw_final.xls',sep="\t", na_rep='NA')
                flagdata.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Flag_final.xls',sep="\t", na_rep='NA')
                qualitydata.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Quality_final.xls',sep="\t", na_rep='NA')
            else:
                normdata_nonmouse = (normdata[~normdata.AB_name.str.contains('_M_')])
                qualitydata_nonmouse = (qualitydata[~qualitydata.AB_name.str.contains('_M_')])
                rawdata_nonmouse = (rawdata[~rawdata.AB_name.str.contains('_M_')])
                normdata_mouse = (normdata.iloc[0:1,:]).append(normdata[normdata.AB_name.str.contains('_M_')])
                qualitydata_mouse = (qualitydata[qualitydata.AB_name.str.contains('_M_')])
                rawdata_mouse = (rawdata.iloc[0:1,:]).append(rawdata[rawdata.AB_name.str.contains('_M_')])
        
                normdata_nonmouse.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Norm_nonmouse_final.xls',sep="\t", na_rep='NA')
                qualitydata_nonmouse.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Quality_nonmouse_final.xls',sep="\t", na_rep='NA')
                rawdata_nonmouse.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Raw_nonmouse_final.xls',sep="\t", na_rep='NA')
                normdata_mouse.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Norm_mouse_final.xls',sep="\t", na_rep='NA')
                qualitydata_mouse.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Quality_mouse_final.xls',sep="\t", na_rep='NA')
                rawdata_mouse.to_csv(self.pi_dir_path+pi+'_'+experiment+'/'+pi+'_'+experiment+'_Raw_mouse_final.xls',sep="\t", na_rep='NA')


    def reports(self, experiment,pi_list,mouse_pi_list,extra_term):

        batch_groups = self.slidedata.groupby([(self.slidedata.ix[:,0]).values],sort=False)
        ordered_batch = np.unique((self.slidedata.ix[:,0]).values)
        batched_data = {}
        neg_list = {}
        batch_no = 0
        for batch in ordered_batch:
            # print batch
            neg_list[batch_no] = []
            batched_data[batch_no] = batch_groups.get_group(batch)
            for item in (batched_data[batch_no]).index:
                cols = ((batched_data[batch_no]).ix[item,7]).split("_")
                neg_id = "Ne"+cols[3]+cols[2]
                if neg_id not in (neg_list[batch_no]):
                    (neg_list[batch_no]).append(neg_id) 
            print (batched_data[batch_no]).index
            batch_no += 1
        print neg_list

        if extra_term:
            ext = extra_term+"_"
        else:
            ext = ""
        for pi in pi_list:
            print pi
            path1 = self.pi_dir_path+pi+'_'+experiment+'/'
            df = {}
            dfv = {}
            dfn = {}
            worksheet = {}
            if pi in mouse_pi_list:
                df[1] = pd.read_table(path1+pi+'_'+experiment+'_Norm_nonmouse_final.xls', sep="\t",index_col=0)
                dfv[1] = df[1][((df[1])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[1]['AB_name'].str.contains("Tot"))|(df[1]['AB_name'].str.contains("Neg"))]
                dfn[1] = df[1][~(df[1])['AB_name'].map(lambda x: x.endswith('_V'))]        
                df[2] = pd.read_table(path1+pi+'_'+experiment+'_Quality_nonmouse_final.xls', sep="\t",index_col=0)
                dfv[2] = df[2][((df[2])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[2]['AB_name'].str.contains("Tot"))|(df[2]['AB_name'].str.contains("Neg"))]
                dfn[2] = df[2][~(df[2])['AB_name'].map(lambda x: x.endswith('_V'))]  
                df[3] = pd.read_table(path1+pi+'_'+experiment+'_Raw_nonmouse_final.xls', sep="\t",index_col=0)
                dfv[3] = df[3][((df[3])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[3]['AB_name'].str.contains("Tot"))|(df[3]['AB_name'].str.contains("Neg"))]
                dfn[3] = df[3][~(df[3])['AB_name'].map(lambda x: x.endswith('_V'))]
                df[4] = pd.read_table(path1+pi+'_'+experiment+'_Norm_mouse_final.xls', sep="\t",index_col=0)
                dfv[4] = df[4][((df[4])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[4]['AB_name'].str.contains("Tot"))|(df[4]['AB_name'].str.contains("Neg"))]
                dfn[4] = df[4][~(df[4])['AB_name'].map(lambda x: x.endswith('_V'))]
                df[5] = pd.read_table(path1+pi+'_'+experiment+'_Quality_mouse_final.xls', sep="\t",index_col=0)
                dfv[5] = df[5][((df[5])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[5]['AB_name'].str.contains("Tot"))|(df[5]['AB_name'].str.contains("Neg"))]
                dfn[5] = df[5][~(df[5])['AB_name'].map(lambda x: x.endswith('_V'))]
                df[6] = pd.read_table(path1+pi+'_'+experiment+'_Raw_mouse_final.xls', sep="\t",index_col=0)
                dfv[6] = df[6][((df[6])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[6]['AB_name'].str.contains("Tot"))|(df[6]['AB_name'].str.contains("Neg"))]
                dfn[6] = df[6][~(df[6])['AB_name'].map(lambda x: x.endswith('_V'))]
            else:
                df[1] = pd.read_table(path1+pi+'_'+experiment+'_Norm_final.xls', sep="\t",index_col=0)
                dfv[1] = df[1][((df[1])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[1]['AB_name'].str.contains("Tot"))|(df[1]['AB_name'].str.contains("Neg"))]
                dfn[1] = df[1][~(df[1])['AB_name'].map(lambda x: x.endswith('_V'))]
                df[2] = pd.read_table(path1+pi+'_'+experiment+'_Quality_final.xls', sep="\t",index_col=0)
                dfv[2] = df[2][((df[2])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[2]['AB_name'].str.contains("Tot"))|(df[2]['AB_name'].str.contains("Neg"))]
                dfn[2] = df[2][~(df[2])['AB_name'].map(lambda x: x.endswith('_V'))]
                df[3] = pd.read_table(path1+pi+'_'+experiment+'_Raw_final.xls', sep="\t",index_col=0)
                dfv[3] = df[3][((df[3])['AB_name'].map(lambda x: x.endswith('_V')))|((df[1])['AB_name'].map(lambda x: x.startswith('AB_')))|
                               (df[3]['AB_name'].str.contains("Tot"))|(df[3]['AB_name'].str.contains("Neg"))]
                dfn[3] = df[3][~(df[3])['AB_name'].map(lambda x: x.endswith('_V'))]


            # Create a workbook and add a worksheet.
            workbook = xlsxwriter.Workbook(path1+pi+'_'+experiment+'_'+ext+'final_report.xlsx')
            bold = workbook.add_format({'bold': True})
            bgcolbold = workbook.add_format({'bold': True})
            bgcolbold.set_pattern(1)  # This is optional when using a solid fill.
            bgcolbold.set_bg_color('yellow')
            bgcol = workbook.add_format({'bold': True})
            bgcol.set_pattern(1)  # This is optional when using a solid fill.
            bgcol.set_bg_color('yellow')

            worksheet[0] = workbook.add_worksheet("Legend")
            worksheet[1] = workbook.add_worksheet("Norm")
            worksheet[2] = workbook.add_worksheet("QI")
            worksheet[3] = workbook.add_worksheet("Raw")

            if pi in mouse_pi_list:
                worksheet[4] = workbook.add_worksheet("Mouse_Norm")
                worksheet[5] = workbook.add_worksheet("Mouse_QI")
                worksheet[6] = workbook.add_worksheet("Mouse_Raw")
        
            if pi in mouse_pi_list:
                worksheet[0].set_column(0, 0, 110)
                worksheet[0].insert_image('A1', experiment+'_Legend_mouse_Page_1.png')
                worksheet[0].set_column(0, 2, 110)
                worksheet[0].insert_image('B1', experiment+'_Legend_mouse_Page_2.png')
            else:
                worksheet[0].set_column(0, 0, 110)
                worksheet[0].insert_image('A1', experiment+'_Legend_human_Page_1.png')
                worksheet[0].set_column(0, 2, 110)
                worksheet[0].insert_image('B1', experiment+'_Legend_human_Page_2.png')
        
            raw_list = [3,6]
            
            for num in range(1,len(dfv)+1):
                if num not in raw_list:
                    table = dfv[num]
                    worksheet[num].set_column(0, 0, 10)
                    worksheet[num].set_column(1, 1, 20)
                    if num in [2,5]:
                        worksheet[num].write(0, 0, "AB_ID", bold)
                    for j in range(0, len(table.columns)):
                        worksheet[num].write(0, j+1, table.columns[j], bold)
                    for i in range(0, len(table.index)):
                        worksheet[num].write(i+1, 0, table.index[i], bold)
                    for i in range(0, len(table.index)):
                        for j in range(0, len(table.columns)):
                            if str(table.values[i,j])=="nan":
                                #print "yes"
                                worksheet[num].write(i+1, j+1, ("NA"))
                            else:
                                worksheet[num].write(i+1, j+1, (table.values[i,j]))
                                #worksheet[num].write(i+1, j+1, (table.values[i,j]))

                else:
                    table = dfv[num]
                    worksheet[num].set_column(0, 0, 10)
                    worksheet[num].set_column(1, 1, 20)
                    for j in range(0, len(table.columns)):
                        worksheet[num].write(0, j+1, table.columns[j], bold)
                        current_row = 1
                    #row = 0
                    batch_no = 0
                    i = 0
                    while i < len(table.index):
                        rppa_id_count = 0
                        #for i in range(0, len(table.index)):
                        if str(table.index[i]).isdigit():
                            print table.index[i]
                            #print (((batched_data[batch_no])["Ab_ID"]).values).astype(str)
                            if table.index[i] in (((batched_data[batch_no])["Ab_ID"]).values).astype(str):
                                worksheet[num].write(current_row, 0, table.index[i], bold)
                                for j in range(0, len(table.columns)):
                                    if str(table.values[i,j])=="nan":
                                        #print "yes"
                                        worksheet[num].write(current_row, j+1, ("NA"))
                                    else:
                                        worksheet[num].write(current_row, j+1, (table.values[i,j]))
                                current_row += 1
                            else:
                                print "batch_"+str(batch_no)
                                for neg_no in (neg_list[(batch_no)]):
                                    if neg_no in table.index:
                                        worksheet[num].write(current_row, 0, neg_no, bgcolbold)
                                        for j in range(0, len(table.columns)):
                                            if str(table.ix[neg_no, j]) == "nan":
                                                # print "yes"
                                                worksheet[num].write(current_row, j+1, ("NA"))
                                            else:
                                                worksheet[num].write(current_row, j+1, table.ix[neg_no,j], bgcolbold)
                                        current_row += 1
                                        rppa_id_count += 1
                                batch_no += 1
                                if rppa_id_count > 0:
                                    for j in range(0, len(table.columns)):
                                        worksheet[num].write(current_row, j+1, table.columns[j], bold)
                                    current_row += 1
                                i -= 1
                                #Printing the 1st ab from next batch
                                """worksheet[num].write(current_row, 0, table.index[i], bold)
                                for j in range(0, len(table.columns)):
                                    if str(table.values[i,j])=="nan":
                                        #print "yes"
                                        worksheet[num].write(current_row, j+1, ("NA"))
                                    else:
                                        worksheet[num].write(current_row, j+1, (table.values[i,j]))
                                current_row += 1"""
                        elif str(table.index[i]).startswith("AB_ID"):
                            worksheet[num].write(current_row, 0, table.index[i], bold)
                            for j in range(0, len(table.columns)):
                                if str(table.values[i, j])=="nan":
                                    # print "yes"
                                    worksheet[num].write(current_row, j+1, ("NA"))
                                else:
                                    worksheet[num].write(current_row, j+1, (table.values[i,j]))
                            current_row += 1
                        elif str(table.index[i]).startswith("Pr"):
                            worksheet[num].write(current_row, 0, table.index[i], bold)
                            for j in range(0, len(table.columns)):
                                if str(table.values[i,j])=="nan":
                                    #print "yes"
                                    worksheet[num].write(current_row, j+1, ("NA"))
                                else:
                                    worksheet[num].write(current_row, j+1, (table.values[i,j]))
                            current_row += 1
                        i += 1
                    for neg_no in (neg_list[(batch_no)]):
                        if neg_no in table.index:
                            worksheet[num].write(current_row, 0, neg_no, bgcolbold)
                            for j in range(0, len(table.columns)):
                                worksheet[num].write(current_row, j+1, table.ix[neg_no,j], bgcolbold)
                            current_row += 1
                                                
                    """
                    for i in range(0, len(table.index)):
                        if ("Neg" in table.values[i,1]):
                                worksheet[num].write(i+1+row, 0, table.index[i], bgcolbold)
                                if (i < len(table.index)-1):
                                    if ("Neg" not in table.values[i+1,1]):
                                        row += 1
                                        for j in range(0, len(table.columns)):
                                            worksheet[num].write(i+1+row, j+1, table.columns[j], bold)
                        else:
                            worksheet[num].write(i+1+row, 0, table.index[i], bold)
                    row = 0
                    for i in range(0, len(table.index)):
                        if ("Neg" in table.values[i,1]):
                            for j in range(0, len(table.columns)):
                                if str(table.values[i,j])=="nan":
                                    #print "yes"
                                    worksheet[num].write(i+1+row, j+1, ("NA"),bgcol)
                                else:
                                    worksheet[num].write(i+1+row, j+1, (table.values[i,j]),bgcol)
                            if (i < len(table.index)-1):
                                if ("Neg" not in table.values[i+1,1]):
                                    row += 1
                        else:
                            for j in range(0, len(table.columns)):
                                if str(table.values[i,j])=="nan":
                                    #print "yes"
                                    worksheet[num].write(i+1+row, j+1, ("NA"))
                                else:
                                    worksheet[num].write(i+1+row, j+1, (table.values[i,j]))"""     
        
            workbook.close() 
            
            #Un-validated antibodies
            workbook = xlsxwriter.Workbook(path1+pi+'_'+experiment+'_'+ext+'_final_report_testing_antibodies.xlsx')
            bold = workbook.add_format({'bold': True})
            bgcolbold = workbook.add_format({'bold': True})
            bgcolbold.set_pattern(1)  # This is optional when using a solid fill.
            bgcolbold.set_bg_color('yellow')
            bgcol = workbook.add_format({'bold': True})
            bgcol.set_pattern(1)  # This is optional when using a solid fill.
            bgcol.set_bg_color('yellow')
            
            worksheet[0] = workbook.add_worksheet("Legend")
            worksheet[1] = workbook.add_worksheet("Norm")
            worksheet[2] = workbook.add_worksheet("QI")
            worksheet[3] = workbook.add_worksheet("Raw")
            
            if pi in mouse_pi_list:
                worksheet[4] = workbook.add_worksheet("Mouse_Norm")
                worksheet[5] = workbook.add_worksheet("Mouse_QI")
                worksheet[6] = workbook.add_worksheet("Mouse_Raw")
        
            if pi in mouse_pi_list:
                worksheet[0].set_column(0, 0, 110)
                worksheet[0].insert_image('A1', experiment+'_Legend_mouse_Page_1.png')
                worksheet[0].set_column(0, 2, 110)
                worksheet[0].insert_image('B1', experiment+'_Legend_mouse_Page_2.png')
            else:
                worksheet[0].set_column(0, 0, 110)
                worksheet[0].insert_image('A1', experiment+'_Legend_human_Page_1.png')
                worksheet[0].set_column(0, 2, 110)
                worksheet[0].insert_image('B1', experiment+'_Legend_human_Page_2.png')
        
            raw_list = [3,6]
            
            for num in range(1,len(dfn)+1):
                if num not in raw_list:
                    table = dfn[num]
                    worksheet[num].set_column(0, 0, 10)
                    worksheet[num].set_column(1, 1, 20)
                    if num in [2,5]:
                        worksheet[num].write(0, 0, "AB_ID", bold)
                    for j in range(0, len(table.columns)):
                        worksheet[num].write(0, j+1, table.columns[j], bold)
                    for i in range(0, len(table.index)):
                        worksheet[num].write(i+1, 0, table.index[i], bold)
                    for i in range(0, len(table.index)):
                        for j in range(0, len(table.columns)):
                            if str(table.values[i,j])=="nan":
                                #print "yes"
                                worksheet[num].write(i+1, j+1, ("NA"))
                            else:
                                worksheet[num].write(i+1, j+1, (table.values[i,j]))
                                #worksheet[num].write(i+1, j+1, (table.values[i,j]))
                                
                else:
                    table = dfn[num]
                    worksheet[num].set_column(0, 0, 10)
                    worksheet[num].set_column(1, 1, 20)
                    for j in range(0, len(table.columns)):
                        worksheet[num].write(0, j+1, table.columns[j], bold)
                        current_row = 1
                    #row = 0
                    batch_no = 0
                    i = 0
                    while i < len(table.index):
                        #for i in range(0, len(table.index)):
                        rppa_id_count = 0
                        if str(table.index[i]).isdigit():
                            print table.index[i]
                            #print (((batched_data[batch_no])["Ab_ID"]).values).astype(str)
                            if table.index[i] in (((batched_data[batch_no])["Ab_ID"]).values).astype(str):
                                worksheet[num].write(current_row, 0, table.index[i], bold)
                                for j in range(0, len(table.columns)):
                                    if str(table.values[i,j])=="nan":
                                        #print "yes"
                                        worksheet[num].write(current_row, j+1, ("NA"))
                                    else:
                                        worksheet[num].write(current_row, j+1, (table.values[i,j]))
                                current_row += 1
                            else:
                                print "batch_"+str(batch_no)
                                for neg_no in (neg_list[(batch_no)]):
                                    if neg_no in table.index:
                                        worksheet[num].write(current_row, 0, neg_no, bgcolbold)
                                        for j in range(0, len(table.columns)):
                                            if str(table.ix[neg_no,j])=="nan":
                                                #print "yes"
                                                worksheet[num].write(current_row, j+1, ("NA"))
                                            else:
                                                worksheet[num].write(current_row, j+1, table.ix[neg_no,j], bgcolbold)
                                        current_row += 1
                                        rppa_id_count += 1
                                batch_no += 1
                                if rppa_id_count > 0:
                                    for j in range(0, len(table.columns)):
                                        worksheet[num].write(current_row, j+1, table.columns[j], bold)
                                    current_row += 1
                                i -= 1
                                """#Printing the 1st ab from next batch
                                worksheet[num].write(current_row, 0, table.index[i], bold)
                                for j in range(0, len(table.columns)):
                                    if str(table.values[i,j])=="nan":
                                        #print "yes"
                                        worksheet[num].write(current_row, j+1, ("NA"))
                                    else:
                                        worksheet[num].write(current_row, j+1, (table.values[i,j]))
                                current_row += 1"""
                        elif str(table.index[i]).startswith("AB_ID"):
                            worksheet[num].write(current_row, 0, table.index[i], bold)
                            for j in range(0, len(table.columns)):
                                if str(table.values[i,j])=="nan":
                                    #print "yes"
                                    worksheet[num].write(current_row, j+1, ("NA"))
                                else:
                                    worksheet[num].write(current_row, j+1, (table.values[i,j]))
                            current_row += 1
                        elif str(table.index[i]).startswith("Pr"):
                            worksheet[num].write(current_row, 0, table.index[i], bold)
                            for j in range(0, len(table.columns)):
                                if str(table.values[i,j])=="nan":
                                    #print "yes"
                                    worksheet[num].write(current_row, j+1, ("NA"))
                                else:
                                    worksheet[num].write(current_row, j+1, (table.values[i,j]))
                            current_row += 1
                        i += 1
                    for neg_no in (neg_list[(batch_no)]):
                        if neg_no in table.index:
                            worksheet[num].write(current_row, 0, neg_no, bgcolbold)
                            for j in range(0, len(table.columns)):
                                if str(table.ix[neg_no,j])=="nan":
                                    #print "yes"
                                    worksheet[num].write(current_row, j+1, ("NA"))
                                else:
                                    worksheet[num].write(current_row, j+1, table.ix[neg_no,j], bgcolbold)
                            current_row += 1

            workbook.close()

        """
        # Bar charts with images and both the ctrls
        pp = PdfPages(self.masterfile[:-4]+'_barcharts_with_max_and_images.pdf')
        imagedata = pd.read_table(experiment+'_image_table.xls', sep="\t",index_col=0)
        print masterfile

        for i in (self.masterdata.index.values):
            print i
            try:
                ind = np.arange(0,len(qc_mean_data.ix[i]),1)
                fig = plt.figure()
                #plt.gcf().subplots_adjust(bottom=0.5)
                #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
                
                ax1 = plt.subplot2grid((2,4), (1,0), colspan=3)
                #ax = fig.add_subplot(1,1,1)
                ax1.bar(ind, qc_mean_data.ix[i].values,yerr=qc_sd_data.ix[i].values,label=qc_mean_data.columns.values,color='b',edgecolor = "none",align='center',log=0,ecolor="black")
                ax1.set_ylabel('Intensity')
                ax1.set_ylim(0)
                #ax.set_yscale('log')
                ax1.set_xticks(ind)
                group_labels = qc_mean_data.columns.values
                ax1.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()
                
                ax2 = plt.subplot2grid((2,4),(0,3), rowspan=2)
                #ax2 = fig.add_subplot(1,2,2)
                #fname = '2014-04-30_GBL 9041374_390PMT_4M_ALDH_Ab5_A1_W1.jpg'
                if "Pr" in str(i):
                    row_id = str(i)
                else:
                    row_id = str(i)+"_"+pmt
                fname = imagedata.loc[row_id,"Image_file"]
                parts = str(fname).split("/")
                image = Image.open(str(fname))
                ax2.imshow(image,aspect='equal')
                ax2.set_title("Scanned_image",fontweight='medium',fontsize=6)
                print parts[-1]
                ax2.set_xlabel(parts[-1],fontsize=3)
                ax2.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
                ax2.set_xticks([]) ### Hides ticks
                ax2.set_yticks([])
                #image.close()
                
                ind = np.arange(0,2*len((bio_rep["Cellmix1"]).ix[i]),2)
                ax3 = plt.subplot2grid((2,4),(0,0))
                #width = 0.5
                bar1 = ax3.bar(ind, (bio_rep["Cellmix1"]).ix[i].values,width=0.5, label=(bio_rep["Cellmix1"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax3.errorbar(ind, (bio_rep["Cellmix1"]).ix[i].values,yerr=(bio_rep_std["Cellmix1"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                bar2 = ax3.bar(ind+0.5, (bio_rep["Cellmix2"]).ix[i].values,width=0.5, label=(bio_rep["Cellmix2"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                ax3.errorbar(ind+0.5, (bio_rep["Cellmix2"]).ix[i].values,yerr=(bio_rep_std["Cellmix2"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                ax3.set_ylabel('Intensity')
                ax3.set_ylim(0)
                ax3.legend( (bar1[0], bar2[0]), ('Cellmix1', 'Cellmix2'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5 )
                #ax.set_yscale('log')
                ax3.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cellmix1"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax3.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax3.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate())
        
                ind = np.arange(0,2*len((bio_rep["Cal_1_HP"]).ix[i]),2)
                ax4 = plt.subplot2grid((2,4),(0,1))
                #ax = fig.add_subplot(1,1,1)
                bar1 = ax4.bar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,width=0.5, label=(bio_rep["Cal_1_HP"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax4.errorbar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,yerr=(bio_rep_std["Cal_1_HP"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                bar2 = ax4.bar(ind+0.5, (bio_rep["Cal_2_HP"]).ix[i].values,width=0.5, label=(bio_rep["Cal_2_HP"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                ax4.errorbar(ind+0.5, (bio_rep["Cal_2_HP"]).ix[i].values,yerr=(bio_rep_std["Cal_2_HP"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #ax4.set_ylabel('Intensity')
                ax4.set_ylim(0)
                ax4.legend( (bar1[0], bar2[0]), ('Cal_1_HP', 'Cal_2_HP'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5 )
                #ax.set_yscale('log')
                ax4.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cal_1_HP"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax4.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax4.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()
                
                ind = np.arange(0,2*len((bio_rep["Cal_1_JC"]).ix[i]),2)
                ax5 = plt.subplot2grid((2,4),(0,2))
                #ax = fig.add_subplot(1,1,1)
                bar1 = ax5.bar(ind, (bio_rep["Cal_1_JC"]).ix[i].values, width=0.5, label=(bio_rep["Cal_1_JC"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax5.errorbar(ind, (bio_rep["Cal_1_JC"]).ix[i].values,yerr=(bio_rep_std["Cal_1_JC"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                bar2 = ax5.bar(ind+0.5, (bio_rep["Cal_2_JC"]).ix[i].values, width=0.5, label=(bio_rep["Cal_2_JC"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                ax5.errorbar(ind+0.5, (bio_rep["Cal_2_JC"]).ix[i].values,yerr=(bio_rep_std["Cal_2_JC"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #ax4.set_ylabel('Intensity')
                ax5.set_ylim(0)
                ax5.legend( (bar1[0], bar2[0]), ('Cal_1_JC', 'Cal_2_JC'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5)
                #ax.set_yscale('log')
                ax5.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cal_1_JC"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax5.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax5.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()
                
                try:
                    plt.suptitle(self.antibodydata.ix[int(i),"PI_name"]+"_Ab-ID = "+str(i)+"_Dilutions = "+str(self.antibodydata.ix[int(i),"Dilution"]),fontweight='bold')
                except:
                    plt.suptitle("Ab-ID = "+str(i),fontweight='bold')
                
                plt.tight_layout()
                plt.subplots_adjust(top=0.90)
                
                pp.savefig(dpi=200,bbox_inches='tight')
                plt.close()
            except:
                print str(i)+" error"            
        pp.close()
        
        #Bar charts with images and one the ctrls
        pp = PdfPages(self.masterfile[:-4]+'_QI_barcharts.pdf')
        imagedata = pd.read_table(experiment+'_image_table.xls', sep="\t",index_col=0)
        print masterfile
        
        for i in (self.masterdata.index.values):
            print i
            try:
                ind = np.arange(0,len(qc_mean_data1.ix[i]),1)
                fig = plt.figure()
                #plt.gcf().subplots_adjust(bottom=0.5)
                #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
                
                #ax1 = plt.subplot2grid((2,4), (1,0), colspan=3)
                ax1 = fig.add_subplot(1,1,1)
                ax1.bar(ind, qc_mean_data1.ix[i].values,yerr=qc_sd_data1.ix[i].values,label=qc_mean_data1.columns.values,color='b',edgecolor = "none",align='center',log=0,ecolor="black")
                ax1.set_ylabel('Intensity')
                ax1.set_ylim(0)
                #ax.set_yscale('log')
                ax1.set_xticks(ind)
                group_labels = qc_mean_data1.columns.values
                ax1.set_xticklabels(group_labels,fontsize=8,rotation=45,ha='right')
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(8)
                #fig.autofmt_xdate()
                try:
                    ax1.set_title(self.antibodydata.ix[int(i),"PI_name"]+"_Ab-ID = "+str(i)+"_Dilutions = "+str(self.antibodydata.ix[int(i),"Dilution"]),fontweight='bold')
                except:
                    ax1.set_title("Ab-ID = "+str(i),fontweight='bold')
                plt.tight_layout()
                plt.subplots_adjust(top=0.90)
                pp.savefig(dpi=200,bbox_inches='tight')
                plt.close()
            except:
                print str(i)+" error"            
        pp.close()

        #Bar charts with images and one ctrls
        pp = PdfPages(self.masterfile[:-4]+'_barcharts_with_max_and_images_Ctrl1.pdf.pdf')
        imagedata = pd.read_table(experiment+'_image_table.xls', sep="\t",index_col=0)
        print masterfile

        for i in (self.masterdata.index.values):
            print i
            try:
                ind = np.arange(0,len(qc_mean_data1.ix[i]),1)
                fig = plt.figure()
                #plt.gcf().subplots_adjust(bottom=0.5)
                #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

                ax1 = plt.subplot2grid((2,4), (1,0), colspan=3)
                #ax = fig.add_subplot(1,1,1)
                ax1.bar(ind, qc_mean_data1.ix[i].values,yerr=qc_sd_data1.ix[i].values,label=qc_mean_data1.columns.values,color='b',edgecolor = "none",align='center',log=0,ecolor="black")
                ax1.set_ylabel('Intensity')
                ax1.set_ylim(0)
                #ax.set_yscale('log')
                ax1.set_xticks(ind)
                group_labels = qc_mean_data1.columns.values
                ax1.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()

                ax2 = plt.subplot2grid((2,4),(0,3), rowspan=2)
                #ax2 = fig.add_subplot(1,2,2)
                #fname = '2014-04-30_GBL 9041374_390PMT_4M_ALDH_Ab5_A1_W1.jpg'
                if "Pr" in str(i):
                    row_id = str(i)
                else:
                    row_id = str(i)+"_"+pmt
                fname = imagedata.loc[row_id,"Image_file"]
                parts = str(fname).split("/")
                image = Image.open(str(fname))
                ax2.imshow(image,aspect='equal')
                ax2.set_title("Scanned_image",fontweight='medium',fontsize=6)
                print parts[-1]
                ax2.set_xlabel(parts[-1],fontsize=3)
                ax2.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
                ax2.set_xticks([]) ### Hides ticks
                ax2.set_yticks([])
                #image.close()

                ind = np.arange(0,2*len((bio_rep["Cellmix1"]).ix[i]),2)
                ax3 = plt.subplot2grid((2,4),(0,0))
                #width = 0.5
                bar1 = ax3.bar(ind, (bio_rep["Cellmix1"]).ix[i].values,width=0.5, label=(bio_rep["Cellmix1"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax3.errorbar(ind, (bio_rep["Cellmix1"]).ix[i].values,yerr=(bio_rep_std["Cellmix1"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #bar2 = ax3.bar(ind+0.5, (bio_rep["Cellmix2"]).ix[i].values,width=0.5, label=(bio_rep["Cellmix2"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                #ax3.errorbar(ind+0.5, (bio_rep["Cellmix2"]).ix[i].values,yerr=(bio_rep_std["Cellmix2"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                ax3.set_ylabel('Intensity')
                ax3.set_ylim(0)
                #ax3.legend( (bar1[0], bar2[0]), ('Cellmix1', 'Cellmix2'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5 )
                #ax.set_yscale('log')
                ax3.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cellmix1"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax3.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax3.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate())

                ind = np.arange(0,2*len((bio_rep["Cal_1_HP"]).ix[i]),2)
                ax4 = plt.subplot2grid((2,4),(0,1))
                #ax = fig.add_subplot(1,1,1)
                bar1 = ax4.bar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,width=0.5, label=(bio_rep["Cal_1_HP"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax4.errorbar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,yerr=(bio_rep_std["Cal_1_HP"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #bar2 = ax4.bar(ind+0.5, (bio_rep["Cal_2_HP"]).ix[i].values,width=0.5, label=(bio_rep["Cal_2_HP"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                #ax4.errorbar(ind+0.5, (bio_rep["Cal_2_HP"]).ix[i].values,yerr=(bio_rep_std["Cal_2_HP"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #ax4.set_ylabel('Intensity')
                ax4.set_ylim(0)
                #ax4.legend( (bar1[0], bar2[0]), ('Cal_1_HP', 'Cal_2_HP'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5 )
                #ax.set_yscale('log')
                ax4.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cal_1_HP"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax4.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax4.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()

                ind = np.arange(0,2*len((bio_rep["Cal_1_JC"]).ix[i]),2)
                ax5 = plt.subplot2grid((2,4),(0,2))
                #ax = fig.add_subplot(1,1,1)
                bar1 = ax5.bar(ind, (bio_rep["Cal_1_JC"]).ix[i].values, width=0.5, label=(bio_rep["Cal_1_JC"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax5.errorbar(ind, (bio_rep["Cal_1_JC"]).ix[i].values,yerr=(bio_rep_std["Cal_1_JC"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #bar2 = ax5.bar(ind+0.5, (bio_rep["Cal_2_JC"]).ix[i].values, width=0.5, label=(bio_rep["Cal_2_JC"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                #ax5.errorbar(ind+0.5, (bio_rep["Cal_2_JC"]).ix[i].values,yerr=(bio_rep_std["Cal_2_JC"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #ax4.set_ylabel('Intensity')
                ax5.set_ylim(0)
                #ax5.legend( (bar1[0], bar2[0]), ('Cal_1_JC', 'Cal_2_JC'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5)
                #ax.set_yscale('log')
                ax5.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cal_1_JC"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax5.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax5.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()

                try:
                    plt.suptitle(self.antibodydata.ix[int(i),"PI_name"]+"_Ab-ID = "+str(i)+"_Dilutions = "+str(self.antibodydata.ix[int(i),"Dilution"]),fontweight='bold')
                except:
                    plt.suptitle("Ab-ID = "+str(i),fontweight='bold')

                plt.tight_layout()
                plt.subplots_adjust(top=0.90)

                pp.savefig(dpi=200,bbox_inches='tight')
                plt.close()
            except:
                print str(i)+" error"
        pp.close()

        #Bar charts without images and ctrls-1
        pp = PdfPages(self.masterfile[:-4]+'_barcharts_without_max_Ctrl1.pdf.pdf')
        imagedata = pd.read_table(experiment+'_image_table.xls', sep="\t",index_col=0)
        print masterfile

        for i in (self.masterdata.index.values):
            print i
            try:
                ind = np.arange(0,len(qc_mean_data1.ix[i]),1)
                fig = plt.figure()
                #plt.gcf().subplots_adjust(bottom=0.5)
                #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

                ax1 = plt.subplot2grid((2,3), (1,0), colspan=3)
                #ax = fig.add_subplot(1,1,1)
                ax1.bar(ind, qc_mean_data1.ix[i].values,yerr=qc_sd_data1.ix[i].values,label=qc_mean_data1.columns.values,color='b',edgecolor = "none",align='center',log=0,ecolor="black")
                ax1.set_ylabel('Intensity')
                ax1.set_ylim(0)
                #ax.set_yscale('log')
                ax1.set_xticks(ind)
                group_labels = qc_mean_data1.columns.values
                ax1.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()

                ind = np.arange(0,2*len((bio_rep["Cellmix1"]).ix[i]),2)
                ax3 = plt.subplot2grid((2,3),(0,0))
                #width = 0.5
                bar1 = ax3.bar(ind, (bio_rep["Cellmix1"]).ix[i].values,width=0.5, label=(bio_rep["Cellmix1"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax3.errorbar(ind, (bio_rep["Cellmix1"]).ix[i].values,yerr=(bio_rep_std["Cellmix1"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #bar2 = ax3.bar(ind+0.5, (bio_rep["Cellmix2"]).ix[i].values,width=0.5, label=(bio_rep["Cellmix2"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                #ax3.errorbar(ind+0.5, (bio_rep["Cellmix2"]).ix[i].values,yerr=(bio_rep_std["Cellmix2"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                ax3.set_ylabel('Intensity')
                ax3.set_ylim(0)
                #ax3.legend( (bar1[0], bar2[0]), ('Cellmix1', 'Cellmix2'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5 )
                #ax.set_yscale('log')
                ax3.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cellmix1"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix2"]).columns.values
                group_labels = fac
                ax3.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax3.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate())

                ind = np.arange(0,2*len((bio_rep["Cal_1_HP"]).ix[i]),2)
                ax4 = plt.subplot2grid((2,3),(0,1))
                #ax = fig.add_subplot(1,1,1)
                bar1 = ax4.bar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,width=0.5, label=(bio_rep["Cal_1_HP"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax4.errorbar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,yerr=(bio_rep_std["Cal_1_HP"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #bar2 = ax4.bar(ind+0.5, (bio_rep["Cal_2_HP"]).ix[i].values,width=0.5, label=(bio_rep["Cal_2_HP"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                #ax4.errorbar(ind+0.5, (bio_rep["Cal_2_HP"]).ix[i].values,yerr=(bio_rep_std["Cal_2_HP"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #ax4.set_ylabel('Intensity')
                ax4.set_ylim(0)
                #ax4.legend( (bar1[0], bar2[0]), ('Cal_1_HP', 'Cal_2_HP'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5 )
                #ax.set_yscale('log')
                ax4.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cal_1_HP"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax4.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax4.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()

                ind = np.arange(0,2*len((bio_rep["Cal_1_JC"]).ix[i]),2)
                ax5 = plt.subplot2grid((2,3),(0,2))
                #ax = fig.add_subplot(1,1,1)
                bar1 = ax5.bar(ind, (bio_rep["Cal_1_JC"]).ix[i].values, width=0.5, label=(bio_rep["Cal_1_JC"]).columns.values,color='b',edgecolor = "none",align='center',log=0)
                ax5.errorbar(ind, (bio_rep["Cal_1_JC"]).ix[i].values,yerr=(bio_rep_std["Cal_1_JC"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #bar2 = ax5.bar(ind+0.5, (bio_rep["Cal_2_JC"]).ix[i].values, width=0.5, label=(bio_rep["Cal_2_JC"]).columns.values,color='r',edgecolor = "none",align='center',log=0)
                #ax5.errorbar(ind+0.5, (bio_rep["Cal_2_JC"]).ix[i].values,yerr=(bio_rep_std["Cal_2_JC"]).ix[i].values,fmt='none',ecolor='black',elinewidth=0.5, capsize=1,capthick=0.5)
                #ax4.set_ylabel('Intensity')
                ax5.set_ylim(0)
                #ax5.legend( (bar1[0], bar2[0]), ('Cal_1_JC', 'Cal_2_JC'),prop={'size':4},loc='best',fancybox=True,framealpha=0.5)
                #ax.set_yscale('log')
                ax5.set_xticks(ind+0.25)
                fac = []
                for sam in (bio_rep["Cal_1_JC"]).columns:
                    #print i
                    cols = sam.split("_")
                    fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])
                #group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax5.set_xticklabels(group_labels,fontsize=4,rotation=45,ha='right')
                for tick in ax5.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                #fig.autofmt_xdate()

                try:
                    plt.suptitle(self.antibodydata.ix[int(i),"PI_name"]+"_Ab-ID = "+str(i)+"_Dilutions = "+str(self.antibodydata.ix[int(i),"Dilution"]),fontweight='bold')
                except:
                    plt.suptitle("Ab-ID = "+str(i),fontweight='bold')

                plt.tight_layout()
                plt.subplots_adjust(top=0.90)

                pp.savefig(dpi=200,bbox_inches='tight')
                plt.close()
            except:
                print str(i)+" error"
        pp.close()"""

    def cv(self, experiment, masterfile, pmt):

        masterdata1 = pd.read_table(masterfile, sep="\t", index_col=0)
        masterdata = masterdata1.drop('Slide_file', axis=1)

        fac3 = []
        for i in masterdata.columns.values:
            a = i.split('.')
            b = (a[-1]).split('.')
            fac3.append((a[0]))

        grouped2 = masterdata.groupby([fac3], axis=1)

        gr_mean = grouped2.aggregate(np.mean)
        gr_sd = grouped2.aggregate(np.std)
        gr_median = grouped2.aggregate(np.median)
        gr_cv = gr_sd*100/gr_mean

        f1 = open('cv_protein.txt','w')
        f1.write("Sample\tantiboy\tgr_mean\tgr_median\tCV\n")
        for i in gr_mean.columns.values:
            for j in gr_mean.index.values:
                f1.write(str(i)+"\t"+str(j)+"\t"+str(gr_mean.loc[j,i])+"\t"+str(gr_median.loc[j,i])+"\t"+str(gr_cv.loc[j,i])+"\n")
                # print i,j,gr_mean.loc[j,i],gr_cv.loc[j,i]
        f1.close()

        cvdata = pd.read_table("cv_protein.txt", sep="\t", index_col=0)

        f1 = open(masterfile[:-4]+"_CV_protein_per_antibody.xls",'w')
        f1.write("AB_name\ttotal no. of samples\tmean_signal\tmedian_signal\tsample(signal<200)\tCV<=5\tCV<=10\tCV<=20\tCV>=20")
        f1.write("\tsample(200<signal<500)\tCV<=5\tCV<=10\tCV<=20\tCV>=20")
        f1.write("\tsample(500<signal<1000)\tCV<=5\tCV<=10\tCV<=20\tCV>=20")
        f1.write("\tsample(1000<signal)\tCV<=5\tCV<=10\tCV<=20\tCV>=20")
        f1.write("\tsample(200<signal)\tCV<=5\tCV<=10\tCV<=20\tCV>=20")
        f1.write("\tsample(500<signal)\tCV<=5\tCV<=10\tCV<=20\tCV>=20\n")
        for item in np.unique(cvdata.antiboy.values):
            total_samples = (len(cvdata[cvdata.antiboy==item].index))
            mean_signal = np.mean(cvdata[cvdata.antiboy==item].gr_mean)
            median_signal = np.median(cvdata[cvdata.antiboy==item].gr_mean)
            sample_less_200 = (len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean<=200)].index))
            if sample_less_200 > 0:
                cv_less_5_less_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean<=200)&(cvdata.CV<=5)].index)*1.00/sample_less_200)
                cv_less_10_less_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean<=200)&(cvdata.CV<=10)].index)*1.00/sample_less_200)
                cv_less_20_less_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean<=200)&(cvdata.CV<=20)].index)*1.00/sample_less_200)
                cv_great_20_less_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean<=200)&(cvdata.CV>20)].index)*1.00/sample_less_200)
            else:
                cv_less_5_less_200 = "NA"
                cv_less_10_less_200 = "NA"
                cv_less_20_less_200 = "NA"
                cv_great_20_less_200 = "NA"

            sample_200_500 = (len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.gr_mean<=500)].index))
            if sample_200_500 > 0:
                cv_less_5_200_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.gr_mean<=500)&(cvdata.CV<=5)].index)*1.00/sample_200_500)
                cv_less_10_200_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.gr_mean<=500)&(cvdata.CV<=10)].index)*1.00/sample_200_500)
                cv_less_20_200_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.gr_mean<=500)&(cvdata.CV<=20)].index)*1.00/sample_200_500)
                cv_great_20_200_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.gr_mean<=500)&(cvdata.CV>20)].index)*1.00/sample_200_500)
            else:
                cv_less_5_200_500 = "NA"
                cv_less_10_200_500 = "NA"
                cv_less_20_200_500 = "NA"
                cv_great_20_200_500 = "NA"
            sample_500_1000 = (len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.gr_mean<=1000)].index))
            if sample_500_1000 > 0:
                cv_less_5_500_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.gr_mean<=1000)&(cvdata.CV<=5)].index)*1.00/sample_500_1000)
                cv_less_10_500_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.gr_mean<=1000)&(cvdata.CV<=10)].index)*1.00/sample_500_1000)
                cv_less_20_500_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.gr_mean<=1000)&(cvdata.CV<=20)].index)*1.00/sample_500_1000)
                cv_great_20_500_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.gr_mean<=1000)&(cvdata.CV>20)].index)*1.00/sample_500_1000)
            else:
                cv_less_5_500_1000 = "NA"
                cv_less_10_500_1000 = "NA"
                cv_less_20_500_1000 = "NA"
                cv_great_20_500_1000 = "NA"
            sample_great_1000 = (len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>1000)].index))
            if sample_great_1000 > 0:
                cv_less_5_great_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>1000)&(cvdata.CV<=5)].index)*1.00/sample_great_1000)
                cv_less_10_great_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>1000)&(cvdata.CV<=10)].index)*1.00/sample_great_1000)
                cv_less_20_great_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>1000)&(cvdata.CV<=20)].index)*1.00/sample_great_1000)
                cv_great_20_great_1000 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>1000)&(cvdata.CV>20)].index)*1.00/sample_great_1000)
            else:
                cv_less_5_great_1000 = "NA"
                cv_less_10_great_1000 = "NA"
                cv_less_20_great_1000 = "NA"
                cv_great_20_great_1000 = "NA"

            sample_great_200 = (len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)].index))
            if sample_great_200 > 0:
                cv_less_5_great_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.CV<=5)].index)*1.00/sample_great_200)
                cv_less_10_great_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.CV<=10)].index)*1.00/sample_great_200)
                cv_less_20_great_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.CV<=20)].index)*1.00/sample_great_200)
                cv_great_20_great_200 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>200)&(cvdata.CV>20)].index)*1.00/sample_great_200)
            else:
                cv_less_5_great_200 = "NA"
                cv_less_10_great_200 = "NA"
                cv_less_20_great_200 = "NA"
                cv_great_20_great_200 = "NA"

            sample_great_500 = (len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)].index))
            if sample_great_500 > 0:
                cv_less_5_great_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.CV<=5)].index)*1.00/sample_great_500)
                cv_less_10_great_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.CV<=10)].index)*1.00/sample_great_500)
                cv_less_20_great_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.CV<=20)].index)*1.00/sample_great_500)
                cv_great_20_great_500 =str(len(cvdata[(cvdata.antiboy==item)&(cvdata.gr_mean>500)&(cvdata.CV>20)].index)*1.00/sample_great_500)
            else:
                cv_less_5_great_500 = "NA"
                cv_less_10_great_500 = "NA"
                cv_less_20_great_500 = "NA"
                cv_great_20_great_500 = "NA"

            f1.write(str(item)+"\t"+str(total_samples)+"\t"+str(mean_signal)+"\t"+str(median_signal)+
                     "\t"+str(sample_less_200*1.00/total_samples)+"\t"+cv_less_5_less_200+"\t"+cv_less_10_less_200+
                     "\t"+cv_less_20_less_200+"\t"+cv_great_20_less_200+
                     "\t"+str(sample_200_500*1.00/total_samples)+"\t"+cv_less_5_200_500+"\t"+cv_less_10_200_500+
                     "\t"+cv_less_20_200_500+"\t"+cv_great_20_200_500+
                     "\t"+str(sample_500_1000*1.00/total_samples)+"\t"+cv_less_5_500_1000+"\t"+cv_less_10_500_1000+
                     "\t"+cv_less_20_500_1000+"\t"+cv_great_20_500_1000+
                     "\t"+str(sample_great_1000*1.00/total_samples)+"\t"+cv_less_5_great_1000+"\t"+cv_less_10_great_1000+
                     "\t"+cv_less_20_great_1000+"\t"+cv_great_20_great_1000+
                     "\t"+str(sample_great_200*1.00/total_samples)+"\t"+cv_less_5_great_200+"\t"+cv_less_10_great_200+
                     "\t"+cv_less_20_great_200+"\t"+cv_great_20_great_200+
                     "\t"+str(sample_great_500*1.00/total_samples)+"\t"+cv_less_5_great_500+"\t"+cv_less_10_great_500+
                     "\t"+cv_less_20_great_500+"\t"+cv_great_20_great_500+
                     "\n")
        f1.close()

    def imagedf(self, experiment, image_path):
        f1 = open(experiment+'_image_table.xls', "w")
        f1.write("Ab_ID_and_PMT\tSlide_number\tImage_file\n")
        # path = os.getcwd()
        self.path = image_path
        dirs = sorted(os.listdir(self.path))
        for folder in dirs:
            if folder.startswith(experiment+'_'):
                negative_slide = {}
                folder_with_path = os.path.join(self.path, folder)
                path2 = folder_with_path
                files = sorted(os.listdir(path2))
                no_of_files = 0
                for filename in files:
                    if (filename.endswith('W1.jpg')) or (filename.endswith(
                                                         'W2.jpg')):
                        print filename
                        cols = filename.rstrip().split('_')
                        if (cols[5])[0:2] == "Ne":
                            no_of_files += 1
                            im_pmt = cols[2]
                            slide_number = (cols[3])[:-1]
                            ant_id = "Ne"+cols[3]+im_pmt
                            f1.write(str(ant_id)+"_"+(im_pmt)+"\t"+
                                     slide_number + "\t" +
                                     os.path.join(path2, filename)+"\n")
                        elif (cols[5])[0:2] == "Pr":
                            no_of_files += 1
                            im_pmt = cols[2]
                            slide_number = (cols[3])[:-1]
                            ant_id = cols[5]+"T"
                            f1.write(str(ant_id)+"\t"+slide_number+"\t" +
                                     os.path.join(path2, filename)+"\n")
                        else:
                            no_of_files += 1
                            im_pmt = cols[2]
                            slide_number = (cols[3])[:-1]
                            ant_id = int((cols[5])[2:])
                            pi_name = self.antibodydata.loc[ant_id, "PI_name"]
                            f1.write(str(ant_id)+"_"+(im_pmt)+"\t"+slide_number
                                     + "\t"+os.path.join(path2, filename)+"\n")
                print("No. of files in the "+folder+" = "+str(no_of_files))
        f1.close()
        imagedata = pd.read_table(experiment+'_image_table.xls',
                                  sep="\t", index_col=0)
        return imagedata

    def ctrl_bardata(self, experiment, masterfile, conf_file, drop_items,
                     pi=None, master_data_dic=None, datatype=None):

        if pi is None:
            self.masterfile = masterfile
            self.masterdata_all = pd.read_table(self.masterfile, sep="\t",
                                                index_col=0)
            self.masterdata_raw = self.masterdata_all.ix[:, 1:]
            
            drop_list = []
            if drop_items is not None:
                for item in (drop_items[0]):
                    for i in self.masterdata_raw.columns:
                        #print drop_items[0],item,i
                        if re.match(item, i):
                            drop_list.append(i)
                            print i
            
            self.masterdata_clean = self.masterdata_raw.drop(drop_list, axis=1)
                            
            drop_max_list = []
            drop_items_for_max = drop_items[1]
            if drop_items_for_max is not None:
                for item in drop_items_for_max:
                    for i in self.sampledata_clean.columns:
                        if re.match(item, i):
                            drop_max_list.append(i)
                            print i
            
            self.masterdata_clean_for_max = self.masterdata_raw.drop(drop_list, axis=1)
                        
            pi = ""
        else:
            ctrl_data_all = pd.read_table(
                experiment+"_results/pi_dirs/Ctrl_"+experiment+"/Ctrl_" +
                experiment+"_"+datatype+"_final.xls", sep="\t", index_col=0, skiprows=1)
            ctrl_data = ctrl_data_all[(ctrl_data_all)['AB_name'].map(
                                      lambda x: x.endswith('_V'))]


            i = 0
            for ab in ctrl_data.index:
                # print "ab = "+str(ab)
                cols = ctrl_data.loc[ab, "Slide_file"].split("_")
                # print ctrl_data.loc[ab, "Slide_file"]
                # print cols[2]
                if i == 0:
                    final_data = pd.DataFrame(index=ctrl_data.index,
                                              columns=(master_data_dic[cols[2]]).columns,
                                              dtype=float)
                    # print master_data_dic[cols[2]]
                    # final_data = (master_data_dic[cols[2]]).iloc[0:1, :]
                # else:
                final_data.loc[ab] = (master_data_dic[cols[2]]).loc[ab]
                i += 1

            # print "test"
            # print len((final_data).index)
            # final_data.convert_objects(convert_numeric=True)
            self.masterdata_raw = final_data
            drop_list = []
            if drop_items is not None:
                for item in (drop_items[0]):
                    for i in self.masterdata_raw.columns:
                        #print drop_items[0],item,i
                        if re.match(item, i):
                            drop_list.append(i)
                            print i
            self.masterdata_clean = self.masterdata_raw.drop(drop_list, axis=1)

            drop_max_list = []
            drop_items_for_max = drop_items[1]
            if drop_items_for_max is not None:
                for item in drop_items_for_max:
                    for i in self.sampledata_clean.columns:
                        if re.match(item, i):
                            drop_max_list.append(i)
                            print i

            self.masterdata_clean_for_max = self.masterdata_raw.drop(drop_list, axis=1)
            # print "test"
            # print final_data.index
            # print final_data.dtypes
            pi = "pi_"

        qc_mean_data = pd.DataFrame(index=self.masterdata_clean.index)
        qc_sd_data = pd.DataFrame(index=self.masterdata_clean.index)
        qc_mean_data['Slide_Mean'] = scipy.stats.nanmean(self.masterdata_clean,
                                                         axis=1)
        qc_sd_data['Slide_Mean'] = 0
        qc_mean_data['Slide_Median'] = scipy.stats.nanmedian(self.masterdata_clean,
                                                             axis=1)
        qc_sd_data['Slide_Median'] = 0
        qc_mean_data['Slide_Max'] = np.nanmax(self.masterdata_clean_for_max, axis=1)
        qc_sd_data['Slide_Max'] = 0
        qc_ind_list = []
        with open(conf_file) as f:
            while True:
                text1 = f.readline()
                if text1 == "":
                    break
                cols = text1.strip().split('\t')
                print len(cols)
                if cols[0] == pi+"qc_list":
                    qc_col_num = -1*((len(cols)-1)+4)
                    for i in range(1, len(cols)):
                        qc_ind_list = []
                        terms = ((cols[i])[1:-1]).strip().split(',')
                        # print (terms)
                        for col in self.masterdata_clean.columns:
                            # print col
                            if re.match(terms[1], col):
                                # print col
                                qc_ind_list.append(col)
                        qc_mean_data[terms[0]] = scipy.stats.nanmedian(
                            self.masterdata_clean.ix[:, qc_ind_list], axis=1)
                        qc_sd_data[terms[0]] = scipy.stats.nanstd(
                            self.masterdata_clean.ix[:, qc_ind_list], axis=1)

        """qc_mean_data1 = pd.DataFrame(index=self.masterdata.index)
        qc_sd_data1 = pd.DataFrame(index=self.masterdata.index)
        qc_mean_data1['Slide_Mean'] = scipy.stats.nanmean(self.masterdata,
                                                          axis=1)
        qc_sd_data1['Slide_Mean'] = 0
        qc_mean_data1['Slide_Median'] = scipy.stats.nanmedian(self.masterdata,
                                                              axis=1)
        qc_sd_data1['Slide_Median'] = 0
        # qc_mean_data1['Slide_Max'] = np.nanmax(self.masterdata, axis=1)
        # qc_sd_data1['Slide_Max'] = 0
        qc_ind_list1 = []
        with open(conf_file) as f:
            while True:
                text1 = f.readline()
                if text1 == "":
                    break
                cols = text1.strip().split('\t')
                print len(cols)
                if cols[0] == "qc_list":
                    qc_col_num = -1*((len(cols)-1)+4)
                    for i in range(1, len(cols)):
                        qc_ind_list1 = []
                        terms = ((cols[i])[1:-1]).strip().split(',')
                        # print (terms)
                        for col in self.masterdata.columns:
                            # print col
                            if re.match(terms[1], col):
                                # print col
                                qc_ind_list1.append(col)
                        qc_mean_data1[terms[0]] = scipy.stats.nanmedian(
                            self.masterdata.ix[:, qc_ind_list1], axis=1)
                        qc_sd_data1[terms[0]] = scipy.stats.nanstd(
                            self.masterdata.ix[:, qc_ind_list1], axis=1)"""

        cellmix_data = {}
        cal_data = {}
        bio_rep = {}
        bio_rep_std = {}

        with open(conf_file) as f:
            while True:
                text1 = f.readline()
                if text1 == "":
                    break
                cols = text1.strip().split('\t')
                print len(cols)
                if cols[0] == "cellmix_conc":
                    conc_list = ((cols[1])[1:-1]).strip().split(",")

        with open(conf_file) as f:
            while True:
                text1 = f.readline()
                if text1 == "":
                    break
                cols = text1.strip().split('\t')
                print len(cols)
                if cols[0] == pi+"cellmix":
                    # qc_col_num = -1*((len(cols)-1)+4)
                    for i in range(1, len(cols)):
                        cellmix_list = []
                        terms = ((cols[i])[1:-1]).strip().split(',')
                        print (terms)
                        for col in self.masterdata_clean.columns:
                            # print col
                            if re.match(terms[1], col):
                                # print col
                                cellmix_list.append(col)
                        cellmix_data[terms[0]] = self.masterdata_clean.ix[
                            :, cellmix_list]
                        fac = []
                        for i in (cellmix_data[terms[0]]).columns.values:
                            a = re.split(r"(\..$)", i)
                            fac.append(a[0])
                        grouped = (cellmix_data[terms[0]]).groupby([fac],
                                                                   axis=1,
                                                                   sort=False)
                        bio_rep[terms[0]] = (grouped.aggregate(np.median))
                        (bio_rep[terms[0]]) = (bio_rep[terms[0]]).reindex_axis(
                            sorted((bio_rep[terms[0]]).columns), axis=1)
                        (bio_rep[terms[0]]).columns = conc_list
                        bio_rep_std[terms[0]] = grouped.aggregate(np.std)
                        (bio_rep_std[terms[0]]) = (bio_rep_std[terms[0]]).reindex_axis(
                            sorted((bio_rep_std[terms[0]]).columns), axis=1)
                        (bio_rep_std[terms[0]]).columns = conc_list

                if (cols[0] == pi+"Cal_HP") or (cols[0] == pi+"Cal_JC"):
                    # qc_col_num = -1*((len(cols)-1)+4)
                    for i in range(1, len(cols)):
                        cal_list = []
                        terms = (cols[i]).strip().split(',')
                        print (terms)
                        for col in self.masterdata_clean.columns:
                            # print col
                            if re.match(terms[0], col):
                                # print col
                                cal_list.append(col)
                        cal_data[terms[0]] = self.masterdata_clean.ix[:, cal_list]
                        fac = []
                        for i in (cal_data[terms[0]]).columns.values:
                            a = re.split(r"(\..$)", i)
                            fac.append(a[0])
                        grouped = (cal_data[terms[0]]).groupby([fac], axis=1,
                                                               sort=False)
                        bio_rep[terms[0]] = (grouped.aggregate(np.median))
                        bio_rep_std[terms[0]] = grouped.aggregate(np.std)

        if pi == "pi_":
            with open(conf_file) as f:
                while True:
                    text1 = f.readline()
                    if text1 == "":
                        break
                    cols = text1.strip().split('\t')
                    print len(cols)
                    if cols[0] == "drop_pi_Cal_HP":
                        for i in range(2, len(cols)):
                            drop_list = []
                            terms = (cols[i]).strip().split(',')
                            print (terms)
                            for col in (bio_rep[cols[1]]).columns:
                                if terms[0] in col:
                                    drop_list.append(col)
                                    print col
                        (bio_rep[cols[1]]) = (bio_rep[cols[1]]).drop(drop_list,
                                                                     axis=1)
                        (bio_rep_std[cols[1]]) = (bio_rep_std[cols[1]]).drop(drop_list,
                                                                             axis=1)
                    if cols[0] == "drop_pi_Cal_JC":
                        for i in range(2, len(cols)):
                            drop_list = []
                            terms = (cols[i]).strip().split(',')
                            print (terms)
                            for col in (bio_rep[cols[1]]).columns:
                                if terms[0] in col:
                                    drop_list.append(col)
                                    print col
                        (bio_rep[cols[1]]) = (bio_rep[cols[1]]).drop(drop_list,
                                                                     axis=1)
                        (bio_rep_std[cols[1]]) = (bio_rep_std[cols[1]]).drop(drop_list,
                                                                             axis=1)

        return qc_mean_data, qc_sd_data, bio_rep, bio_rep_std

    def combined_plots(self, qc_mean_data, qc_sd_data, bio_rep, bio_rep_std,
                       pmt, pdf_name, cellmix_ctrls=None, calhp_ctrls=None,
                       caljc_ctrls=None, imagedata=None, bars_for_pi=None):
        qc_mean_bar_data = qc_mean_data
        qc_sd_bar_data = qc_sd_data
        color1 = "b"
        color2 = "b"
        print pdf_name

        if cellmix_ctrls is None:
            if ("Cellmix1" in bio_rep) & ("Cellmix2" in bio_rep):
                cellmix_ctrls = "both"
                color2 = "r"
            elif "Cellmix1" in bio_rep:
                    cellmix_ctrls = "1"
            elif "Cellmix2" in bio_rep:
                    cellmix_ctrls = "2"
        print "cellmix_ctrls="+cellmix_ctrls

        if calhp_ctrls is None:
            if ("Cal_1_HP" in bio_rep) & ("Cal_2_HP" in bio_rep):
                calhp_ctrls = "both"
            elif "Cal_1_HP" in bio_rep:
                    calhp_ctrls = "1"
            elif "Cal_2_HP" in bio_rep:
                    calhp_ctrls = "2"
        print "calhp_ctrls="+str(calhp_ctrls)

        if caljc_ctrls is None:
            if ("Cal_1_JC" in bio_rep) & ("Cal_2_JC" in bio_rep):
                caljc_ctrls = "both"
                color2 = "r"
            elif "Cal_1_JC" in bio_rep:
                    caljc_ctrls = "1"
            elif "Cal_2_JC" in bio_rep:
                    caljc_ctrls = "2"
        print "caljc_ctrls="+str(caljc_ctrls)

        pp = PdfPages(str(pdf_name))
        for i in (qc_mean_bar_data.index.values):
            if imagedata is not None:
                    n_row = 2
                    n_col = 4
            else:
                    n_row = 2
                    n_col = 3
            if cellmix_ctrls == "both":
                    adj = 0.5
            else:
                    adj = 0

            print i
            # try:
            ind = np.arange(0, len(qc_mean_bar_data.ix[i]), 1)
            fig = plt.figure()
            # plt.gcf().subplots_adjust(bottom=0.5)
            # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

            ax1 = plt.subplot2grid((n_row, n_col), (1, 0), colspan=3)
            # ax = fig.add_subplot(1,1,1)
            ax1.bar(ind, qc_mean_bar_data.ix[i].values, yerr=qc_sd_bar_data.ix[i].values,
                    label=qc_mean_bar_data.columns.values, color='b',
                    edgecolor="none", align='center', log=0, ecolor="black")
            ax1.set_ylabel('Intensity')
            ax1.set_ylim(0)
            # ax.set_yscale('log')
            ax1.set_xticks(ind)
            group_labels = qc_mean_bar_data.columns.values
            ax1.set_xticklabels(group_labels, fontsize=4,
                                rotation=45, ha='right')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(4)
            # fig.autofmt_xdate()

            if imagedata is not None:
                try:
                    ax2 = plt.subplot2grid((n_row, n_col), (0, 3), rowspan=2)
                    # ax2 = fig.add_subplot(1,2,2)
                    # fname = '2014-04-30_GBL 9041374_390PMT_4M_ALDH_Ab5_A1_W1.jpg'
                    row_id = str(i)+"_"+pmt
                    fname = imagedata.loc[row_id, "Image_file"]
                    parts = fname.split("\\")
                    image = Image.open(fname)
                    ax2.imshow(image, aspect='equal')
                    ax2.set_title("Scanned_image", fontweight='medium', fontsize=6)
                    print parts[-1]
                    ax2.set_xlabel(parts[-1], fontsize=3)
                    ax2.tick_params(labelcolor='none', top='off', bottom='off',
                                    left='off', right='off')
                    ax2.set_xticks([])  # Hides ticks
                    ax2.set_yticks([])
                except:
                    print "No image for "+str(i)

            ax3 = plt.subplot2grid((n_row, n_col), (0, 0))
            # width = 0.5
            if (cellmix_ctrls == "both") or (cellmix_ctrls == "1"):
                ind = np.arange(0, 2*len((bio_rep["Cellmix1"]).ix[i]), 2)
                bar1 = ax3.bar(ind, (bio_rep["Cellmix1"]).ix[i].values,
                               width=0.5,
                               label=(bio_rep["Cellmix1"]).columns.values,
                               color=color1, edgecolor = "none", align='center',
                               log=0)
                ax3.errorbar(ind, (bio_rep["Cellmix1"]).ix[i].values,
                             yerr=(bio_rep_std["Cellmix1"]).ix[i].values,
                             fmt='none', ecolor='black', elinewidth=0.5,
                             capsize=1, capthick=0.5)
                group_labels = (bio_rep["Cellmix1"]).columns.values
            if (cellmix_ctrls == "both") or (cellmix_ctrls == "2"):
                ind = np.arange(0, 2*len((bio_rep["Cellmix2"]).ix[i]), 2)
                bar2 = ax3.bar(ind+adj, (bio_rep["Cellmix2"]).ix[i].values,
                               width=0.5,
                               label=(bio_rep["Cellmix2"]).columns.values,
                               color=color2, edgecolor = "none", align='center',
                               log=0)
                ax3.errorbar(ind+adj, (bio_rep["Cellmix2"]).ix[i].values,
                             yerr=(bio_rep_std["Cellmix2"]).ix[i].values,
                             fmt='none', ecolor='black', elinewidth=0.5,
                             capsize=1, capthick=0.5)
                group_labels = (bio_rep["Cellmix2"]).columns.values
            ax3.set_ylabel('Intensity')
            ax3.set_ylim(0)
            if (cellmix_ctrls == "both"):
                ax3.legend((bar1[0], bar2[0]), ('Cellmix1', 'Cellmix2'),
                           prop={'size': 4}, loc='best', fancybox=True,
                           framealpha=0.5)
            # ax.set_yscale('log')
            ax3.set_xticks(ind+adj/2)
            """fac = []
            for sam in (bio_rep["Cellmix1"]).columns:
                # print i
                cols = sam.split("_")
                fac.append(cols[0]+"_"+cols[1]+"_"+cols[2])"""
            # group_labels = conc_list
            ax3.set_xticklabels(group_labels, fontsize=4, rotation=45,
                                ha='right')
            for tick in ax3.yaxis.get_major_ticks():
                tick.label.set_fontsize(4)
            # fig.autofmt_xdate())

            try:
                plt.suptitle(self.antibodydata.ix[int(i), "PI_name"]+"_Ab-ID = " +
                             str(i)+"_Dilutions = " +
                             str(self.antibodydata.ix[int(i), "Dilution"]),
                             fontweight='bold')
                ab_name = self.antibodydata.ix[int(i), "PI_name"]
            except:
                plt.suptitle("Ab-ID = "+str(i), fontweight='bold')
                ab_name = str(i)

            if bars_for_pi == "yes":
                if (ab_name[0:2] == "p-"):
                    calibarator = "yes"
                else:
                    calibarator = "no"
            else:
                calibarator = "yes"

            if calibarator == "yes":

                ax4 = plt.subplot2grid((n_row, n_col), (0, 1))
                # ax = fig.add_subplot(1,1,1)
                if (calhp_ctrls == "both") or (calhp_ctrls == "1"):
                    ind = np.arange(0, 2*len((bio_rep["Cal_1_HP"]). ix[i]), 2)
                    bar1 = ax4.bar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,
                                   width=0.5,
                                   label=(bio_rep["Cal_1_HP"
                                                  ]).columns.values,
                                   color=color1, edgecolor = "none",
                                   align='center', log=0)
                    ax4.errorbar(ind, (bio_rep["Cal_1_HP"]).ix[i].values,
                                 yerr=(bio_rep_std["Cal_1_HP"]).ix[i].values,
                                 fmt='none', ecolor='black', elinewidth=0.5,
                                 capsize=1, capthick=0.5)
                    fac = []
                    for sam in (bio_rep["Cal_1_HP"]).columns:
                        # print i
                        cols = sam.split("_")
                        fac.append(cols[2])
                if (calhp_ctrls == "both") or (calhp_ctrls == "2"):
                    ind = np.arange(0, 2*len((bio_rep["Cal_2_HP"]). ix[i]), 2)
                    bar2 = ax4.bar(ind+adj, (bio_rep["Cal_2_HP"]).ix[i].values,
                                   width=0.5,
                                   label=(bio_rep["Cal_2_HP"]).columns.values,
                                   color=color2, edgecolor = "none",
                                   align='center', log=0)
                    ax4.errorbar(ind+adj, (bio_rep["Cal_2_HP"]).ix[i].values,
                                 yerr=(bio_rep_std["Cal_2_HP"]).ix[i].values,
                                 fmt='none', ecolor='black', elinewidth=0.5,
                                 capsize=1, capthick=0.5)
                    fac = []
                    for sam in (bio_rep["Cal_2_HP"]).columns:
                        # print i
                        cols = sam.split("_")
                        fac.append(cols[2])
                # ax4.set_ylabel('Intensity')
                ax4.set_ylim(0)
                if (calhp_ctrls == "both"):
                    ax4.legend((bar1[0], bar2[0]), ('Cal_1_HP', 'Cal_2_HP'),
                               prop={'size': 4}, loc='best', fancybox=True,
                               framealpha=0.5)
                # ax.set_yscale('log')
                ax4.set_xticks(ind+adj/2)

                # group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax4.set_xticklabels(group_labels, fontsize=4, rotation=45,
                                    ha='right')
                for tick in ax4.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                # fig.autofmt_xdate()

                ax5 = plt.subplot2grid((n_row, n_col), (0, 2))
                # ax = fig.add_subplot(1,1,1)
                if (caljc_ctrls == "both") or (caljc_ctrls == "1"):
                    ind = np.arange(0, 2*len((bio_rep["Cal_1_JC"]).ix[i]), 2)
                    bar1 = ax5.bar(ind, (bio_rep["Cal_1_JC"]).ix[i].values,
                                   width=0.5,
                                   label=(bio_rep["Cal_1_JC"]).columns.values,
                                   color=color1, edgecolor = "none",
                                   align='center', log=0)
                    ax5.errorbar(ind, (bio_rep["Cal_1_JC"]).ix[i].values,
                                 yerr=(bio_rep_std["Cal_1_JC"]).ix[i].values,
                                 fmt='none', ecolor='black', elinewidth=0.5,
                                 capsize=1, capthick=0.5)
                    fac = []
                    for sam in (bio_rep["Cal_1_JC"]).columns:
                        # print i
                        cols = sam.split("_")
                        fac.append(cols[2])

                if (caljc_ctrls == "both") or (caljc_ctrls == "2"):
                    ind = np.arange(0, 2*len((bio_rep["Cal_2_JC"]).ix[i]), 2)
                    bar2 = ax5.bar(ind+adj, (bio_rep["Cal_2_JC"]).ix[i].values,
                                   width=0.5,
                                   label=(bio_rep["Cal_2_JC"]).columns.values,
                                   color=color2, edgecolor = "none",
                                   align='center', log=0)
                    ax5.errorbar(ind+adj, (bio_rep["Cal_2_JC"]).ix[i].values,
                                 yerr=(bio_rep_std["Cal_2_JC"]).ix[i].values,
                                 fmt='none', ecolor='black', elinewidth=0.5,
                                 capsize=1, capthick=0.5)
                    fac = []
                    for sam in (bio_rep["Cal_2_JC"]).columns:
                        # print i
                        cols = sam.split("_")
                        fac.append(cols[2])
                # ax4.set_ylabel('Intensity')
                ax5.set_ylim(0)
                if (caljc_ctrls == "both"):
                    ax5.legend((bar1[0], bar2[0]), ('Cal_1_JC', 'Cal_2_JC'),
                               prop={'size': 4}, loc='best', fancybox=True,
                               framealpha=0.5)
                # ax.set_yscale('log')
                ax5.set_xticks(ind+adj/2)

                # group_labels = (bio_rep["Cellmix1"]).columns.values
                group_labels = fac
                ax5.set_xticklabels(group_labels, fontsize=4, rotation=45,
                                    ha='right')
                for tick in ax5.yaxis.get_major_ticks():
                    tick.label.set_fontsize(4)
                # fig.autofmt_xdate()


                """try:
                    plt.suptitle(database.ix[int(i), "PI_name"]+"_Ab-ID = " +
                                 str(i)+"_Dilutions = " +
                                 str(database.ix[int(i), "Dilution"]),
                                 fontweight='bold')
                except:
                    plt.suptitle("Ab-ID = "+str(i), fontweight='bold')"""

            plt.tight_layout()
            plt.subplots_adjust(top=0.90)
            # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

            """plt.suptitle(data.ix[i,"ant_name"]+"_Ab-ID = "+str(i)+
            "_Dilutions = "+str(database.ix[i,"RPPA0013_Dilution"]),
            fontweight='bold')"""

            # savefig("test.pdf",dpi=200)
            pp.savefig(dpi=200, bbox_inches='tight')
            plt.close()
            # except:
            # print str(i)+" error"

        pp.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculation of signature score for samples')
    parser.add_argument('-sl', '--slide_table', type=str, required=True,
                        help='Slide table')
    parser.add_argument('-a', '--ant_dbase', type=str, required=True,
                        help='Antibody database')
    parser.add_argument('-p1', '--path', type=str, required=True,
                        help='path')
    parser.add_argument('-qp1', '--qpath', type=str, required=False,
                        help='Path for Quincy files')
    parser.add_argument('-p2', '--pmt', nargs='+', type=str, required=True,
                        help='PMT settings')
    parser.add_argument('-st', '--start', type=str, required=True,
                        help='Start of folders')
    parser.add_argument('-qst', '--qstart', type=str, required=False,
                        help='Start of folders for Quincy files')
    parser.add_argument('-e', '--end', type=str, required=True,
                        help='End of files')
    parser.add_argument('-pr', '--prot', type=str, required=True,
                        help='Protein data')
    parser.add_argument('-ex', '--experiment', type=str, required=True,
                        help='Experiment number')
    parser.add_argument('-sa', '--sample', type=str, required=True,
                        help='Sample list')
    parser.add_argument('-co', '--conf', type=str, required=True,
                        help='Configuration file')
    parser.add_argument('-ip', '--imagepath', type=str, required=True,
                        help='Image file')
    parser.add_argument('-d', '--debug', type=str, required=False,
                        help='Image file')
    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    slide_table = args.slide_table
    ant_dbase = args.ant_dbase
    path = args.path
    pmt_settings = args.pmt
    start = args.start
    end = args.end
    experiment = args.experiment
    conf_file = args.conf
    image_path = args.imagepath
    if args.debug:
        debug = args.debug
    else:
        debug = 0

    prot_data = pd.read_table(args.prot, sep="\t", skiprows=32)

    rppa = RPPA(slide_table, ant_dbase)

    norm_data = rppa.normalize(args.path,start,end,pmt_settings,prot_data,debug)
    #qnorm_data = rppa.norm_data(args.qpath, args.qstart, pmt_settings)
    raw_data = rppa.raw_data(args.path, start, end, pmt_settings)
    flag_data = rppa.flag_data(args.path, start, end, pmt_settings)

    for pmt in pmt_settings:
        norm_data[pmt].to_csv(experiment+"_Norm_"+pmt+".xls", sep="\t", na_rep='NA')
        raw_data[pmt].to_csv(experiment+"_Raw_"+pmt+".xls", sep="\t", na_rep='NA')
        flag_data[pmt].to_csv(experiment+"_Flag_"+pmt+".xls", sep="\t", na_rep='NA')

    for pmt in pmt_settings:
        print pmt

    sample_data = pd.read_table(experiment+"_Norm_"+pmt+".xls", sep="\t", index_col=0)
    sample_list = pd.read_table(args.sample, sep="\t", index_col=0)

    sample_table_file = rppa.sample_table(sample_data, sample_list, experiment)

    list1 = [["Ctrl_IgGmix", "0.", "blank.", "Ctrl_GridHP"], []]

    for pmt in pmt_settings:
        rppa.pi_data(experiment+"_Norm_"+pmt+".xls", sample_table_file,
                     experiment, list1, conf_file)
        rppa.pi_data(experiment+"_Raw_"+pmt+".xls", sample_table_file,
                     experiment, list1, conf_file)
        rppa.pi_data(experiment+"_Flag_"+pmt+".xls", sample_table_file,
                     experiment, list1, conf_file)

    pi_list = []
    mouse_pi_list = []
    # pi_list = ["YL"]
    # mouse_pi_list = [""]
    with open(conf_file) as f:
        while True:
            text1 = f.readline()
            if text1 == "":
                break
            cols = text1.strip().split('\t')
            if cols[0] == "pi_list":
                pi_list_len = len(cols)
                for i in range(1, pi_list_len):
                    pi_list.append(cols[i])
            if cols[0] == "mouse_pi_list":
                mouse_pi_list_len = len(cols)
                for i in range(1, mouse_pi_list_len):
                    mouse_pi_list.append(cols[i])

    rppa.pi_directory(experiment, pmt_settings, pi_list, mouse_pi_list)

    rppa.reports(experiment, pi_list, mouse_pi_list, "")

    # return imagedata, qc_mean_data, qc_sd_data, bio_rep, bio_rep_std
    imagedata = rppa.imagedf(experiment, image_path)

    """for pmt in pmt_settings:
        bar_data = rppa.ctrl_bardata(experiment,
                                     experiment+"_Norm_"+pmt+".xls",
                                     conf_file, list1)
        rppa.combined_plots(bar_data[0], bar_data[1], bar_data[2], bar_data[3],
                            pmt,
                            experiment+"_Norm_"+pmt+"_Barcharts_with_max.pdf",
                            imagedata=imagedata, bars_for_pi=None)
        rppa.combined_plots((bar_data[0]).drop("Slide_Max", axis=1),
                            (bar_data[1]).drop("Slide_Max", axis=1),
                            bar_data[2], bar_data[3],
                            pmt,
                            experiment+"_Norm_"+pmt+"_Barcharts_without_max.pdf",
                            imagedata=imagedata, bars_for_pi=None)

    for pmt in pmt_settings:
        bar_data = rppa.ctrl_bardata(experiment,
                                     experiment+"_Raw_"+pmt+".xls",
                                     conf_file, list1)
        rppa.combined_plots(bar_data[0], bar_data[1], bar_data[2], bar_data[3],
                            pmt,
                            experiment+"_Raw_"+pmt+"_Barcharts_with_max.pdf",
                            imagedata=imagedata, bars_for_pi=None)
        rppa.combined_plots((bar_data[0]).drop("Slide_Max", axis=1),
                            (bar_data[1]).drop("Slide_Max", axis=1),
                            bar_data[2], bar_data[3],
                            pmt,
                            experiment+"_Raw_"+pmt+"_Barcharts_without_max.pdf",
                            imagedata=imagedata, bars_for_pi=None)"""

    master_rawdata_dic = {}
    master_normdata_dic = {}
    for pmt in pmt_settings:
        data_org = pd.read_table(experiment+"_Raw_"+pmt+".xls",
                                 sep="\t", index_col=0)
        master_rawdata_dic[pmt] = data_org.ix[:, 1:]

        data_org = pd.read_table(experiment+"_Norm_"+pmt+".xls",
                                 sep="\t", index_col=0)
        master_normdata_dic[pmt] = data_org.ix[:, 1:]
        # print (master_normdata_dic[pmt]).iloc[0:1, :]

    pi_bar_data = rppa.ctrl_bardata(experiment,
                                    experiment+"_Norm_"+pmt+".xls",
                                    conf_file, list1, pi="yes",
                                    master_data_dic=master_normdata_dic,
                                    datatype="Norm")
    # print (pi_bar_data[0]).columns
    # print (pi_bar_data[1]).columns
    # print (pi_bar_data[2]).keys()
    # print (pi_bar_data[3]).keys()
    rppa.combined_plots(pi_bar_data[0], pi_bar_data[1],
                        pi_bar_data[2], pi_bar_data[3],
                        pmt, experiment+"_Norm_PI_barcharts_with_max.pdf",
                        imagedata=None, bars_for_pi="yes")
    rppa.combined_plots((pi_bar_data[0]).drop("Slide_Max", axis=1),
                        (pi_bar_data[1]).drop("Slide_Max", axis=1),
                        pi_bar_data[2], pi_bar_data[3],
                        pmt, experiment+"_Norm_PI_barcharts_without_max.pdf",
                        imagedata=None, bars_for_pi="yes")

    pi_bar_data = rppa.ctrl_bardata(experiment,
                                    experiment+"_Raw_"+pmt+".xls",
                                    conf_file, list1, pi="yes",
                                    master_data_dic=master_rawdata_dic,
                                    datatype="Raw")

    rppa.combined_plots(pi_bar_data[0], pi_bar_data[1],
                        pi_bar_data[2], pi_bar_data[3],
                        pmt, experiment+"_Raw_PI_barcharts_with_max.pdf",
                        imagedata=None, bars_for_pi="yes")
    rppa.combined_plots((pi_bar_data[0]).drop("Slide_Max", axis=1),
                        (pi_bar_data[1]).drop("Slide_Max", axis=1),
                        pi_bar_data[2], pi_bar_data[3],
                        pmt, experiment+"_Raw_PI_barcharts_without_max.pdf",
                        imagedata=None, bars_for_pi="yes")

    for pmt in pmt_settings:
        rppa.cv(experiment, experiment+"_Raw_"+pmt+".xls", pmt)
        rppa.cv(experiment, experiment+"_Norm_"+pmt+".xls", pmt)

if __name__ == '__main__':
    try:
        main()
    except:
        print "An unknown error occurred.\n"
        raise