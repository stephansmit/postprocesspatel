import re
import os, os.path
import math
class PipeCase():
    def __init__(self, caselocation, casemap):
        self.caselocation = caselocation
        self.casename = self.get_casename()
        self.rpow = self.get_rpow(casemap)
        self.mpow = self.get_mpow(casemap)
        self.lpow = self.get_lpow(casemap)
        self.mainfilename = 'main_comp.f'
        self.postprocesslocation = '/vonkarman/apatel/Re395_f_cons_ys'
        self.postprocessfoldername = 'post'
        self.datafoldername = 'DNS'
        self.mainfile  = "/".join([self.caselocation,self.mainfilename])
        self.mainfilelines = self.get_mainfile_lines()
        self.stretch_factor = self.get_stretch_factor()
        self.volumetric_heatsource = self.get_volumetric_heatsource()
        self.streamwise_length = self.get_streamwise_length()
        self.spanwise_length = self.get_spanwise_length()
        self.reynolds = self.get_reynolds()
        self.prandtl = self.get_prandtl()
        self.numberfiles = self.get_number_files()
        self.startnumber = self.get_startnumber()
        #self.copyoldpostprocessfolder()
        #self.copypostprocessfolder()
        self.write_postprocessparamfile()
        self.numberofcores = self.get_number_of_cores()
        self.write_slurmjob()
	#self.compile_postprocess()
        self.submit_slurmjob()

        
    def get_number_files(self):
        data_dir = "/".join([self.caselocation, self.datafoldername])
        files = [name for name in os.listdir(data_dir) if os.path.isfile(data_dir + '/'+name)]
        
        numbers = []
        for i in files:
            regex = "[0-9]+"
            matches = re.search(regex, i)
            numbers.append(matches.group(0))
        startnumber = sorted([ int(x) for x in numbers])[0]
        endnumber = sorted([ int(x) for x in numbers])[-1]
        return endnumber-startnumber
     
    def get_startnumber(self):
        data_dir = "/".join([self.caselocation, self.datafoldername])
        files = [name for name in os.listdir(data_dir) if os.path.isfile(data_dir + '/'+name)]
        numbers = []
        for i in files:
            regex = "[0-9]+"
            matches = re.search(regex, i)
            numbers.append(matches.group(0))
        startnumber = sorted([ int(x) for x in numbers])[0]
        return startnumber
        

    def get_mainfile_lines(self):
        with open (self.mainfile, "r") as myfile:
            data=myfile.readlines()
        return data
    
    
    def get_volumetric_heatsource(self):
        try:
            volumetric = [x for x in self.mainfilelines if "dfc_n =" in x]
            string = volumetric[0]
        except:
            volumetric = [x for x in self.mainfilelines if "dcdt =" in x]
            string = volumetric[0]
            
        regex = "[0-9]+\.[0-9]+"
        matches = re.search(regex, string)
        volumetric_heatsource = matches.group(0)
        return volumetric_heatsource
    
    def get_stretch_factor(self):
        stretch = [x for x in self.mainfilelines if "dx = 0.5" in x]
        regex = "-[0-9]+\.[0-9]+\*\("
        matches = re.search(regex, stretch[0])
        stretch_factor = re.search('[0-9]+\.[0-9]+',matches.group(0)).group(0)
        return stretch_factor
    
    def get_streamwise_length(self):
        string = [x for x in self.mainfilelines if "Lz =" in x][0]
        streamwise_length = eval(
            string.strip('\n')
            .replace('Lz =','')
            .replace(' ','')
            .replace(' ','')
            .replace('*', "*math.")
        )/math.atan(1)
        return streamwise_length
    def get_reynolds(self):
        string = [x for x in self.mainfilelines if "Re =" in x][0]
        reynolds = eval(
            string.strip('\n')
            .replace('Re =','')
            .replace(' ','')
            .replace(' ','').split('!')[0]
        )
        return reynolds
    def get_prandtl(self):
        string = [x for x in self.mainfilelines if "Pr =" in x][0]
        reynolds = eval(
            string.strip('\n')
            .replace('Pr =','')
            .replace(' ','')
            .replace(' ','').split('!')[0]
        )
        return reynolds
    
    
    def get_spanwise_length(self):
        try:
            string = [x for x in self.mainfilelines if "Lt = " in x][0]
            spanwise_length = eval(
                string.strip('\n')
                .replace('Lt =','')
                .replace(' ','')
                .replace(' ','')
                .replace('*', "*math.")
            )/math.atan(1)
        except:
            string = [x for x in self.mainfilelines if "dtheta =" in x][0]
            regex = "[0-9]+\.[0-9]*"
            matches = re.search(regex, string)
            spanwise_length = re.search('[0-9]+\.[0-9]*',matches.group(0)).group(0)

        return spanwise_length

    def get_casename(self):
        return self.caselocation.split('/')[-1]
    def get_rpow(self,casemap):
        return casemap[self.casename]['rpow']
    def get_mpow(self,casemap):
        return casemap[self.casename]['mpow']
    def get_lpow(self,casemap):
        return casemap[self.casename]['lpow']
    
    def copyoldpostprocessfolder(self):
        print("copying the old post to post_old for case:"+ self.casename)


        string = "/".join([self.caselocation,self.postprocessfoldername])
        if os.path.exists(string + "_old"):       
            os.system('cp -r ' + string + ' ' + string+"_old2")
            #print('cp -r ' + string + ' ' + string+"_old2")
        else:
            os.system('cp -r ' + string + ' ' + string+"_old")
            #print('cp -r ' + string + ' ' + string+"_old")

    def copypostprocessfolder(self):
        print("copying the new post to post for case:"+ self.casename)
        os.system('cp -r ' + "/".join([self.postprocesslocation,self.postprocessfoldername]) + ' ' + self.caselocation + '/')
        #print('cp -r ' + "/".join([self.postprocesslocation,self.postprocessfoldername]) + ' ' + self.caselocation + '/')
    
    
    def write_postprocessparamfile(self):
        filename = "/".join([self.caselocation ,self.postprocessfoldername, 'parameters'])
        print("writing the parameter file for case:"+ self.casename)
        file = open(filename,"w") 
        file.write(" ".join([str(self.numberfiles),str(self.startnumber)])+'\n')
        file.write(" ".join([str(self.reynolds),str(self.prandtl)])+'\n')
        file.write(" ".join([str(self.streamwise_length),str(self.spanwise_length),str(self.stretch_factor)])+'\n')
        file.write(" ".join([str(self.rpow), str(self.mpow), str(self.lpow)])+'\n')
        file.write(str(self.volumetric_heatsource))
        file.close()
    def get_number_of_cores(self):
        return 24
    
    def write_slurmjob(self):
        print("writing job file")
        # Read in the template file
        with open('JOB_template', 'r') as file :
            filedata = file.read()
        # Replace the target string
        filedata2 = filedata.replace('REPLACE_CORES', str(self.numberofcores))
        self.jobname = "JOB_"+self.casename
        self.jobfile = "/".join([self.caselocation ,self.postprocessfoldername, self.jobname])
        # Write the file out again
        with open(self.jobfile, 'w') as file:
            file.write(filedata2)
        return
    
    def submit_slurmjob(self):
        print("submitting job for case:" + self.casename)
        #print('cd '+ '/'.join([self.caselocation, self.postprocessfoldername]) + ' && sbatch '+ self.jobname)
        os.system('cd '+ self.caselocation + ' && sbatch '+ self.jobname)


    def compile_postprocess(self):
        print('compiling post processing for case:' +self.casename)
        #print('cd '+ '/'.join([self.caselocation, self.postprocessfoldername]) + ' && make -f makepost')
        os.system('cd '+ '/'.join([self.caselocation, self.postprocessfoldername]) + ' && make -f makepost')

    




if __name__=="__main__":
    cases = [
    '/tennekes/apatel/Re395_f_cons',
    '/vonkarman/apatel/Re395_f_cons_ys',
    '/vonkarman/apatel/Re395_f_gas',
    '/vonkarman/apatel/Re395_f_gl',
    #'/vonkarman/apatel/Re395_f_liq',
    '/vonkarman/apatel/Re395_f_cons_nu2',
    '/vonkarman/apatel/Re395_f_varmusqrt2',
    '/vonkarman/apatel/sca_gl',
    '/vonkarman/apatel/sca_vl']
    
    casemap = {
        'Re395_f_cons':       {'rpow': 0, 'mpow': 0,'lpow': 0},
        'Re395_f_cons_ys':    {'rpow': -1, 'mpow': -0.5,'lpow': 0},
        'Re395_f_gas':        {'rpow': -1,'mpow': 0.7,'lpow': 0},
        'Re395_f_gl':         {'rpow': 0, 'mpow': 1.2,'lpow': 0},
        'Re395_f_liq':        {'rpow': 0, 'mpow': -1,'lpow': 0},
        'Re395_f_cons_nu2':   {'rpow': -1, 'mpow': -1,'lpow': 0},
        'Re395_f_varmusqrt2': {'rpow': 0, 'mpow': -0.5 ,'lpow': 0},
        'sca_gl':             {'rpow': -1, 'mpow': 0.7,'lpow': 0.7},
        'sca_vl':             {'rpow': 0, 'mpow': 0,'lpow': 1}
                  }
    Cases = [PipeCase(x, casemap) for x in cases]


