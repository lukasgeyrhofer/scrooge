# ************************************************************************* #
# **         SCROOGE                                                     ** #
# **         estimate genome size from single copy gene coverage         ** #
# **                                                                     ** #
# **                                                                     ** #
# **                                                                     ** #
# **                                                                     ** #
# **                                                                     ** #
# ************************************************************************* #


import xml.etree.ElementTree as ET
import os,sys
import subprocess
from math import sqrt



# ***************************************************************** #
# **         helper routines                                     ** #
# ***************************************************************** #


def file_exists(filenames):
  r = True
  for f in filenames:
    if os.path.isfile(str(f).strip()):
      r = False
  return r

# ***************************************************************** #
# **         class to load parameters from external XML file     ** #
# **         and execute those programs (after checks)           ** #
# ***************************************************************** #
class externalprogram:
    def __init__(self,optionfile,step,verbose=True):
        try:
            self.__options = ET.parse(optionfile)
        except:
            raise IOError
        self.__xmlroot = self.__options.getroot()
        self.__haveexecutable = False
        self.__executable = ""
        self.__parameters = []
        self.__kwargs = {}
        self.__flags = []
        self.__fnamestdout = None
        self.__fnamestderr = None
        self.__verbose = verbose
        self.__filetypes = {}

        self.__hit = self.__xmlroot.findall("./*[@step='"+step+"']")
        if len(self.__hit) >= 1:
            for t in self.__hit[0].getchildren():
                if t.tag == "executable":
                    self.__executable = t.text
                    self.__haveexecutable = True
                elif t.tag == "flag":                          # usage in cmdline "--flag"        (prepend --)
		    if t.attrib.has_key('name'):
			self.__flags.append(t.attrib['name'])
                elif t.tag == "option":                        # usage in cmdline "-option value" (prepend -)
		    if t.attrib.has_key('name'):
			self.__kwargs[t.attrib['name']] = t.text
                elif t.tag == "parameter":                     # usage in cmdline "paramter"      (no prepend)
                    self.__parameters.append(t.text)
                elif t.tag == "filetype":
		    if t.attrib.has_key('name'):
			if not self.__filetypes.has_key(t.attrib['name']):
			    self.__filetypes[t.attrib['name']] = []
			self.__filetypes[t.attrib['name']].append(t.text)
        else:
            raise ValueError
        

    def __str__(self):
        if self.__haveexecutable:
            return " ".join(self.cmdlineparameters())
        else:
            return None
    
    def get_filetypes(self,filetype):
	if self.__filetypes.has_key(filetype):
	    return self.__filetypes[filetype]
	else:
	    return None
    
    def get_files(self,option):
	if self.__kwargs.has_key(option):
	    if self.get_filetypes(option) != None:
		return [self.get_option(option) + ft for ft in self.get_filetypes(option)]
	    else:
		return [self.get_option(option)]
	else:
	    return None
    
    def check_existence(self):
        if self.__haveexecutable:
            def is_exe(fpath):
                return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
            fpath, fname = os.path.split(self.__executable)
            if fpath:
                if is_exe(self.__executable):
                    return True
            else: 
                for path in os.environ["PATH"].split(os.pathsep):
                    path = path.strip('"')
                    path = os.path.expanduser(path)
                    exe_file = os.path.join(path, self.__executable)
                    if is_exe(exe_file):
                        return True
            return False    
        else:
            return False
    
    def execute(self,wait=True):
        if self.check_existence():
            if self.__fnamestdout == None:
                self.__pso = sys.stdout
            else:
                self.__pso = open(self.__fnamestdout,'w')
            if self.__fnamestderr == None:
                self.__pse = sys.stderr
            else:
                self.__pse = open(self.__fnamestderr,'w')
            if self.__verbose: print "EXECUTE: ",
            if self.__verbose: print self
            self.__pid = subprocess.Popen(self.cmdlineparameters(),stdout = self.__pso,stderr = self.__pse)
            if wait:
                self.__pid.wait()
            return self.__pid
        else:
	    if self.__verbose: print "EXECUTABLE NOT FOUND"
            return None

    def __del__(self):
        try:
            self.__pse.close()
        except:
            pass
        try:
            self.__pso.close()
        except:
            pass
        #if self.__fnamestderr != None:
            #self.__pse.close()
        #if self.__fnamestdout != None:
            #self.__pso.close()

    def cmdlineparameters(self):
        tmp = [self.__executable]
        for a in self.__kwargs.iterkeys():
            tmp.append("-%s"%a)
            if self.__kwargs[a]:
		tmp.append(self.__kwargs[a])
        for a in self.__flags:
            tmp.append("--%s"%a)
        for a in self.__parameters:
            tmp.append(a)
        return tmp

    def set_option(self,option,value,outfile=False):
        self.__kwargs[option] = value
        
    def add_flag(self,flag):
        self.__flags.append(flag)
        
    def add_executable(self,executable):
        self.__executable = executable
        
    def del_option(self,option):
        try:
            del self.__kwargs[option]
        except:
            pass
    
    def del_flag(self,flag):
        try:
            self.__flags.remove(flag)
        except:
            pass
                
    def get_option(self,option):
        if self.__kwargs.has_key(option):
            return self.__kwargs[option]
        else:
            return None
    def get_executable(self):
        if self.__haveexecutable:
            return self.__executable
        else:
            return None

    def set_stderr(self,fname):
        self.__fnamestderr = fname
    def set_stdout(self,fname):
        self.__fnamestdout = fname

    def add_parameter(self,param):
        self.__parameters.append(param)
        
    def del_parameter(self,param):
        try:
            self.__parameters.remove(param)
        except:
            pass
    
    def get_parameters(self):
        return self.__parameters
    def get_flags(self):
        return self.__flags





# ***************************************************************** #
# **         single copy gene sequence                           ** #
# **         with information on coverage                        ** #
# ***************************************************************** #
class SequenceRecord():
    
    def __init__(self,sid,sequence):
        self.__name = sid
        self.__sequence = sequence
        self.__length = len(sequence)
        self.__coverage = [0 for i in range(self.__length)]
        self.__countreads = 0
    
    def __len__(self):
        return self.__length
    
    def __str__(self):
        return self.__name
    
    def get_sequence(self):
        return self.__sequence
    
    def add_coverage(self,mi,ma):
        for i in range(mi,ma):
            self.__coverage[i] += 1
        self.__countreads += 1

    def get_count_reads(self):
	return self.__countreads

    def get_coverage(self):
	return self.__coverage

    def get_coverage_mean(self):
        if self.__countreads > 0:
            return float(sum(self.__coverage))/float(self.__length)
        else:
            return None
    
    def get_coverage_stddev(self):
        if self.__countreads > 0:
            c2 = float(sum([ci**2 for ci in self.__coverage]))
            c1 = float(sum(self.__coverage))
            n = float(self.__length)
            return sqrt(n*c2-c1*c1)/sqrt(n*n-n)
        else:
            return None



# ***************************************************************** #
# **         list of all single copy genes                       ** #
# ***************************************************************** #
class SingleCopyGeneList():
    def __init__(self,scglength = None, readlength = None):
        self.__seqid = []
        self.__seq = {}
        self.__readminlenght = readlength
        self.__scgminlength = scglength
    
    def __getitem__(self,sid):
        return self.__seq[sid]
    
    def __iter__(self):
        for sid in self.__seqid:
            yield self.__seq[sid]
    
    def add_sequence(self,sid,sequence):
	if self.__scgminlength:
	    if self.__scgminlength > len(str(sequence)):
		return None
        self.__seqid.append(sid)
        self.__seq[sid] = SequenceRecord(sid,str(sequence))
        
    def add_coverage(self,sid,start=None,end=None):
        if sid in self.__seqid:
            mi = min(start,end)
            if mi<0:mi=0
            ma = max(start,end)
            if ma > len(self.__seq[sid]):ma = len(self.__seq[sid])
            if self.__readminlenght:
		if ma-mi > self.__readminlenght:
		    self.__seq[sid].add_coverage(mi,ma)
	    else:
		self.__seq[sid].add_coverage(mi,ma)
        else:
            print >> sys.stderr,"did not find sID"

    def get_ids(self):
        return self.__seqid

    def write_sequence_file(self,filename):
        """" write in FASTA format """
        
        f = open(filename,"w")
        for sid in self.__seqid:
            f.write(">"+sid+"\n")
            f.write(self.__seq[sid].get_sequence()+"\n")
        f.close()
        del f
    
    def write_coverage_file(self,filename):
	f = open(filename,"w")
	for sid in self.__seqid:
	    print >> f,"#",sid,self.__seq[sid].get_count_reads()
	    c = self.__seq[sid].get_coverage()
	    s = self.__seq[sid].get_sequence()
	    for i in range(len(self.__seq[sid])):
		print >>f,i,c[i],s[i]
	    print >> f
	f.close()
	


