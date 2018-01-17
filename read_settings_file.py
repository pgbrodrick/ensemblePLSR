




import ConfigParser
import sys,os,re
import numpy as np


class settings():
  
  def clean_comments(self,instring):
     orig_instring = instring
     while ('(' in instring):
       start_ind = instring.index('(')
       if (')' not in instring):
         print 'unmatched ( in input::::         ',orig_instring
         print 'terminating'
         quit()
       end_ind = instring.split('=')[0].rindex(')')
       instring = instring.replace(instring[start_ind:end_ind+1],'')
     return instring
  
  
  def get_all_settings(self,set_obj):
      out = {}
      for section in set_obj.sections():
        options = set_obj.options(section)
        for option in options:
          read = set_obj.get(section,option)
          option = self.clean_comments(option).lower().strip()
          read = self.clean_comments(read)

          if (',' in read):
            read = re.split(',',read)

          if (read == ''):
            read = None
            
          
          out[option] = read
      return out
  
  def get_setting(self,ref):
    if (ref not in self.settings_dict):
      print 'could not locate setting\'' + ref + '\' in  settings file',self.settings_filename
      print 'terminating'
      quit()
    return self.settings_dict[ref]


  def clean_and_check(self):

    if ('csv file' not in self.settings_dict):
      print 'no csv input provided.  terminating'
      quit()
    if (self.settings_dict['csv file'] == None):
      print 'no csv input provided.  terminating'
      quit()
    if (os.path.isfile(self.settings_dict['csv file']) == False):
      print 'csv file: \'',self.settings_dict['csv file'],'\' is invalid.  terminating'
      quit()

    
    ##### bad bands list ( list converted to ints)
    if (self.settings_dict['bad bands'] == None):
      self.settings_dict['bad bands'] = np.array([-1])
    else:
      print 'n bad bands',len(self.settings_dict['bad bands'])
      self.settings_dict['bad bands'] = np.array([int(x) for x in self.settings_dict['bad bands']])

    ###### string lists
    for lk in ['chems','chem transforms','ignore columns']:
      if (self.settings_dict[lk] == None):
        self.settings_dict[lk] = [-1]
      else:
        if (not isinstance(self.settings_dict[lk],str)):
          self.settings_dict[lk] = [x.replace('\'','').replace('\"','').strip() for x in self.settings_dict[lk]]
        else:
          self.settings_dict[lk] = [self.settings_dict[lk].strip()]


    bool_vars = ['brightness normalize','find bad bands','scale response',\
                 'scale features']
    for bool_var in bool_vars:
      if (self.settings_dict[bool_var] == 'false' or self.settings_dict[bool_var] == 'False'):
        self.settings_dict[bool_var] = False
      if (self.settings_dict[bool_var] == 'true' or self.settings_dict[bool_var] == 'True'):
        self.settings_dict[bool_var] = True

    int_vars = ['iterations','samples per crown','max band',\
                'min pixel per crown','max components','test set value']
    for int_var in int_vars:
     if (self.settings_dict[int_var] != None):
      self.settings_dict[int_var] = int(self.settings_dict[int_var].strip())

    float_vars = ['ndvi minimum','ndvi maximum','brightness minimum','brightness maximum','test set holdout fraction',\
                  'iteration holdout fraction','iteration fraction used', 'lowest wavelength','wavelength interval']
    for float_var in float_vars:
     if (self.settings_dict[float_var] != None):
      self.settings_dict[float_var] = float(self.settings_dict[float_var])

    string_vars = ['band preface','crown col','version name']
    for string_var in string_vars:
     if (self.settings_dict[string_var] != None):
      self.settings_dict[string_var] = self.settings_dict[string_var].replace('\'','')
      self.settings_dict[string_var] = self.settings_dict[string_var].replace('\"','')
    
    if (self.settings_dict['iteration holdout fraction'] in [None,0] and \
        self.settings_dict['iteration fraction used'] != 1):
       print 'if iteration fraction used is specified, an iteration holdout fraction greater than 0 must be specified as well'
       quit()
  
  

  def __init__(self,filename):
    self.settings_filename = ''
    self.settings_obj = ConfigParser.ConfigParser()
    self.settings_dict = ''
    self.settings_filename = filename
    self.settings_obj.read(filename)
    self.settings_dict = self.get_all_settings(self.settings_obj)










