from CRABClient.UserUtilities import config, ClientException, getUsernameFromCRIC
#from input_crab_data import dataset_files
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser

production_tag = datetime.date.today().strftime('%Y%b%d')

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
#config.General.workArea = 'B0ToXKsLivia_%s' % production_tag
#config.General.workArea = 'BParkingUL2018_%s' % production_tag
config.General.workArea = 'BParkingUL2018_test'
#config.General.workArea = 'CharmoniumUL2018_TrgTest'#_%s' % production_tag

config.section_('Data')
config.Data.publication = False
# To save on /eos at CERN !!!!!!phys_bphys dir full!!!!!!
#####config.Data.outLFNDirBase = '/store/group/phys_bphys/crovelli/nanoaod_X/%s' % (config.General.workArea)
# To save on Rome T2
#config.Data.outLFNDirBase = '/store/user/cquarant/CharmoniumUL2018/%s' % (config.General.workArea)
config.Data.outLFNDirBase = '/store/user/cquarant/BParkUL2018/%s' % (config.General.workArea)
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../test/run_nano_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_('User')
config.section_('Site')
# chiara: scommenta x salvare al cern
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T2_IT_Rome'

if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from httplib import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config = config, dryrun=False)
      except HTTPException as hte:
          print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)


  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'samples.yml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}
    
    # loop over samples
    for sample, info in doc['samples'].iteritems():
      # Given we have repeated datasets check for different parts
      parts = info['parts'] if 'parts' in info else [None]
      for part in parts:
        name = sample % part if part is not None else sample
        
        # filter names according to what we need
        if not fnmatch(name, args.filter): continue
        print 'submitting', name

        isMC = info['isMC']

        config.Data.inputDataset = info['dataset'] % part \
                                   if part is not None else \
                                      info['dataset']

        config.General.requestName = name
        common_branch = 'mc' if isMC else 'data'
        config.Data.splitting = 'FileBased' if isMC else 'LumiBased'
        if not isMC:
            config.Data.lumiMask = info.get(
                'lumimask', 
                common[common_branch].get('lumimask', None)
            )
        else:
            config.Data.lumiMask = ''

        config.Data.unitsPerJob = info.get(
            'splitting',
            common[common_branch].get('splitting', None)
        )
        globaltag = info.get(
            'globaltag',
            common[common_branch].get('globaltag', None)
        )
        
        config.JobType.pyCfgParams = [
            'isMC=%s' % isMC, 'reportEvery=1000',
            'tag=%s' % production_tag,
            'globalTag=%s' % globaltag,
        ]
        
        config.JobType.outputFiles = ['_'.join(['xNANO', 'mc' if isMC else 'data', production_tag])+'.root']
        
        print config
        submit(config)

