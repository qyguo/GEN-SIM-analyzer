# How To Use

```bash
cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src
cmsenv
git clone git@github.com:qyguo/GEN-SIM-analyzer.git
scramv1 b -j 8
cd GEN-SIM-analyzer/GenAnalyzer
# Befor running, check the config file and comment/uncomment line for MINIAOD/LHE-GEN
cmsRun python/ConfFile_cfg.py
```
