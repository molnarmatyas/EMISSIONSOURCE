#!/bin/bash
cp urqmd-3.4/drho_analyze/converter_f19.cc ./converter_f19_backup.cc
cp urqmd-3.4/drho_analyze/pairsource_urqmd.cc ./pairsource_urqmd_backup.cc
git add ./converter_f19_backup.cc ./pairsource_urqmd_backup.cc
#git commit -m "Backup of modified UrQMD to Drho source files"
#git push origin master
