### MASTER

#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")


sbatch elr_fitQTL.sh
sbatch brp_fitQTL.sh
sbatch new_fitQTL.sh
sbatch nbh_fitQTL.sh


sbatch fitQTL_long.sh NBH
sbatch fitQTL_long.sh BRP
sbatch fitQTL_long.sh NEW
sbatch fitQTL_long.sh ELR
