#PBS -lselect=1:ncpus=100:mem=1200gb
#PBS -lwalltime=24:00:00
#PBS -o /rdsgpfs/general/user/aalmeid2/projects/cardiac_single_cell_biology/live/Pdgfra_project/cellchat.out

cd $Home

source activate r_403

Rscript --vanilla cellchat_script_control_rv.R