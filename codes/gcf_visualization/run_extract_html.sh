cd /home/mlk442/FIND/MicroMGx_Pilot/genomics/antismash_outputs
for folder in DB*
do 
echo $folder
python3 /home/mlk442/FIND/mlk442/BGC_table/extract_html.py -i /home/mlk442/FIND/MicroMGx_Pilot/genomics/antismash_outputs -o /home/mlk442/FIND/mlk442/BGC_table -s $folder -g /home/mlk442/FIND/mlk442/gcfs_data.txt
done
