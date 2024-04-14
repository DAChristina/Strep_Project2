for file in *.fasta; do
  name=${file%.fasta} 
  sed -i "s/>/>$name:/" "$file"
  
done

cat *.fasta > AllTarget_db.fasta

# Previously I create 2 folders to separate GPSC2 (n=1) & GPSC31 (n=74).
# This script is created to summarise all 74 targeted fasta files for SKA2 feed.
