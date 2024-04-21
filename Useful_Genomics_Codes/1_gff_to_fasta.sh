for file in *.gff; do
  name=${file%.gff} 
  grep -A999999 "##FASTA" $file | grep -v "##FASTA" >> ./fa/${name}.fasta
done


# For managing gff in some folders & creating many FASTA files from a gff, see:
# https://www.microbialsystems.cn/en/post/gff3tofasta/

# Seqret by EMBOSS:
# SOURCE: http://emboss.open-bio.org/html/use/ch06s01.html
# SOURCE: https://askubuntu.com/questions/1342080/convert-multiple-files-in-one-command

# Some notes:
# gff files (General Feature Format) has already contains FASTA with annotations.
# The only thing that is reasonable is to extract the FASTA-part of the gff files and collect it in a single FASTA file output.
