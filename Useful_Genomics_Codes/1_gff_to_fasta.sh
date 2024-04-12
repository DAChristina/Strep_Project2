for file in *.gff; do
        seqret -sequence "$file" -outseq "${file%.*}.fasta" -osformat fasta
done

# SOURCE: http://emboss.open-bio.org/html/use/ch06s01.html
# SOURCE: https://askubuntu.com/questions/1342080/convert-multiple-files-in-one-command
