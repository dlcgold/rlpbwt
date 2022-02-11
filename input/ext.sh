while IFS=$'\t' read -r -a myArray
do
 echo "${myArray[4]}"
done < test.txt
