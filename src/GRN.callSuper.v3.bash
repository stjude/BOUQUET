module load R/4.0.2
file=${1:-'network.plusIN.promoter.18858.PE_nodes.03172010_614ATAAXX_B5_mm9.sorted.bam.txt'}
ScriptDir=${2:-'script'}
#genes=${3:-'Irf2bpl|NM_145836.2'}
genes=${3:-''}
cols=${4:-'15,16,23,24,25,26,31,32,42,43,49'}

# Call Super based on total enhancer signal (either MED1 or H3K27ac)

a=($(head -1 $file))
cols=(${cols//,/ })

echo -e "#### Selected Measurment For Enhancer Signal Loading ####"
for i in  ${cols[@]} 
do
    echo ${a[i-1]}
done
echo -e "Hilighted genes"
echo $genes

for i in  ${cols[@]} 
do 
    #Arrays in Bash are indexed from zero, and in zsh (mac OS default) they're indexed from one.
    cut -f 1,$i $file|tail -n +2 > $file.${a[$i-1]}.txt
    R --save $file.${a[$i-1]} $file.${a[$i-1]}.txt  $genes < $ScriptDir/PriMROSE_callSuper.flexibleID.R
done


