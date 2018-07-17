BASE=$PWD  ###This should be a temporary directory
CONFIG=config_Ath_TAIR10.sh

{
$CONFIG
MSG="==== Inspect and make sure all [TEST] worked"
prompt_yn "$MSG" && . $CONFIG
}

{
echo "[index]"
find /home/feng/syno3/PW_HiSeq_data/RNA-seq/Raw_data/187R-2/S* -type d | tee ath.index \
  | head 
echo "[command]"
./mapper.template ath.index | head
prompt_yn "#### Inspect the command and make sure they make sense"
}

{
   CMD="parallel --gnu -j5 <ath.sh 2>&1 | tee ath.log &"
   MSG="====You are now ready to start your pipeline with: \n\
   $CMD & \n\
   "
   printf "$MSG"
   prompt_yn "" && eval $CMD
   printf "====pipeline started, monitor by : \n\
   watch tree -h -L 2 . 
"
}
